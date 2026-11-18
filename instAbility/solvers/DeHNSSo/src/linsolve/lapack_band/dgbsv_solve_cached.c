/* License (GNU GENERAL PUBLIC LICENSE v3)
 *                  Delft Harmonic Navier-Stokes Solver
 *     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *     Please cite this code and the paper if you have used this for your
 *     publication:
 *
 * DOI:
 */

/* dgbsv_solve_cached.c
 *
 * Banded LU (LU factorization) SOLVE step (LAPACK *gbtrs), using a factor
 * already produced by dgbsv_factor. Reuses the cached band matrix across
 * right-hand sides (RHS) and Picard iterations, with no re-factorization.
 *
 *   x = dgbsv_solve_cached(ABf, IPIV, KL, KU, rhs)
 *
 * ABf : [LDAB x N] double matrix from dgbsv_factor (real or complex); LDAB
 *       (leading dimension of the banded array) = 2*KL+KU+1, with KL and KU
 *       the lower and upper bandwidth.
 * IPIV: [N x 1] int64 pivot indices from dgbsv_factor
 * rhs : dense N-by-1 double vector (real or complex)
 *
 * If ABf is real but rhs is complex (the MFD (0,0) case), the real factor is
 * applied to the real and imaginary parts separately.
 *
 * Compile via src/linsolve/lapack_band/mexbuild_dgbsv.m (-lmwlapack).
 */

#include "mex.h"
#include "lapack.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 5)
        mexErrMsgIdAndTxt("dgbsv_solve_cached:nrhs",
                          "Usage: x = dgbsv_solve_cached(ABf, IPIV, KL, KU, rhs)");
    const mxArray *ABm = prhs[0];
    const mxArray *IPm = prhs[1];
    ptrdiff_t KL = (ptrdiff_t)mxGetScalar(prhs[2]);
    ptrdiff_t KU = (ptrdiff_t)mxGetScalar(prhs[3]);
    const mxArray *RHS = prhs[4];

    if (!mxIsDouble(ABm))
        mexErrMsgIdAndTxt("dgbsv_solve_cached:ABf", "ABf must be double.");
    if (!mxIsInt64(IPm))
        mexErrMsgIdAndTxt("dgbsv_solve_cached:IPIV", "IPIV must be int64.");

    mwSize N = mxGetN(ABm);                 /* ABf is [LDAB x N] */
    ptrdiff_t LDAB = 2*KL + KU + 1;
    if ((ptrdiff_t)mxGetM(ABm) != LDAB)
        mexErrMsgIdAndTxt("dgbsv_solve_cached:dim", "size(ABf,1) must equal 2*KL+KU+1.");
    if (mxGetM(RHS)*mxGetN(RHS) != N)
        mexErrMsgIdAndTxt("dgbsv_solve_cached:rhs", "rhs length must equal size(ABf,2).");

    ptrdiff_t *IPIV = (ptrdiff_t*)mxGetInt64s(IPm);
    ptrdiff_t Nl = (ptrdiff_t)N, NRHS = 1, LDB = (ptrdiff_t)N, INFO = 0;
    char trans = 'N';

    int ab_cplx  = mxIsComplex(ABm);
    int rhs_cplx = mxIsComplex(RHS);

    if (ab_cplx) {
        /* complex factor -> zgbtrs; ensure rhs is complex */
        plhs[0] = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
        mxComplexDouble *out = mxGetComplexDoubles(plhs[0]);
        const mxArray *Rin = RHS; mxArray *Rc = NULL;
        if (!rhs_cplx) {
            mxArray *a[1] = { (mxArray*)RHS };
            mexCallMATLAB(1, &Rc, 1, a, "complex");
            Rin = Rc;
        }
        memcpy(out, mxGetComplexDoubles(Rin), (size_t)N*sizeof(mxComplexDouble));
        if (Rc) mxDestroyArray(Rc);
        zgbtrs(&trans, &Nl, &KL, &KU, &NRHS,
               (double*)mxGetComplexDoubles(ABm), &LDAB, IPIV,
               (double*)out, &LDB, &INFO);
    } else if (rhs_cplx) {
        /* real factor, complex rhs -> solve real & imaginary parts separately */
        plhs[0] = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
        mxComplexDouble *rz = mxGetComplexDoubles(RHS);
        mxComplexDouble *out = mxGetComplexDoubles(plhs[0]);
        double *AB = mxGetDoubles(ABm);
        double *br = (double*)mxMalloc((size_t)N*sizeof(double));
        double *bi = (double*)mxMalloc((size_t)N*sizeof(double));
        for (mwSize i = 0; i < N; i++) { br[i] = rz[i].real; bi[i] = rz[i].imag; }
        dgbtrs(&trans, &Nl, &KL, &KU, &NRHS, AB, &LDAB, IPIV, br, &LDB, &INFO);
        if (INFO == 0)
            dgbtrs(&trans, &Nl, &KL, &KU, &NRHS, AB, &LDAB, IPIV, bi, &LDB, &INFO);
        for (mwSize i = 0; i < N; i++) { out[i].real = br[i]; out[i].imag = bi[i]; }
        mxFree(br); mxFree(bi);
    } else {
        /* real factor, real rhs -> dgbtrs */
        plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
        double *out = mxGetDoubles(plhs[0]);
        memcpy(out, mxGetDoubles(RHS), (size_t)N*sizeof(double));
        dgbtrs(&trans, &Nl, &KL, &KU, &NRHS, mxGetDoubles(ABm), &LDAB, IPIV,
               out, &LDB, &INFO);
    }

    if (INFO != 0)
        mexErrMsgIdAndTxt("dgbsv_solve_cached:gbtrs",
            "*gbtrs returned INFO=%td.", INFO);
}
