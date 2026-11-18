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

/* dgbsv_factor.c
 *
 * Banded LU (LU factorization) FACTORIZATION step (LAPACK *gbtrf). Produces a
 * factor that is cached and reused by dgbsv_solve_cached. Companion to
 * dgbsv_solve.c, which factors and solves in one call. Splitting factor from
 * solve lets the nonlinear HNS iteration factor each per-mode operator once
 * and reuse it across Picard iterations, instead of re-factoring every time.
 *
 *   [ABf, IPIV, KL, KU] = dgbsv_factor(ML)            % auto-detect kl, ku
 *   [ABf, IPIV, KL, KU] = dgbsv_factor(ML, kl, ku)    % user-supplied bandwidth
 *
 * ML  : sparse N-by-N double matrix (real or complex)
 * ABf : [LDAB x N] double matrix (real if ML is real, complex if ML is
 *       complex), holding the LAPACK banded LU factors, ready for *gbtrs.
 *       KL and KU are the lower and upper bandwidth; LDAB, the leading
 *       dimension of the banded array, equals 2*KL+KU+1.
 * IPIV: [N x 1] int64 vector of pivot indices.
 *
 * Compile via src/linsolve/lapack_band/mexbuild_dgbsv.m (-lmwlapack).
 */

#include "mex.h"
#include "lapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1 && nrhs != 3)
        mexErrMsgIdAndTxt("dgbsv_factor:nrhs",
                          "Usage: [ABf,IPIV,KL,KU] = dgbsv_factor(ML[,kl,ku])");
    const mxArray *ML = prhs[0];
    if (!mxIsSparse(ML) || !mxIsDouble(ML))
        mexErrMsgIdAndTxt("dgbsv_factor:ML", "ML must be sparse double.");
    mwSize N = mxGetN(ML);
    if (mxGetM(ML) != N)
        mexErrMsgIdAndTxt("dgbsv_factor:square", "ML must be square.");

    int is_cplx = mxIsComplex(ML);
    const mwIndex *Jc = mxGetJc(ML);
    const mwIndex *Ir = mxGetIr(ML);

    ptrdiff_t KL, KU;
    if (nrhs == 3) {
        KL = (ptrdiff_t)mxGetScalar(prhs[1]);
        KU = (ptrdiff_t)mxGetScalar(prhs[2]);
    } else {
        ptrdiff_t klm = 0, kum = 0;
        for (mwSize j = 0; j < N; j++)
            for (mwIndex p = Jc[j]; p < Jc[j+1]; p++) {
                ptrdiff_t d = (ptrdiff_t)j - (ptrdiff_t)Ir[p];
                if (d > kum) kum = d;
                if (-d > klm) klm = -d;
            }
        KL = klm; KU = kum;
    }
    ptrdiff_t LDAB = 2*KL + KU + 1;

    /* Allocate the factored band as a MATLAB array (so it can be cached). */
    mxArray *ABm = mxCreateDoubleMatrix((mwSize)LDAB, N, is_cplx ? mxCOMPLEX : mxREAL);

    /* Fill the band: AB[KL+KU+i-j, j] = A[i,j] (column-major storage, 0-indexed). */
    if (is_cplx) {
        mxComplexDouble *AB = mxGetComplexDoubles(ABm);   /* zeroed */
        mxComplexDouble *Pr = mxGetComplexDoubles(ML);
        for (mwSize j = 0; j < N; j++)
            for (mwIndex p = Jc[j]; p < Jc[j+1]; p++) {
                ptrdiff_t row = KL + KU + (ptrdiff_t)Ir[p] - (ptrdiff_t)j;
                AB[(size_t)j*(size_t)LDAB + (size_t)row] = Pr[p];
            }
    } else {
        double *AB = mxGetDoubles(ABm);
        double *Pr = mxGetDoubles(ML);
        for (mwSize j = 0; j < N; j++)
            for (mwIndex p = Jc[j]; p < Jc[j+1]; p++) {
                ptrdiff_t row = KL + KU + (ptrdiff_t)Ir[p] - (ptrdiff_t)j;
                AB[(size_t)j*(size_t)LDAB + (size_t)row] = Pr[p];
            }
    }

    mxArray *IPm = mxCreateNumericMatrix(N, 1, mxINT64_CLASS, mxREAL);
    ptrdiff_t *IPIV = (ptrdiff_t*)mxGetInt64s(IPm);

    ptrdiff_t M = (ptrdiff_t)N, INFO = 0;
    if (is_cplx)
        zgbtrf(&M, &M, &KL, &KU, (double*)mxGetComplexDoubles(ABm), &LDAB, IPIV, &INFO);
    else
        dgbtrf(&M, &M, &KL, &KU, mxGetDoubles(ABm), &LDAB, IPIV, &INFO);

    if (INFO != 0) {
        mxDestroyArray(ABm); mxDestroyArray(IPm);
        mexErrMsgIdAndTxt("dgbsv_factor:gbtrf",
            "*gbtrf returned INFO=%td (>0 means singular U at that pivot).", INFO);
    }

    plhs[0] = ABm;
    plhs[1] = IPm;
    if (nlhs > 2) plhs[2] = mxCreateDoubleScalar((double)KL);
    if (nlhs > 3) plhs[3] = mxCreateDoubleScalar((double)KU);
}
