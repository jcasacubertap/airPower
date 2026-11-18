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

/* dgbsv_solve.c
 *
 * MEX (MATLAB executable, compiled C) wrapper for LAPACK *gbsv: a banded LU
 * (LU factorization) solve with partial pivoting.
 *   Real  ML -> dgbsv
 *   Complex ML -> zgbsv
 *
 *   x = dgbsv_solve(ML, rhs)         % auto-detect kl, ku
 *   x = dgbsv_solve(ML, rhs, kl, ku) % user-supplied bandwidth (kl = lower, ku = upper)
 *
 * ML must be a sparse N-by-N double matrix (real or complex).
 * rhs (right-hand side) must be a dense N-by-1 vector (real or complex;
 * complex if ML is).
 *
 * Compile via src/linsolve/lapack_band/mexbuild_dgbsv.m (-R2018a -lmwlapack).
 */

#include "mex.h"
#include "lapack.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2 && nrhs != 4) {
        mexErrMsgIdAndTxt("dgbsv_solve:nrhs",
                          "Usage: x = dgbsv_solve(ML, rhs[, kl, ku])");
    }
    const mxArray *ML  = prhs[0];
    const mxArray *RHS = prhs[1];

    if (!mxIsSparse(ML) || !mxIsDouble(ML)) {
        mexErrMsgIdAndTxt("dgbsv_solve:ML", "ML must be sparse double.");
    }
    if (mxIsSparse(RHS) || !mxIsDouble(RHS)) {
        mexErrMsgIdAndTxt("dgbsv_solve:rhs", "rhs must be dense double.");
    }
    mwSize N = mxGetN(ML);
    if (mxGetM(ML) != N) {
        mexErrMsgIdAndTxt("dgbsv_solve:square", "ML must be square.");
    }
    if (mxGetM(RHS) * mxGetN(RHS) != N) {
        mexErrMsgIdAndTxt("dgbsv_solve:dim", "rhs length must equal size(ML,1).");
    }

    /* Use complex solve if EITHER ML or RHS is complex. This handles the
     * common case in nonlinear HNS where the MFD (0,0) ML is real but the
     * NLT-projected RHS is stored as complex. */
    int is_cplx = mxIsComplex(ML) || mxIsComplex(RHS);
    int ML_cplx = mxIsComplex(ML);   /* needed when filling AB */

    const mwIndex *Jc = mxGetJc(ML);
    const mwIndex *Ir = mxGetIr(ML);

    /* --- Detect bandwidth --- */
    ptrdiff_t KL, KU;
    if (nrhs == 4) {
        KL = (ptrdiff_t)mxGetScalar(prhs[2]);
        KU = (ptrdiff_t)mxGetScalar(prhs[3]);
    } else {
        ptrdiff_t kl_max = 0, ku_max = 0;
        for (mwSize j = 0; j < N; j++) {
            for (mwIndex p = Jc[j]; p < Jc[j+1]; p++) {
                ptrdiff_t i = (ptrdiff_t)Ir[p];
                ptrdiff_t d = (ptrdiff_t)j - i;
                if (d >  ku_max) ku_max = d;
                if (-d > kl_max) kl_max = -d;
            }
        }
        KL = kl_max;
        KU = ku_max;
    }

    ptrdiff_t LDAB = 2*KL + KU + 1;
    size_t AB_elem = (size_t)LDAB * (size_t)N;
    size_t elem_bytes = is_cplx ? 16 : 8;
    double GB = (double)AB_elem * (double)elem_bytes / 1e9;

    mexPrintf("[dgbsv] %s N=%td  KL=%td  KU=%td  LDAB=%td  AB=%.2f GB\n",
              is_cplx ? "ZGBSV(complex)" : "DGBSV(real)",
              (ptrdiff_t)N, KL, KU, LDAB, GB);
    mexEvalString("drawnow;");  /* flush stdout before the big allocation below */

    /* --- Allocate banded matrix (zeroed) --- */
    double *AB = (double*)mxCalloc(AB_elem, elem_bytes);
    if (!AB) {
        mexErrMsgIdAndTxt("dgbsv_solve:oom",
            "Failed to allocate AB (%.2f GB).", GB);
    }

    /* --- Fill AB ---
     * Banded format (column-major): AB[KL+KU+i-j, j] = A[i,j], 0-indexed.
     * For complex: AB is stored interleaved (real, imaginary, real, imaginary, ...),
     * as 2*LDAB*N doubles, where LDAB (leading dimension of the banded array) = 2*KL+KU+1.
     */
    if (is_cplx) {
        /* AB is complex (interleaved). If ML is complex, copy re/im;
         * if ML is real (because RHS forced is_cplx), copy real, imag=0
         * (AB is already zeroed by mxCalloc). */
        if (ML_cplx) {
            mxComplexDouble *Pr_z = mxGetComplexDoubles(ML);
            if (!Pr_z) {
                mxFree(AB);
                mexErrMsgIdAndTxt("dgbsv_solve:cplx", "mxGetComplexDoubles failed.");
            }
            for (mwSize j = 0; j < N; j++) {
                for (mwIndex p = Jc[j]; p < Jc[j+1]; p++) {
                    ptrdiff_t i = (ptrdiff_t)Ir[p];
                    ptrdiff_t row_AB = KL + KU + i - (ptrdiff_t)j;
                    size_t off = ((size_t)j * (size_t)LDAB + (size_t)row_AB) * 2;
                    AB[off]   = Pr_z[p].real;
                    AB[off+1] = Pr_z[p].imag;
                }
            }
        } else {
            const double *Pr = mxGetDoubles(ML);
            for (mwSize j = 0; j < N; j++) {
                for (mwIndex p = Jc[j]; p < Jc[j+1]; p++) {
                    ptrdiff_t i = (ptrdiff_t)Ir[p];
                    ptrdiff_t row_AB = KL + KU + i - (ptrdiff_t)j;
                    size_t off = ((size_t)j * (size_t)LDAB + (size_t)row_AB) * 2;
                    AB[off]   = Pr[p];
                    /* AB[off+1] = 0 already */
                }
            }
        }
    } else {
        const double *Pr = mxGetDoubles(ML);
        for (mwSize j = 0; j < N; j++) {
            for (mwIndex p = Jc[j]; p < Jc[j+1]; p++) {
                ptrdiff_t i = (ptrdiff_t)Ir[p];
                ptrdiff_t row_AB = KL + KU + i - (ptrdiff_t)j;
                AB[(size_t)j * (size_t)LDAB + (size_t)row_AB] = Pr[p];
            }
        }
    }

    ptrdiff_t *IPIV = (ptrdiff_t*)mxMalloc((size_t)N * sizeof(ptrdiff_t));
    if (!IPIV) {
        mxFree(AB);
        mexErrMsgIdAndTxt("dgbsv_solve:oom", "Failed to allocate IPIV.");
    }

    /* --- Output / RHS buffer --- */
    if (is_cplx) {
        plhs[0] = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
    } else {
        plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    }

    ptrdiff_t NRHS = 1, INFO = 0, N_lp = (ptrdiff_t)N, LDB = (ptrdiff_t)N;

    if (is_cplx) {
        /* RHS: ensure complex */
        const mxArray *RHS_in = RHS;
        mxArray *RHS_cplx = NULL;
        if (!mxIsComplex(RHS)) {
            mxArray *args[1] = { (mxArray*)RHS };
            mexCallMATLAB(1, &RHS_cplx, 1, args, "complex");
            RHS_in = RHS_cplx;
        }
        mxComplexDouble *rhs_z = mxGetComplexDoubles(RHS_in);
        mxComplexDouble *out_z = mxGetComplexDoubles(plhs[0]);
        /* Copy rhs into output buffer (zgbsv operates in-place on B) */
        memcpy(out_z, rhs_z, N * sizeof(mxComplexDouble));
        if (RHS_cplx) mxDestroyArray(RHS_cplx);
        zgbsv(&N_lp, &KL, &KU, &NRHS, AB, &LDAB, IPIV,
              (double*)out_z, &LDB, &INFO);
    } else {
        const double *rhs_data = mxGetDoubles(RHS);
        double *B = mxGetDoubles(plhs[0]);
        memcpy(B, rhs_data, N * sizeof(double));
        dgbsv(&N_lp, &KL, &KU, &NRHS, AB, &LDAB, IPIV, B, &LDB, &INFO);
    }

    mxFree(AB);
    mxFree(IPIV);

    if (INFO != 0) {
        mexErrMsgIdAndTxt("dgbsv_solve:gbsv",
            "*gbsv returned INFO=%td (>0 means singular at that row).", INFO);
    }
}
