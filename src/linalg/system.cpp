#include "linalg/system.hpp"
#include <cmath>
#include <iostream>

namespace linalg
{

/*     void cholesky_factorization(const SKRMatrix &M, SKRMatrix &L)
    {
        L = M;
        for(int i = 0; i < L.n_rows; i++)
        {
            std::cout << i << std::endl;
            double sum = 0;
            for(int j = 0 ; j<=i-L.offset[i];j++) // primo punto strano
            {
                sum = 0;
                for(int k = 0; k<j;k++)
                    sum+= L.data[L.offset[i]+k]*L.data[L.offset[i]+k];
                L.data[L.offset[i]+j] = (L.data[L.offset[i]+j]-sum)/L.data[L.offset[j]+j];
            }
            sum = 0;
            for(int k = 0; k<=i-L.offset[i];k++) // primo punto strano
                sum+=L.data[L.offset[i]+k]*L.data[L.offset[i]+k];
            L.data[L.offset[i]+i] = std::sqrt(L.data[L.offset[i]+i] -sum);
        }
    } */

    void cholesky_factorization(const SKRMatrix &M, SKRMatrix &L)
    {
        L = M;
        for(int i = 0; i < L.n_rows; i++)
        {
            std::cout << i << std::endl;

            int row_i_start  = L.offset[i];
            int nnz_row_i = i - L.col_index[i]; // + 1 
            
            double sum;
            for(int row_offset = 0; row_offset < nnz_row_i ; row_offset++) // primo punto strano
            {    
                int row_j_start  = L.offset[L.col_index[i] + row_offset]; // j 

                sum = 0;
                for(int k = 0; k < row_offset; k++)
                    // L_ik * L_jk
                    sum += L.data[row_i_start + k] * L.data[row_j_start + k];

                L[row_i_start + row_offset] = (L[row_i_start + row_offset] - sum) / L[row_j_start + row_offset];
            }
            
            sum = 0;
            // for(int k = 0; k <= i-L.offset[i]; k++) // primo punto strano
            //     sum+=L.data[L.offset[i]+k]*L.data[L.offset[i]+k];
            // L.data[L.offset[i]+i] = std::sqrt(L.data[L.offset[i]+i] -sum);
            for(int row_offset=0; row_offset < nnz_row_i; row_offset++)
                sum += L[row_i_start + row_offset] * L[row_i_start + row_offset];
            L[row_i_start + i] = std::sqrt(L[row_i_start + i] - sum);
        }
    }


    /******************************************************************************
     * Computes the product M * v when M is a SparseMatrix
     *****************************************************************************/
    void matrix_vector_product(const struct SparseMatrix *M, const double *v,
                               double *Mv)
    {

        memset(Mv, 0, M->rows * sizeof(double));

        for (size_t index = 0; index < M->nnz; index++)
        {
            // MV[i] += M_ij * v[j]
            Mv[M->coeffs[index].i] += M->coeffs[index].val * v[M->coeffs[index].j];
        }
    }

    /* Create an (unitialized) array of N double precision floating point values */
    double *array(int N) { return (double *)malloc(N * sizeof(double)); }

    double two_product(double x, double y, double &err)
    {
        double p = x * y;
        err = x * y - p;
        return p; // product
    }

    double two_sum(double x, double y, double &err)
    {
        double s = x + y;
        double z = s - x;
        err = (x - (s - z)) + (y - z);
        return s;
    }

    /* Vector product between two vectors in dim N */
    double compensated_dot(const double *A, const double *B, int N)
    {

        double res = 0;
        double err_tot = 0;
        double err_p, err_s;
        double p_temp;

        double product = two_product(A[0], B[0], err_tot); // initialize the product and the error

        for (int i = 1; i < N; i++)
        {
            p_temp = two_product(A[i], B[i], err_p);

            product = two_sum(p_temp, product, err_s);

            err_tot += err_p + err_s;
        }

        return res + err_tot;
    }

    double blas_dot(const double *A, const double *B, int N)
    {
        double res = 0;
#pragma omp parallel for reduction(+ : res)
        for (int i = 0; i < N; i++)
        {
            res += A[i] * B[i];
        }
        return res;
    }

    /* aX + bY -> Y  (axpby reads as aX plus bY)
     * a and b are scalar, X and Y are vectors in dim N
     */
    void blas_axpby(double a, const double *X, double b, double *Y, int N)
    {
#pragma omp parallel for
        for (int i = 0; i < N; i++)
        {
            Y[i] = a * X[i] + b * Y[i];
        }
    }

    /******************************************************************************
     * Solving AU=B where A is SPD of size NxN using steepest descent method.
     * Find it back on a sheet of paper, not on Google !
     * One minimizes the functional 1/2 <AU,U> - <B,U>.
     * The minor peculiarity here is that A = S + M and we do not wish to
     * add these two sparse matrices up-front but simply compute AU as SU + MU
     * wherever needed.
     *****************************************************************************/
    int gradient_system_solve(const struct SparseMatrix *S,
                              const struct SparseMatrix *M, const double *B,
                              double *U, int N)
    {
        /* Set U_0 to the initial guess */
        memset(U, 0, N * sizeof(double));

        /* computing the first residual r_0 = B - (M+S)U_0 = B */
        double *r = array(N);
        memcpy(r, B, N * sizeof(double));

        /* l^2 (squared) norm of the residue */
        double error2 = blas_dot(r, r, N);

        /*
         * Temporaries to store matrix-vector mul and avoid
         * multiple computation of the same product
         */
        double *Mr = array(N);
        double *Ar = array(N);

        double tol2 = 1e-6;
        int max_iter = 1000;
        int iterate = 0;
        while (error2 > tol2 && iterate < max_iter)
        {

            iterate++;

            /* Compute Ar_k=Sr_k+Mr_k */
            matrix_vector_product(S, r, Ar);
            matrix_vector_product(M, r, Mr);
            blas_axpby(1, Mr, 1, Ar, N);

            /* Compute alpha_k */
            double alpha = blas_dot(r, r, N) / blas_dot(Ar, r, N);

            /* Update the new solution */
            blas_axpby(alpha, r, 1, U, N);

            /* Update the new residual using r= r-a_k*Ar which allows to
             * avoid one matrix-vector product */
            blas_axpby(-alpha, Ar, 1, r, N);

            /* Update error2 */
            error2 = blas_dot(r, r, N);
        }
        /* Release system memory */
        free(Ar);
        free(Mr);
        free(r);
        return iterate;
    }

    int CG_system_solve(const struct SparseMatrix *S,
                        const struct SparseMatrix *M, const double *B,
                        double *U, int N)
    {
        /* Set U_0 to the initial guess */
        memset(U, 0, N * sizeof(double));

        double *r = array(N);
        memcpy(r, B, N * sizeof(double));

        double *p = array(N);
        memcpy(p, r, N * sizeof(double));

        double *Ap = array(N);
        double *Mp = array(N);

        double error2 = blas_dot(r, r, N);

        double tol2 = 1e-6;
        int max_iter = 1000;
        int iterate = 0;

        while (error2 > tol2 && iterate < max_iter)
        {
            matrix_vector_product(S, p, Ap);
            matrix_vector_product(M, p, Mp);
            blas_axpby(1, Mp, 1, Ap, N);

            double alpha = blas_dot(r, r, N) / blas_dot(p, Ap, N); // alpha_k = (r_k, r_k) / (p_k, A * p_k)

            blas_axpby(alpha, p, 1, U, N); // u_k+1 = u_k + alpha_k * p_k

            blas_axpby(-alpha, Ap, 1, r, N); // r_k+1 = r_k - alpha_k * A * p_k

            // if(blas_dot(r, r, N) < tol2) break;

            double beta = blas_dot(r, r, N) / error2; // beta_k = (r_k+1, r_k+1) / (r_k, r_k)

            blas_axpby(1, r, beta, p, N); // p_k+1 = r_k+1 + beta_k * p_k

            error2 = blas_dot(r, r, N);

            iterate++;
        }

        free(r);
        free(p);
        free(Ap);
        free(Mp);

        return iterate;
    }

} // namespace linalg
