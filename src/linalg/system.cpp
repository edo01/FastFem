#include "linalg/system.hpp"

namespace linalg
{
/******************************************************************************
 * Computes the product M * v when M is a SparseMatrix
 *****************************************************************************/
void matrix_vector_product(const struct SparseMatrix *M, const double *v,
			   double *Mv)
{

	// initialize Mv to 0
	for(size_t i=0; i < M->rows; i++){
		Mv[i] = 0;
	}

	for(size_t index=0; index < M->nnz; index++){
		// MV[i] += M_ij * v[j]
		Mv[M->coeffs[index].i] += M->coeffs[index].val * v[M->coeffs[index].j]; 
	}
}

/* Create an (unitialized) array of N double precision floating point values */
double *array(int N) { return (double *)malloc(N * sizeof(double)); }


double two_product(double x, double y, double &err)
{
	double p = x * y;
	err = x*y - p;
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
double compensated_dot(const double *A, const double *B, int N){

	double res = 0;
	double err_tot = 0;
	double err_p, err_s;
	double p_temp;

	double product = two_product(A[0], B[0], err_tot); // initialize the product and the error

	for(int i=1; i<N; i++){
		p_temp = two_product(A[i], B[i], err_p);
		
		product = two_sum(p_temp, product, err_s);

		err_tot += err_p + err_s;
	}

	return res + err_tot;
}

double blas_dot(const double *A, const double *B, int N)
{
	double res=0;
	for(int i=0; i<N; i++){
		res += A[i]*B[i];
	}
	return res;
}

/* aX + bY -> Y  (axpby reads as aX plus bY)
 * a and b are scalar, X and Y are vectors in dim N
 */
void blas_axpby(double a, const double *X, double b, double *Y, int N)
{
	for(int i=0; i<N; i++){
		Y[i] = a*X[i]+b*Y[i];
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
	while (error2 > tol2 && iterate < max_iter) {

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


}