#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <stdlib.h>
#include <string.h>
#include "linalg/sparse_matrix.hpp"
#include "linalg/vector.hpp"

namespace linalg
{

/******************************************************************************
 * Computes the product M * v when M is a SparseMatrix
 *****************************************************************************/
void matrix_vector_product(const struct SparseMatrix *M, const double *v,
			   double *Mv);

/* Create an (unitialized) array of N double precision floating point values */
double *array(int N);

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
			  double *U, int N);

int CG_system_solve(const struct SparseMatrix *S,
        const struct SparseMatrix *M, const double *B,
        double *U, int N);

void blas_axpby(double a, const double *X, double b, double *Y, int N);

double blas_dot(const double *A, const double *B, int N);


}

#endif // SYSTEM_HPP