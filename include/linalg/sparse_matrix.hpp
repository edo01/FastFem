#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

namespace linalg
{
/* A coefficient of a sparse matrix */
struct Coeff 
{
	int i;
	int j;
	double val;
};

/* A Sparse matrix is an (order independent) array of coeffs,
 * the non zero entries of the matrix.
 * That for is called a COO sparse matrixe (COO for coordinates).
 * It is the simplest form and not the most efficient.
 * In the main repo we use a (probably) better approach.
 */
struct SparseMatrix 
{
	int rows;
	int cols;
	int nnz;
	struct Coeff *coeffs;
};

} // linalg

#endif // SPARSE_MATRIX_HPP