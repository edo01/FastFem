#include "FastFem/linalg/sparseMatrices/COOMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/SymCSRMatrix.hpp"
#include "FastFem/linalg/iterativeSolvers/CGSolver.hpp"

#include <memory>

using namespace fastfem::linalg;

/**
 * @brief This test verifies the correct functioning of the CG algorithm
 * 
 * @return int 
 */
int main()
{

  size_t n = 5;

  COOMatrix M = COOMatrix::random(n, n, 0.5, 1.0, 1.0, true);
  std::cout << "COO nnz: " << M.nnz() << std::endl;

  SymCSRMatrix A = M.to_CSR();
  std::cout << "CSR nnz: " << A.nnz() << std::endl;

  Vector v(n, 1.0);
  Vector Av = A * v;

  Av.print("Av: ");

  CGSolver solver;
  Vector x = solver.solve(A, Av);

  x.print("Solution");

  
  return 1;
}