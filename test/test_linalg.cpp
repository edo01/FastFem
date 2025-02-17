#include "FastFem/linalg/Vector.hpp"
#include "FastFem/linalg/sparseMatrices/COOMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/SymCSRMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/SkylineMatrix.hpp"
#include "FastFem/linalg/iterativeSolvers/CGSolver.hpp"

using namespace fastfem::linalg;

int main()
{

  std::cout << "---------- CSRMatrix ----------" << std::endl;

  Vector v1 = {1.0, 2.0, 3.0};

  CSRPattern pattern = {{0, 2, 3, 5}, {0, 1, 2, 0, 2}};
  CSRMatrix M(3, pattern);

  std::cout << "n_rows: " << M.get_n_rows() << std::endl;
  std::cout << "n_cols: " << M.get_n_cols() << std::endl;
  std::cout << "nnz: " << M.nnz() << std::endl;

  M.print_pattern();

  M.add_entry(0, 2);
  M.add_entry(1, 3);
  M.add_entry(2, -2);
  M.add_entry(3, 1);
  M.add_entry(4, 3);

  M.print();

  Vector Mv = M * v1;
  Mv.print();

  std::cout << "---------- SymCSRMatrix ----------" << std::endl;

  Vector v2 = {1.0, 2.0, 3.0, 4.0};

  CSRPattern pattern2 = {{0, 2, 3, 4, 5}, {1, 3, 1, 2, 3}};
  SymCSRMatrix M2(4, pattern2);

  std::cout << "n_rows: " << M2.get_n_rows() << std::endl;
  std::cout << "n_cols: " << M2.get_n_cols() << std::endl;
  std::cout << "nnz: " << M2.nnz() << std::endl;

  M2.print_pattern();

  M2.add_entry(0, 2);
  M2.add_entry(1, 3);
  M2.add_entry(2, -2);
  M2.add_entry(3, 1);
  M2.add_entry(4, 3);

  M2.print();

  Vector Mv2 = M2 * v2;
  Mv2.print();

  /*std::cout << "---------- SkylineMatrix ----------" << std::endl;

  std::vector<size_t> skyline = {0, 1, 3, 6, 9, 12};
  SkylineMatrix M3(5, skyline);
  M3.add_entry(0, 4);
  M3.add_entry(1, 2);
  M3.add_entry(2, 5);
  M3.add_entry(3, 1);
  M3.add_entry(4, 3);
  M3.add_entry(5, 6);
  M3.add_entry(6, 1);
  M3.add_entry(7, 4);
  M3.add_entry(8, 7);
  M3.add_entry(9, 2);
  M3.add_entry(10, 5);
  M3.add_entry(11, 8);

  std::cout << "n_rows: " << M3.get_n_rows() << std::endl;
  std::cout << "n_cols: " << M3.get_n_cols() << std::endl;
  std::cout << "nnz: " << M3.nnz() << std::endl;
  std::cout << "Printing pattern:" << std::endl << std::endl;
  M3.print_pattern();
  std::cout << std::endl;
  std::cout << "Printing matrix:" << std::endl << std::endl;
  M3.print();
  std::cout << std::endl;

  // Define Right-Hand Side Vector b
  Vector b = {15.0, 10.0, 25.0, 30.0, 40.0};
  std::cout << "Vector b: " << std::endl;
  b.print();
  std::cout << std::endl;

  // Step 1: Factorize A (Cholesky decomposition)
  std::cout << "Performing Cholesky factorization..." << std::endl;
  M3.cholesky_factorize();

  // Step 2: Solve Ax = b
  std::cout << "Solving Ax = b using Cholesky factorization..." << std::endl;
  Vector x = M3.cholesky_solve(b);

  // Print the solution
  std::cout << "Solution vector x:" << std::endl;
  x.print();
  std::cout << std::endl;

  // Step 3: Verify Ax ≈ b
  std::cout << "Verifying solution: Computing A * x..." << std::endl;
  Vector Ax = M3.gemv(x);
  Ax.print();
*/

  std::cout << "---------- Testing Skyline Cholesky ----------" << std::endl;

  // Define the skyline structure for A
  std::vector<size_t> skyline = {0, 1, 3, 5, 7, 9};
  SkylineMatrix A(5, skyline);

  // Fill values (lower triangular part)
  /*std::vector<double> values = {4, 1, 4, 1, 4, 1, 4, 1, 4};
  for (size_t i = 0; i < values.size(); ++i) {
      A.add_entry(i, values[i]);
  }*/
  A.insert_entry(0, 0, 4);
  A.insert_entry(0, 1, 1);
  A.insert_entry(1, 1, 4);
  A.insert_entry(1, 2, 1);
  A.insert_entry(2, 2, 4);
  A.insert_entry(2, 3, 1);
  A.insert_entry(3, 3, 4);
  A.insert_entry(3, 4, 1);
  A.insert_entry(4, 4, 4);

  // Copy a for testing reasona, considring A will be modified
  SkylineMatrix A2 = A;

  std::cout << "Matrix pattern:" << std::endl;
  A.print_pattern();
  std::cout << "Original matrix:" << std::endl;
  A.print();

  // Define right-hand side vector b
  Vector b = {5.0, 6.0, 7.0, 8.0, 9.0};

  // Solve Ax = b using Cholesky
  A.cholesky_factorize();

  // Print the computed lower triangular matrix L
  std::cout << "Cholesky factor L (lower triangular part):" << std::endl;
  A.print();

  Vector x = A.cholesky_solve(b);

  // Print solution
  std::cout << "Solution vector x:" << std::endl;
  x.print();

  // Verify Ax ≈ b
  std::cout << "Computing A * x..." << std::endl;
  Vector Ax = A2.gemv(x);
  Ax.print();

  return 1;
}