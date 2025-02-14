#include "FastFem/linalg/Vector.hpp"
#include "FastFem/linalg/sparseMatrices/COOMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"

using namespace fastfem::linalg;

int main()
{

  Vector v1 = {1.0, 2.0, 3.0};
  Vector v2 = {4.0, 5.0, 6.0};

  CSRPattern pattern = {{0, 2, 3, 5}, {0, 1, 2, 0, 2}};
  CSRMatrix M(3, pattern);

  std::cout << "n_rows: " << M.get_n_rows() << std::endl;
  std::cout << "n_cols: " << M.get_n_cols() << std::endl;
  std::cout << "nnz: " << M.nnz() << std::endl;

  M.add_entry(0, 2);
  M.add_entry(1, 3);
  M.add_entry(2, -2);
  M.add_entry(3, 1);
  M.add_entry(4, 3);

  M.print();

  Vector Mv = M * v1;
  Mv.print();

  return 1;
}