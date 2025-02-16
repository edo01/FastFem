#include "FastFem/linalg/Vector.hpp"
#include "FastFem/linalg/sparseMatrices/COOMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/SymCSRMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/SkylineMatrix.hpp"

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

  std::cout << "---------- SkylineMatrix ----------" << std::endl;

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

  return 1;
}