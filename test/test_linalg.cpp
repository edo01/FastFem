#include "FastFem/linalg/Vector.hpp"
#include "FastFem/linalg/sparseMatrices/COOMatrix.hpp"

using namespace fastfem::linalg;

int main()
{

  Vector v1 = {1.0, 2.0, 3.0};
  Vector v2 = {4.0, 5.0, 6.0};

  Vector v3 = v1 + v2;
  Vector v4 = v1 - v2;
  Vector v5 = v1 * 2.0;

  v3.print();
  v4.print();
  v5.print();

  double dot = v1.dot(v2);
  double norm = v1.norm();

  std::cout << "v1 dot v2 = " << dot << std::endl;
  std::cout << "norm(v1) = " << norm << std::endl;

  COOMatrix M(3, 3, 7);

  std::cout << "nnz = " << M.nnz() << std::endl;

  M.add_entry(0, 0, 1.0);
  M.add_entry(2, 1, 2.0);
  M.add_entry(0, 2, 3.0);

  Vector Mv = M * v1;
  Mv.print();

  std::cout << "nnz = " << M.nnz() << std::endl;

  return 1;
}