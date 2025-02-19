#include "FastFem/linalg/sparseMatrices/COOMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/SymCSRMatrix.hpp"

#include <memory>

using namespace fastfem::linalg;

int main()
{

    COOMatrix M = COOMatrix::random(10, 10, 0.7, -1.0, 1.0, true);
    M.print("COO Matrix");

    CSRMatrix A = M.to_CSR();
    // A.set_row_col_to_zero(1);
    A.print("CSR Matrix");
}