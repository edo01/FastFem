#include "FastFem/linalg/sparseMatrices/SparseMatrix.hpp"

#include <iostream>
#include <iomanip>

namespace fastfem{
namespace linalg{

using types::ff_index;

SparseMatrix::SparseMatrix(ff_index n_rows, ff_index n_cols) : n_rows(n_rows), n_cols(n_cols)
{
    if(n_rows == 0 || n_cols == 0){
        throw std::invalid_argument("SparseMatrix::SparseMatrix(): invalid dimensions");
    }
}

Vector SparseMatrix::operator*(const Vector& x) const
{
    return gemv(x);
}

bool SparseMatrix::check_bounds(ff_index i, ff_index j) const
{
    return i < n_rows && j < n_cols;
}

const double &SparseMatrix::operator()(ff_index i, ff_index j) const
{
    if(!check_bounds(i, j))
    {
        throw std::out_of_range("SparseMatrix::operator(): index out of range");
    }
    return get_entry(i, j);
}

void SparseMatrix::print(const std::string& name) const
{
    if (!name.empty()) {
        std::cout << name << std::endl;
    }

    for (ff_index i = 0; i < n_rows; ++i) {
        for (ff_index j = 0; j < n_cols; ++j) {
            std::cout << std::fixed << std::setprecision(2) << std::setw(5) << get_entry(i, j) << " ";
        }
        std::cout << '\n';  // Use '\n' for better performance over std::endl
    }
}

void SparseMatrix::set_row_col_to_zero(ff_index /*i*/)
{
    throw std::runtime_error("SparseMatrix::set_row_col_to_zero(): not implemented");
}

} // namespace linalg
} // namespace FastFem