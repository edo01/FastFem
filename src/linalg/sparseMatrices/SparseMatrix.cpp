#include "FastFem/linalg/sparseMatrices/SparseMatrix.hpp"

namespace fastfem{
namespace linalg{


SparseMatrix::SparseMatrix(size_t n_rows, size_t n_cols) : n_rows(n_rows), n_cols(n_cols)
{
    if(n_rows == 0 || n_cols == 0){
        throw std::invalid_argument("SparseMatrix::SparseMatrix(): invalid dimensions");
    }
}

Vector SparseMatrix::operator*(const Vector& x) const
{
    return gemv(x);
}

bool SparseMatrix::check_bounds(size_t i, size_t j) const
{
    return i < n_rows && j < n_cols;
}

const double &SparseMatrix::operator()(size_t i, size_t j) const
{
    if(!check_bounds(i, j))
    {
        throw std::out_of_range("SparseMatrix::operator(): index out of range");
    }
    return get_entry(i, j);
}

void SparseMatrix::print() const
{
    for(size_t i = 0; i < n_rows; ++i){
        for(size_t j = 0; j < n_cols; ++j){
            std::cout << get_entry(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

} // namespace linalg
} // namespace FastFem