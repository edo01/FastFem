#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include <algorithm>

namespace fastfem{
namespace linalg{

CSRPattern::CSRPattern(std::initializer_list<size_t> row_ptr, std::initializer_list<size_t> col_indices) : row_ptr(row_ptr), col_indices(col_indices)
{
    if(this->row_ptr.size() < 2){
        throw std::invalid_argument("CSRPattern::CSRPattern(): invalid row_ptr size");
    }
    if(this->col_indices.size() != this->row_ptr.back()){
        throw std::invalid_argument("CSRPattern::CSRPattern(): invalid col_indices size");
    }
}

CSRMatrix::CSRMatrix(size_t n_cols, const CSRPattern& pattern) :
  SparseMatrix(pattern.row_ptr.size() - 1, n_cols),
  base_pattern(std::make_shared<CSRPattern>(pattern)),
  values(pattern.col_indices.size())
{
    if(n_cols == 0 || *std::max_element(pattern.col_indices.begin(), pattern.col_indices.end()) >= n_cols){
        throw std::invalid_argument("CSRMatrix::CSRMatrix(): invalid n_cols");
    }
}

const double &CSRMatrix::get_entry(size_t i, size_t j) const
{
    size_t row_start = base_pattern->row_ptr[i];
    size_t row_end = base_pattern->row_ptr[i + 1];

    for(size_t k = row_start; k < row_end; ++k){
        if(base_pattern->col_indices[k] == j){
            return values[k];
        }
    }

    static double dummy = 0.0;
    return dummy;
}

void CSRMatrix::add_entry(size_t index, double value)
{
    if(index >= base_pattern->col_indices.size()){
        throw std::out_of_range("CSRMatrix::add_entry(): index out of range");
    }
    values[index] = value;
}


Vector CSRMatrix::gemv(const Vector& x) const {

    if(n_cols != x.size()){
        throw std::invalid_argument("CSRMatrix::gemv(): incompatible dimensions");
    }

    Vector y(n_rows);
    const auto& row_ptr = base_pattern->row_ptr;
    const auto& col_indices = base_pattern->col_indices;

    for(size_t i = 0; i < n_rows; ++i){
        size_t row_start = row_ptr[i];
        size_t row_end = row_ptr[i + 1];
        for(size_t k = row_start; k < row_end; ++k){
            y[i] += values[k] * x[col_indices[k]];
        }
    }

    return y;
}

} // namespace linalg
} // namespace FastFem
