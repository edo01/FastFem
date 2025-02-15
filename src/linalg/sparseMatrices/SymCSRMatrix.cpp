#include "FastFem/linalg/sparseMatrices/SymCSRMatrix.hpp"
#include <algorithm>

namespace fastfem{
namespace linalg{

SymCSRMatrix::SymCSRMatrix(size_t n_cols, const CSRPattern& pattern) :
  CSRMatrix(n_cols, pattern){};

// const double &CSRMatrix::get_entry(size_t i, size_t j) const
// {
//     size_t row_start = base_pattern->row_ptr[i];
//     size_t row_end = base_pattern->row_ptr[i + 1];

//     for(size_t k = row_start; k < row_end; ++k){
//         if(base_pattern->col_indices[k] == j){
//             return values[k];
//         }
//     }

//     static double dummy = 0.0;
//     return dummy;
// }

// void CSRMatrix::add_entry(size_t index, double value)
// {
//     if(index >= base_pattern->col_indices.size()){
//         throw std::out_of_range("CSRMatrix::add_entry(): index out of range");
//     }
//     values[index] = value;
// }


// Vector CSRMatrix::gemv(const Vector& x) const {

//     if(n_cols != x.size()){
//         throw std::invalid_argument("CSRMatrix::gemv(): incompatible dimensions");
//     }

//     Vector y(n_rows);
//     const auto& row_ptr = base_pattern->row_ptr;
//     const auto& col_indices = base_pattern->col_indices;

//     for(size_t i = 0; i < n_rows; ++i){
//         size_t row_start = row_ptr[i];
//         size_t row_end = row_ptr[i + 1];
//         for(size_t k = row_start; k < row_end; ++k){
//             y[i] += values[k] * x[col_indices[k]];
//         }
//     }

//     return y;
// }

} // namespace linalg
} // namespace FastFem
