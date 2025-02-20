#include "FastFem/linalg/sparseMatrices/SymCSRMatrix.hpp"
#include <algorithm>

namespace fastfem{
namespace linalg{

using types::ff_index;

SymCSRMatrix::SymCSRMatrix(ff_index n_cols, const CSRPattern& pattern) :
  CSRMatrix(n_cols, pattern){};

const double &SymCSRMatrix::get_entry(ff_index i, ff_index j) const
{
    if(i > j){
        std::swap(i, j);
    }

    ff_index row_start = base_pattern->row_ptr[i];
    ff_index row_end = base_pattern->row_ptr[i + 1];

    for(ff_index k = row_start; k < row_end; ++k){
        if(base_pattern->col_indices[k] == j){
            return values[k];
        }
    }

    static double dummy = 0.0;
    return dummy;
}

void SymCSRMatrix::accumulate_entry(ff_index i, ff_index j, double value)
{
    if(i > j){
        std::swap(i, j);
    }

    ff_index row_start = base_pattern->row_ptr[i];
    ff_index row_end = base_pattern->row_ptr[i + 1];

    for(ff_index k = row_start; k < row_end; ++k){
        if(base_pattern->col_indices[k] == j){
            values[k] += value;
            return;
        }
    }

    throw std::invalid_argument("SymCSRMatrix::accumulate_entry(): entry not found");
}

void SymCSRMatrix::operator=(const double& value){
    std::fill(values.begin(), values.end(), value);
}

Vector SymCSRMatrix::gemv(const Vector& x) const {

    if(n_cols != x.size()){
        throw std::invalid_argument("CSRMatrix::gemv(): incompatible dimensions");
    }

    Vector y(n_rows);
    const auto& row_ptr = base_pattern->row_ptr;
    const auto& col_indices = base_pattern->col_indices;

    for(ff_index i = 0; i < n_rows; ++i){
        ff_index row_start = row_ptr[i];
        ff_index row_end = row_ptr[i + 1];
        for(ff_index k = row_start; k < row_end; ++k){
            y[i] += values[k] * x[col_indices[k]];
        }
        for(ff_index k = row_start; k < row_end; ++k){
            ff_index j = col_indices[k];
            if(i != j){
                y[j] += values[k] * x[i];
            }
        }

    }

    return y;
}

} // namespace linalg
} // namespace FastFem
