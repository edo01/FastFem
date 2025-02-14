#include <FastFem/linalg/sparseMatrices/COOMatrix.hpp>

namespace fastfem{
namespace linalg{

COOMatrix::COOMatrix(size_t n_rows, size_t n_cols, size_t nnz_hint) : SparseMatrix(n_rows, n_cols) {
    row_indices.reserve(nnz_hint);
    col_indices.reserve(nnz_hint);
    values.reserve(nnz_hint);
}

const double &COOMatrix::get_entry(size_t i, size_t j) const {

    static double dummy = 0.0;

    for(size_t k = 0; k < nnz(); ++k){
        if(row_indices[k] == i && col_indices[k] == j){
            return values[k];
        }
    }   
    
    return dummy;
}

Vector COOMatrix::gemv(const Vector& x) const{
    if(x.size() != n_cols){
        throw std::invalid_argument("COOMatrix::gemv(): incompatible dimensions");
    }

    Vector y(n_rows);
    for(size_t i = 0; i < n_rows; ++i){
        y[i] = 0.0;
    }
    for(size_t k = 0; k < nnz(); ++k){
        y[row_indices[k]] += values[k] * x[col_indices[k]];
    }
    return y;
}

void COOMatrix::add_entry(size_t i, size_t j, double value){
    if(!check_bounds(i, j)){
        throw std::out_of_range("COOMatrix::add_entry(): index out of range");
    }
    row_indices.push_back(i);
    col_indices.push_back(j);
    values.push_back(value);
}

} // namespace linalg
} // namespace FastFem