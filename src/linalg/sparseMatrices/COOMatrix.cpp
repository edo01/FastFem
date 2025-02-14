#include <FastFem/linalg/sparseMatrices/COOMatrix.hpp>

namespace fastfem{
namespace linalg{

COOMatrix::COOMatrix(size_t n_rows, size_t n_cols, size_t nnz) : SparseMatrix(n_rows, n_cols, nnz) {
    row_indices.reserve(nnz);
    col_indices.reserve(nnz);
    values.reserve(nnz);
}

double &COOMatrix::operator()(size_t i, size_t j) {
    static double dummy = 0.0;
    if(i >= n_rows || j >= n_cols){
        throw std::out_of_range("COOMatrix::operator(): index out of range");
    }

    for(size_t k = 0; k < nnz; ++k){
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
    for(size_t k = 0; k < nnz; ++k){
        y[row_indices[k]] += values[k] * x[col_indices[k]];
    }
    return y;
}

} // namespace linalg
} // namespace FastFem