#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP
#include <FastFem/linalg/Vector.hpp>

namespace fastfem{
namespace linalg{

class SparseMatrix {
protected:
    size_t n_rows;
    size_t n_cols;
    size_t nnz;

public:
    SparseMatrix(size_t n_rows, size_t n_cols, size_t nnz) : n_rows(n_rows), n_cols(n_cols), nnz(nnz) {}
    virtual ~SparseMatrix() = default;
    virtual Vector gemv(const Vector& x) const = 0;

    Vector operator*(const Vector& x) const {
        return gemv(x);
    }
};

class SymmetricMatrix : virtual public SparseMatrix {
    // Marker class to enforce symmetry
};

} // namespace linalg
} // namespace FastFem

#endif // SPARSEMATRIX_HPP