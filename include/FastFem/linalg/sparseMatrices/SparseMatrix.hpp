#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP
#include <FastFem/linalg/Vector.hpp>

namespace fastfem{
namespace linalg{

class SparseMatrix {
private:
    size_t n_rows;
    size_t n_cols;
    size_t nnz;

public:
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