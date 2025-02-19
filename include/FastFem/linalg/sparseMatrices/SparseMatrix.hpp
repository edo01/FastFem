#ifndef FASTFEM_SPARSEMATRIX_HPP
#define FASTFEM_SPARSEMATRIX_HPP
#include <FastFem/linalg/Vector.hpp>

namespace fastfem{
namespace linalg{

class SparseMatrix {
protected:
    size_t n_rows;
    size_t n_cols;

public:
    SparseMatrix(size_t n_rows, size_t n_cols);
    virtual ~SparseMatrix() = default;

    virtual Vector gemv(const Vector& x) const = 0;
    Vector operator*(const Vector& x) const;

    const double &operator()(size_t i, size_t j) const;
    virtual inline size_t nnz() const = 0;

    void print(const std::string& name = "") const;

    size_t get_n_rows() const { return n_rows; }
    size_t get_n_cols() const { return n_cols; }

    virtual inline bool is_symmetric() const { return false; } // override in symmetric matrices

private:
    virtual const double &get_entry(size_t i, size_t j) const = 0;

protected:
    bool check_bounds(size_t i, size_t j) const;

};

// Marker class to enforce symmetry
class SymmetricMatrix : public SparseMatrix {
    // inherit constructors
    using SparseMatrix::SparseMatrix;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_SPARSEMATRIX_HPP