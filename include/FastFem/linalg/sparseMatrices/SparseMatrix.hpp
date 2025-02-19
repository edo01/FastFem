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

    //read-only access with bounds checking
    const double &operator()(size_t i, size_t j) const;

    //no bound checking, used in assembly
    virtual void set_entry(size_t i, size_t j, double value) = 0;
    virtual void accumulate_entry(size_t i, size_t j, double value) = 0;

    virtual inline size_t nnz() const = 0;

    void print(const std::string& name = "") const;

    size_t get_n_rows() const { return n_rows; }
    size_t get_n_cols() const { return n_cols; }

    virtual inline bool is_symmetric() const { return false; } // override in symmetric matrices

    virtual void set_row_col_to_zero(size_t i);

private:
    virtual const double &get_entry(size_t i, size_t j) const = 0;

    bool check_bounds(size_t i, size_t j) const;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_SPARSEMATRIX_HPP