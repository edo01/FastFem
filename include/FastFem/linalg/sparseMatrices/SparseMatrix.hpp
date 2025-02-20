#ifndef FASTFEM_SPARSEMATRIX_HPP
#define FASTFEM_SPARSEMATRIX_HPP
#include <FastFem/linalg/Vector.hpp>

namespace fastfem{
namespace linalg{

/**
 * @brief Abstract class for sparse matrices
 * 
 */
class SparseMatrix {
protected:
    size_t n_rows;
    size_t n_cols;

public:
    SparseMatrix(size_t n_rows, size_t n_cols);
    virtual ~SparseMatrix() = default;

    /**
     * @brief Matrix-vector product 
     */
    virtual Vector gemv(const Vector& x) const = 0;
    Vector operator*(const Vector& x) const;

    /**
     * @brief read-only access to the matrix entries with bounds checking
     */
    const double &operator()(size_t i, size_t j) const;

    /**
     * @brief set an entry with no bound checking, used in assembly
     */
    virtual void set_entry(size_t i, size_t j, double value) = 0;

    /**
     * @brief accumulate an entry with no bound checking, used in assembly
     */
    virtual void accumulate_entry(size_t i, size_t j, double value) = 0;

    virtual inline size_t nnz() const = 0;
    virtual inline bool is_symmetric() const { return false; } // override in symmetric matrices
    void print(const std::string& name = "") const;

    /**
     * @brief Set a row and a column to zero, used for the application of Dirichlet boundary conditions
     */
    virtual void set_row_col_to_zero(size_t i);

    size_t get_n_rows() const { return n_rows; }
    size_t get_n_cols() const { return n_cols; }

private:
    /**
     * @brief Read-only access to the matrix entries with no bound checking
     */
    virtual const double &get_entry(size_t i, size_t j) const = 0;

    bool check_bounds(size_t i, size_t j) const;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_SPARSEMATRIX_HPP