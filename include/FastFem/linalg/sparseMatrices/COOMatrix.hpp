#ifndef FASTFEM_COOMATRIX_HPP
#define FASTFEM_COOMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>
#include <FastFem/linalg/sparseMatrices/CSRMatrix.hpp>

namespace fastfem{
namespace linalg{

class COOMatrix : public SparseMatrix {
private:
    std::vector<size_t> row_indices;
    std::vector<size_t> col_indices;
    std::vector<double> values;

public:
    COOMatrix(size_t n_rows, size_t n_cols, size_t nnz_hint = 0);
    Vector gemv(const Vector& x) const override;
    inline size_t nnz() const override { return row_indices.size(); }
    void set_entry(size_t i, size_t j, double value) override;
    inline void accumulate_entry(size_t i, size_t j, double value) override { set_entry(i, j, value); }

    static COOMatrix random(size_t n_rows, size_t n_cols, double sparsity, double min_value = -1.0, double max_value = 1.0, bool upper_triangular = false);

    CSRMatrix to_CSR() const;

    void print_pattern() const;
private:
    const double &get_entry(size_t i, size_t j) const override;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_COOMATRIX_HPP