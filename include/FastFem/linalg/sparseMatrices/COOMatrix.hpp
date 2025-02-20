#ifndef FASTFEM_COOMATRIX_HPP
#define FASTFEM_COOMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>
#include <FastFem/linalg/sparseMatrices/CSRMatrix.hpp>

namespace fastfem{
namespace linalg{

using types::ff_index;

class COOMatrix : public SparseMatrix {
private:
    std::vector<ff_index> row_indices;
    std::vector<ff_index> col_indices;
    std::vector<double> values;

public:
    COOMatrix(ff_index n_rows, ff_index n_cols, ff_index nnz_hint = 0);
    Vector gemv(const Vector& x) const override;
    inline ff_index nnz() const override { return row_indices.size(); }
    void set_entry(ff_index i, ff_index j, double value) override;
    inline void accumulate_entry(ff_index i, ff_index j, double value) override { set_entry(i, j, value); }
    void operator=(const double& value) override;

    static COOMatrix random(ff_index n_rows, ff_index n_cols, double sparsity, double min_value = -1.0, double max_value = 1.0, bool upper_triangular = false);

    CSRMatrix to_CSR() const;

    void print_pattern() const;
private:
    const double &get_entry(ff_index i, ff_index j) const override;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_COOMATRIX_HPP