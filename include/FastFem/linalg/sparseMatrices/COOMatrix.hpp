#ifndef FASTFEM_LINALG_SPARSEMATRICES_COOMATRIX_HPP
#define FASTFEM_LINALG_SPARSEMATRICES_COOMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>

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
    void add_entry(size_t i, size_t j, double value);
private:
    const double &get_entry(size_t i, size_t j) const override;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_LINALG_SPARSEMATRICES_COOMATRIX_HPP