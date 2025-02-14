#ifndef CSRMATRIX_HPP
#define CSRMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>
#include <vector>
#include <memory>

namespace fastfem{
namespace linalg{

struct CSRPattern
{
    std::vector<size_t> row_ptr;
    std::vector<size_t> col_indices;

    CSRPattern(std::initializer_list<size_t> row_ptr, std::initializer_list<size_t> col_indices);
};

class CSRMatrix : public SparseMatrix
{
private:
    std::shared_ptr<CSRPattern> base_pattern;
    std::vector<double> values;

public:
    CSRMatrix(size_t n_cols, const CSRPattern& pattern);

    Vector gemv(const Vector& x) const override;

    inline size_t nnz() const override { return base_pattern->col_indices.size(); }

    void add_entry(size_t index, double value);

private:
    const double &get_entry(size_t i, size_t j) const override;
};

} // namespace linalg
} // namespace FastFem

#endif // CSRMATRIX_HPP