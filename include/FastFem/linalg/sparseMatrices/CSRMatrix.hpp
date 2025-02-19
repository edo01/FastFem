#ifndef FASTFEM_CSRMATRIX_HPP
#define FASTFEM_CSRMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>
#include <FastFem/dof/DofHandler.hpp>
#include <vector>
#include <memory>
#include <set>

namespace fastfem{
namespace linalg{

struct CSRPattern
{
    std::vector<size_t> row_ptr;
    std::vector<size_t> col_indices;

private:
    //CSRPattern(std::initializer_list<size_t> row_ptr, std::initializer_list<size_t> col_indices);
    CSRPattern(const std::vector<size_t>& row_ptr, const std::vector<size_t>& col_indices);

public:
    friend class COOMatrix;

    template <unsigned int dim, unsigned int spacedim>
    static CSRPattern create_from_dof_handler(const fastfem::dof::DoFHandler<dim, spacedim>& dof_handler);

    template <unsigned int dim, unsigned int spacedim>
    static CSRPattern create_symmetric_from_dof_handler(const fastfem::dof::DoFHandler<dim, spacedim>& dof_handler);
};

class CSRMatrix : public SparseMatrix
{
protected:
    std::vector<double> values;
    std::shared_ptr<CSRPattern> base_pattern;

public:
    CSRMatrix(size_t n_cols, const CSRPattern& pattern);

    CSRMatrix(const CSRMatrix& A);

    Vector gemv(const Vector& x) const override;
    void set_entry(size_t i, size_t j, double value) override;
    void accumulate_entry(size_t i, size_t j, double value) override;

    inline size_t nnz() const override { return base_pattern->col_indices.size(); }

    void print_pattern() const;

    void set_row_col_to_zero(size_t i) override;
private:
    const double &get_entry(size_t i, size_t j) const override;

    friend class COOMatrix;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_CSRMATRIX_HPP
