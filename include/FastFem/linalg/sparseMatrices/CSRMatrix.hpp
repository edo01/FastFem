#ifndef FASTFEM_CSRMATRIX_HPP
#define FASTFEM_CSRMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>
#include <FastFem/dof/DofHandler.hpp>
#include <FastFem/linalg/MatrixTools.hpp>
#include <vector>
#include <memory>
#include <set>
#include <map>

namespace fastfem{
namespace linalg{

struct CSRPattern
{
    std::vector<size_t> row_ptr;
    std::vector<size_t> col_indices;

private:
    CSRPattern(const std::vector<size_t>& row_ptr, const std::vector<size_t>& col_indices);

public:
    friend class COOMatrix;

    /**
     * @brief Static factory method to create a CSR pattern from a DoFHandler
     */
    template <unsigned int dim, unsigned int spacedim>
    static CSRPattern create_from_dof_handler(const fastfem::dof::DoFHandler<dim, spacedim>& dof_handler);

    /**
     * @brief Static factory method to create a symmetric CSR pattern, storing only the upper triangular part, from a DoFHandler
     */
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

    /**
     * @brief Optimized overloaded method to set a row and a column to zero,
     * exploiting an ordered map that stores the column indices for each row and improves the setting of the column to zero.
     */
    void set_row_col_to_zero(size_t i, std::map<size_t, std::vector<unsigned int>>& col_to_values);

    private:
    const double &get_entry(size_t i, size_t j) const override;

    friend class COOMatrix;
    friend class MatrixTools;
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_CSRMATRIX_HPP
