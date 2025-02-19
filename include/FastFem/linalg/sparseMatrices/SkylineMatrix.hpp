#ifndef FASTFEM_SKYLINEMATRIX_HPP
#define FASTFEM_SKYLINEMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>
#include <FastFem/dof/DofHandler.hpp>
#include <vector>
#include <memory>
#include <set>

namespace fastfem{
namespace linalg{

struct SkylinePattern
{
    std::vector<size_t> skyline_rows;

// private:
public:
    SkylinePattern(const std::vector<size_t>& skyline) : skyline_rows(skyline) {}

public:
    template <unsigned int dim, unsigned int spacedim>
    static SkylinePattern create_from_dof_handler(const fastfem::dof::DoFHandler<dim, spacedim>& dof_handler);
};

/**
 * @brief Stores a symmetric matrix in skyline format, storing only the lower triangular part
 */
class SkylineMatrix : public SparseMatrix {
protected:
    std::vector<double> values;
    std::shared_ptr<SkylinePattern> base_skyline;

public:
    SkylineMatrix(size_t n_cols, const SkylinePattern& skyline);
    Vector gemv(const Vector& x) const override;
    inline size_t nnz() const override { return base_skyline->skyline_rows.back(); } 
    void print_pattern() const;
    inline bool is_symmetric() const override { return true; }
    void cholesky_factorize();
    Vector cholesky_solve(const Vector& b) const;
    void set_entry(size_t i, size_t j, double value) override;
    void accumulate_entry(size_t i, size_t j, double value) override;
    //TO BE CANCELLED
    void set_values(const std::vector<double>& values);

    void set_entry_to_zero(size_t i, size_t j);
    void set_row_col_to_zero(size_t i) override;

private:
    const double &get_entry(size_t i, size_t j) const override;
};


} // namespace linalg
} // namespace FastFem   

#endif // FASTFEM_SKYLINEMATRIX_HPP