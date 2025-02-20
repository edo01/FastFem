#ifndef FASTFEM_SKYLINEMATRIX_HPP
#define FASTFEM_SKYLINEMATRIX_HPP

#include <FastFem/linalg/sparseMatrices/SparseMatrix.hpp>
#include <FastFem/dof/DofHandler.hpp>
#include <vector>
#include <memory>
#include <set>

namespace fastfem{
namespace linalg{

using types::ff_index;

struct SkylinePattern
{
    std::vector<ff_index> skyline_rows;

// private:
public:
    SkylinePattern(const std::vector<ff_index>& skyline) : skyline_rows(skyline) {}

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
    SkylineMatrix(ff_index n_cols, const SkylinePattern& skyline);
    Vector gemv(const Vector& x) const override;
    inline ff_index nnz() const override { return base_skyline->skyline_rows.back(); } 
    void print_pattern(bool values_flag) const;
    inline bool is_symmetric() const override { return true; }

    /**
     * @brief In place factorizatoin of the matrix using the Cholesky method
     */
    void cholesky_factorize();

    /**
     * @brief Solves the linear system Ax = b exploiting a previous Cholesky factorization
     * and performs the forward and backward substitutions
     * */
    Vector cholesky_solve(const Vector& b) const;

    void set_entry(ff_index i, ff_index j, double value) override;
    void accumulate_entry(ff_index i, ff_index j, double value) override;
    void set_entry_to_zero(ff_index i, ff_index j);

    void set_row_col_to_zero(ff_index i) override;

private:
    const double &get_entry(ff_index i, ff_index j) const override;
};


} // namespace linalg
} // namespace FastFem   

#endif // FASTFEM_SKYLINEMATRIX_HPP