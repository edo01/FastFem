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

private:
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
    void add_entry(size_t index, double value);
    void print_pattern() const;
    inline bool is_symmetric() const override { return true; }
    void cholesky_factorize();
    Vector cholesky_solve(const Vector& b) const;
    void insert_entry(size_t i, size_t j, double value);
    //TO BE CANCELLED
    void set_values(const std::vector<double>& values);

private:
    const double &get_entry(size_t i, size_t j) const override;
};

template <unsigned int dim, unsigned int spacedim>
SkylinePattern SkylinePattern::create_from_dof_handler(const fastfem::dof::DoFHandler<dim, spacedim>& dof_handler)
{
    std::vector<std::set<unsigned int>> dof_interactions(dof_handler.get_n_dofs());
    
    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it) {
        const auto& elem = *it;

        const auto& dofs = dof_handler.get_ordered_dofs_on_element(elem);
        for(int j = 0; j < dofs.size(); ++j){
            for(int k = j; k < dofs.size(); ++k){
                dofs[j] > dofs[k] ? dof_interactions[dofs[j]].insert(dofs[k]) : dof_interactions[dofs[k]].insert(dofs[j]); 
            }
        }
    }

    std::vector<size_t> skyline_rows(dof_handler.get_n_dofs() + 1);

    for(unsigned int i = 0; i < dof_interactions.size(); ++i){
        skyline_rows[i + 1] = skyline_rows[i] + (i - *dof_interactions[i].begin() + 1);
    }

    return SkylinePattern(skyline_rows);
}

} // namespace linalg
} // namespace FastFem   

#endif // FASTFEM_SKYLINEMATRIX_HPP