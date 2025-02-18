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
    static CSRPattern create_from_dof_handler(const fastfem::dof::DofHandler<dim, spacedim>& dof_handler);
    template <unsigned int dim, unsigned int spacedim>
    static CSRPattern create_symmetric_from_dof_handler(const fastfem::dof::DofHandler<dim, spacedim>& dof_handler);
};

class CSRMatrix : public SparseMatrix
{
protected:
    std::vector<double> values;
    std::shared_ptr<CSRPattern> base_pattern;

public:
    CSRMatrix(size_t n_cols, const CSRPattern& pattern);

    Vector gemv(const Vector& x) const override;
    void add_entry(size_t index, double value);
    void insert_entry(size_t i, size_t j, double value);
    void accumulate_entry(size_t i, size_t j, double value);

    inline size_t nnz() const override { return base_pattern->col_indices.size(); }
    void print_pattern() const;


private:
    const double &get_entry(size_t i, size_t j) const override;
};

template <unsigned int dim, unsigned int spacedim>
CSRPattern CSRPattern::create_from_dof_handler(const fastfem::dof::DofHandler<dim, spacedim>& dof_handler)
{
    ////////////////
    //std::vector<std::set<unsigned int>> dof_interactions(dof_handler.n_dofs());
    std::vector<std::set<unsigned int>> dof_interactions(4);

    for(unsigned int i = 0; i < dof_handler.n_elements(); ++i){
        const auto& dofs = dof_handler.get_element_dofs(i);
        for(int j = 0; j < dofs.size(); ++j){
            for(int k = j; k < dofs.size(); ++k){
                dof_interactions[dofs[j]].insert(dofs[k]);
                dof_interactions[dofs[k]].insert(dofs[j]);
            }
        }
    }

    //////////
    //std::vector<size_t> row_ptr(dof_handler.n_dofs() + 1);
    std::vector<size_t> row_ptr(5);
    std::vector<size_t> col_indices;

    for(unsigned int i = 0; i < dof_interactions.size(); ++i){
        row_ptr[i + 1] = row_ptr[i] + dof_interactions[i].size();
        col_indices.insert(col_indices.end(), dof_interactions[i].begin(), dof_interactions[i].end());
    }

    return CSRPattern(row_ptr, col_indices);
}

template <unsigned int dim, unsigned int spacedim>
CSRPattern CSRPattern::create_symmetric_from_dof_handler(const fastfem::dof::DofHandler<dim, spacedim>& dof_handler)
{
    ///////n_dofs
    //std::vector<std::set<unsigned int>> dof_interactions(dof_handler.n_dofs());
    std::vector<std::set<unsigned int>> dof_interactions(4);


    for(unsigned int i = 0; i < dof_handler.n_elements(); ++i){
        const auto& dofs = dof_handler.get_element_dofs(i);
        for(int j = 0; j < dofs.size(); ++j){
            for(int k = j; k < dofs.size(); ++k){
                //stores only the lower triangular part
                dofs[j] < dofs[k] ? dof_interactions[dofs[j]].insert(dofs[k]) : dof_interactions[dofs[k]].insert(dofs[j]); 
            }
        }
    }

    /////
    //std::vector<size_t> row_ptr(dof_handler.n_dofs + 1);
    std::vector<size_t> row_ptr(5);
    std::vector<size_t> col_indices;

    for(unsigned int i = 0; i < dof_interactions.size(); ++i){
        row_ptr[i + 1] = row_ptr[i] + dof_interactions[i].size();
        col_indices.insert(col_indices.end(), dof_interactions[i].begin(), dof_interactions[i].end());
    }

    return CSRPattern(row_ptr, col_indices);
}

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_CSRMATRIX_HPP


// 0 1 0 2     0 1 0 2
// 1 1 0 0       1 0 0
// 0 0 1 0 -->     1 0
// 2 0 0 2           2

// row_ptr = {0, 2, 3, 4, 5}
// col_indices = {1, 3, 1, 2, 3}

