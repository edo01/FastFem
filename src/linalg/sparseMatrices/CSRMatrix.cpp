#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include <algorithm>
#include <random>

namespace fastfem{
namespace linalg{

CSRPattern::CSRPattern(const std::vector<size_t>& row_ptr, const std::vector<size_t>& col_indices) : row_ptr(row_ptr), col_indices(col_indices)
{
    if(this->row_ptr.size() < 2){
        throw std::invalid_argument("CSRPattern::CSRPattern(): invalid row_ptr size");
    }
    if(this->col_indices.size() != (!this->row_ptr.empty() ? this->row_ptr.back() : 0)){
        throw std::invalid_argument("CSRPattern::CSRPattern(): invalid col_indices size");
    }
}

CSRMatrix::CSRMatrix(size_t n_cols, const CSRPattern& pattern) :
  SparseMatrix(pattern.row_ptr.size() - 1, n_cols),
  base_pattern(std::make_shared<CSRPattern>(pattern)),
  values(pattern.col_indices.size())
{
    if(nnz() > 0 && *std::max_element(pattern.col_indices.begin(), pattern.col_indices.end()) >= n_cols){
        throw std::invalid_argument("CSRMatrix::CSRMatrix(): invalid n_cols");
    }
}

CSRMatrix::CSRMatrix(const CSRMatrix& A) :
  SparseMatrix(A.n_rows, A.n_cols),
  values(A.values),
  base_pattern(A.base_pattern)
{}

const double &CSRMatrix::get_entry(size_t i, size_t j) const
{
    size_t row_start = base_pattern->row_ptr[i];
    size_t row_end = base_pattern->row_ptr[i + 1];

    for(size_t k = row_start; k < row_end; ++k){
        if(base_pattern->col_indices[k] == j){
            return values[k];
        }
    }

    static double dummy = 0.0;
    return dummy;
}

void CSRMatrix::set_entry(size_t i, size_t j, double value)
{
    size_t row_start = base_pattern->row_ptr[i];
    size_t row_end = base_pattern->row_ptr[i + 1];

    for(size_t k = row_start; k < row_end; ++k){
        if(base_pattern->col_indices[k] == j){
            values[k] = value;
            return;
        }
    }

    throw std::invalid_argument("CSRMatrix::insert_entry(): entry not found");
}

void CSRMatrix::accumulate_entry(size_t i, size_t j, double value)
{
    size_t row_start = base_pattern->row_ptr[i];
    size_t row_end = base_pattern->row_ptr[i + 1];

    for(size_t k = row_start; k < row_end; ++k){
        if(base_pattern->col_indices[k] == j){
            values[k] += value;
            return;
        }
    }

    throw std::invalid_argument("CSRMatrix::accumulate_entry(): entry not found");
}

Vector CSRMatrix::gemv(const Vector& x) const {

    if(n_cols != x.size()){
        throw std::invalid_argument("CSRMatrix::gemv(): incompatible dimensions");
    }

    Vector y(n_rows);
    const auto& row_ptr = base_pattern->row_ptr;
    const auto& col_indices = base_pattern->col_indices;

    #pragma omp parallel for
    for(size_t i = 0; i < n_rows; ++i){
        size_t row_start = row_ptr[i];
        size_t row_end = row_ptr[i + 1];
        for(size_t k = row_start; k < row_end; ++k){
            y[i] += values[k] * x[col_indices[k]];
        }
    }

    return y;
}

void CSRMatrix::set_row_col_to_zero(size_t i)
{
    size_t row_start = base_pattern->row_ptr[i];
    size_t row_end = base_pattern->row_ptr[i + 1];

    #pragma omp parallel
    {
        // set row to zero
        #pragma omp for
        for(size_t k = row_start; k < row_end; ++k){
            values[k] = 0.0;
        }

        // set column to zero
        auto& col_idx = base_pattern->col_indices;
        #pragma omp for
        for(size_t j = 0; j < nnz(); ++j){
            if(col_idx[j] == i){
                values[j] = 0.0;
            }
        }
    }
}

void CSRMatrix::set_row_col_to_zero(size_t i, std::map<size_t, std::vector<unsigned int>>& col_to_values)
{
    size_t row_start = base_pattern->row_ptr[i];
    size_t row_end = base_pattern->row_ptr[i + 1];

    for(size_t k = row_start; k < row_end; ++k){
        values[k] = 0.0;
    }
    
    if(col_to_values.find(i) != col_to_values.end()){
        for(auto& j : col_to_values[i]){
            values[j] = 0.0;
        }
    }
}
        

void CSRMatrix::print_pattern() const
{
    for(size_t i = 0; i < n_rows; ++i)
    {
        for(size_t j = 0; j < n_cols; ++j)
        {
            bool found = false;
            for(size_t k = base_pattern->row_ptr[i]; k < base_pattern->row_ptr[i + 1]; ++k)
            {
                if(base_pattern->col_indices[k] == j)
                {
                    found = true;
                    break;
                }
            }
            std::cout << (found ? "x" : ".") << " ";
        }
        std::cout << std::endl;
    }
}        

/**
 * CSRPattern
 */

template <unsigned int dim, unsigned int spacedim>
CSRPattern CSRPattern::create_from_dof_handler(const fastfem::dof::DoFHandler<dim, spacedim>& dof_handler)
{
    std::vector<std::set<unsigned int>> dof_interactions(dof_handler.get_n_dofs());

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it) {
        const auto& elem = *it;

        const auto& dofs = dof_handler.get_ordered_dofs_on_element(elem);

        for(int j = 0; j < dofs.size(); ++j){
            for(int k = j; k < dofs.size(); ++k){
                dof_interactions[dofs[j]].insert(dofs[k]);
                dof_interactions[dofs[k]].insert(dofs[j]);
            }
        }
    }

    std::vector<size_t> row_ptr(dof_handler.get_n_dofs() + 1);
    std::vector<size_t> col_indices;

    for(unsigned int i = 0; i < dof_interactions.size(); ++i){
        row_ptr[i + 1] = row_ptr[i] + dof_interactions[i].size();
        col_indices.insert(col_indices.end(), dof_interactions[i].begin(), dof_interactions[i].end());
    }

    return CSRPattern(row_ptr, col_indices);
}

template <unsigned int dim, unsigned int spacedim>
CSRPattern CSRPattern::create_symmetric_from_dof_handler(const fastfem::dof::DoFHandler<dim, spacedim>& dof_handler)
{
    std::vector<std::set<unsigned int>> dof_interactions(dof_handler.get_n_dofs());

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it) {
        const auto& elem = *it;
        const auto& dofs = dof_handler.get_ordered_dofs_on_element(elem);
        for(int j = 0; j < dofs.size(); ++j){
            for(int k = j; k < dofs.size(); ++k){
                //stores only the lower triangular part
                dofs[j] < dofs[k] ? dof_interactions[dofs[j]].insert(dofs[k]) : dof_interactions[dofs[k]].insert(dofs[j]); 
            }
        }
    }

    std::vector<size_t> row_ptr(dof_handler.get_n_dofs() + 1);
    std::vector<size_t> col_indices;

    for(unsigned int i = 0; i < dof_interactions.size(); ++i){
        row_ptr[i + 1] = row_ptr[i] + dof_interactions[i].size();
        col_indices.insert(col_indices.end(), dof_interactions[i].begin(), dof_interactions[i].end());
    }

    return CSRPattern(row_ptr, col_indices);
}

// explicit instantiation
template CSRPattern CSRPattern::create_from_dof_handler<1,1>(const fastfem::dof::DoFHandler<1,1>& dof_handler);
template CSRPattern CSRPattern::create_from_dof_handler<1,2>(const fastfem::dof::DoFHandler<1,2>& dof_handler);
template CSRPattern CSRPattern::create_from_dof_handler<1,3>(const fastfem::dof::DoFHandler<1,3>& dof_handler);
template CSRPattern CSRPattern::create_from_dof_handler<2,2>(const fastfem::dof::DoFHandler<2,2>& dof_handler);
template CSRPattern CSRPattern::create_from_dof_handler<2,3>(const fastfem::dof::DoFHandler<2,3>& dof_handler);
template CSRPattern CSRPattern::create_from_dof_handler<3,3>(const fastfem::dof::DoFHandler<3,3>& dof_handler);

template CSRPattern CSRPattern::create_symmetric_from_dof_handler<1,1>(const fastfem::dof::DoFHandler<1,1>& dof_handler);
template CSRPattern CSRPattern::create_symmetric_from_dof_handler<1,2>(const fastfem::dof::DoFHandler<1,2>& dof_handler);
template CSRPattern CSRPattern::create_symmetric_from_dof_handler<1,3>(const fastfem::dof::DoFHandler<1,3>& dof_handler);
template CSRPattern CSRPattern::create_symmetric_from_dof_handler<2,2>(const fastfem::dof::DoFHandler<2,2>& dof_handler);
template CSRPattern CSRPattern::create_symmetric_from_dof_handler<2,3>(const fastfem::dof::DoFHandler<2,3>& dof_handler);
template CSRPattern CSRPattern::create_symmetric_from_dof_handler<3,3>(const fastfem::dof::DoFHandler<3,3>& dof_handler);

} // namespace linalg
} // namespace FastFem
