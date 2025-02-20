#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include <algorithm>
#include <random>

namespace fastfem{
namespace linalg{

using types::ff_index;

CSRPattern::CSRPattern(const std::vector<ff_index>& row_ptr, const std::vector<ff_index>& col_indices) : row_ptr(row_ptr), col_indices(col_indices)
{
    if(this->row_ptr.size() < 2){
        throw std::invalid_argument("CSRPattern::CSRPattern(): invalid row_ptr size");
    }
    if(this->col_indices.size() != (!this->row_ptr.empty() ? this->row_ptr.back() : 0)){
        throw std::invalid_argument("CSRPattern::CSRPattern(): invalid col_indices size");
    }
}

CSRMatrix::CSRMatrix(ff_index n_cols, const CSRPattern& pattern) :
  SparseMatrix(pattern.row_ptr.size() - 1, n_cols),
  values(pattern.col_indices.size()),
  base_pattern(std::make_shared<CSRPattern>(pattern))
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

const double &CSRMatrix::get_entry(ff_index i, ff_index j) const
{
    ff_index row_start = base_pattern->row_ptr[i];
    ff_index row_end = base_pattern->row_ptr[i + 1];

    for(ff_index k = row_start; k < row_end; ++k){
        if(base_pattern->col_indices[k] == j){
            return values[k];
        }
    }

    static double dummy = 0.0;
    return dummy;
}

void CSRMatrix::set_entry(ff_index i, ff_index j, double value)
{
    ff_index row_start = base_pattern->row_ptr[i];
    ff_index row_end = base_pattern->row_ptr[i + 1];

    for(ff_index k = row_start; k < row_end; ++k){
        if(base_pattern->col_indices[k] == j){
            values[k] = value;
            return;
        }
    }

    throw std::invalid_argument("CSRMatrix::insert_entry(): entry not found");
}

void CSRMatrix::accumulate_entry(ff_index i, ff_index j, double value)
{
    ff_index row_start = base_pattern->row_ptr[i];
    ff_index row_end = base_pattern->row_ptr[i + 1];

    for(ff_index k = row_start; k < row_end; ++k){
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

#ifdef HAVE_OPENMP
    #pragma omp parallel for
#endif
    for(ff_index i = 0; i < n_rows; ++i){
        ff_index row_start = row_ptr[i];
        ff_index row_end = row_ptr[i + 1];
        for(ff_index k = row_start; k < row_end; ++k){
            y[i] += values[k] * x[col_indices[k]];
        }
    }

    return y;
}

void CSRMatrix::set_row_col_to_zero(ff_index i)
{
    ff_index row_start = base_pattern->row_ptr[i];
    ff_index row_end = base_pattern->row_ptr[i + 1];

#ifdef HAVE_OPENMP
    #pragma omp parallel
#endif
    {
        // set row to zero
#ifdef HAVE_OPENMP
        #pragma omp for
#endif
        for(ff_index k = row_start; k < row_end; ++k){
            values[k] = 0.0;
        }

        // set column to zero
        auto& col_idx = base_pattern->col_indices;
        
#ifdef HAVE_OPENMP
        #pragma omp for
#endif
        for(ff_index j = 0; j < nnz(); ++j){
            if(col_idx[j] == i){
                values[j] = 0.0;
            }
        }
    }
}

// void CSRMatrix::set_row_col_to_zero(ff_index i, std::map<ff_index, std::vector<unsigned int>>& col_to_values)
// {
//     ff_index row_start = base_pattern->row_ptr[i];
//     ff_index row_end = base_pattern->row_ptr[i + 1];

//     for(ff_index k = row_start; k < row_end; ++k){
//         values[k] = 0.0;
//     }
    
//     if(col_to_values.find(i) != col_to_values.end()){
//         for(auto& j : col_to_values[i]){
//             values[j] = 0.0;
//         }
//     }
// }
        

void CSRMatrix::print_pattern() const
{
    for(ff_index i = 0; i < n_rows; ++i)
    {
        for(ff_index j = 0; j < n_cols; ++j)
        {
            bool found = false;
            for(ff_index k = base_pattern->row_ptr[i]; k < base_pattern->row_ptr[i + 1]; ++k)
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

        for(ff_index j = 0; j < (ff_index) dofs.size(); ++j){
            for(ff_index k = j; k < (ff_index) dofs.size(); ++k){
                dof_interactions[dofs[j]].insert(dofs[k]);
                dof_interactions[dofs[k]].insert(dofs[j]);
            }
        }
    }

    std::vector<ff_index> row_ptr(dof_handler.get_n_dofs() + 1);
    std::vector<ff_index> col_indices;

    for(ff_index i = 0; i < (ff_index) dof_interactions.size(); ++i){
        row_ptr[i + 1] = row_ptr[i] + dof_interactions[i].size();
        col_indices.insert(col_indices.end(), dof_interactions[i].begin(), dof_interactions[i].end());
    }

    return CSRPattern(row_ptr, col_indices);
}

template <unsigned int dim, unsigned int spacedim>
CSRPattern CSRPattern::create_symmetric_from_dof_handler(const fastfem::dof::DoFHandler<dim, spacedim>& dof_handler)
{
    std::vector<std::set<ff_index>> dof_interactions(dof_handler.get_n_dofs());

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it) {
        const auto& elem = *it;
        const auto& dofs = dof_handler.get_ordered_dofs_on_element(elem);
        for(ff_index j = 0; j < (ff_index) dofs.size(); ++j){
            for(ff_index k = j; k < (ff_index) dofs.size(); ++k){
                //stores only the lower triangular part
                dofs[j] < dofs[k] ? dof_interactions[dofs[j]].insert(dofs[k]) : dof_interactions[dofs[k]].insert(dofs[j]); 
            }
        }
    }

    std::vector<ff_index> row_ptr(dof_handler.get_n_dofs() + 1);
    std::vector<ff_index> col_indices;

    for(ff_index i = 0; i < dof_interactions.size(); ++i){
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
