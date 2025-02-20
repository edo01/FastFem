#include "FastFem/linalg/MatrixTools.hpp"

namespace fastfem{
namespace linalg{

template <unsigned int dim, unsigned int spacedim>
void MatrixTools::apply_homogeneous_dirichlet(SparseMatrix& A, Vector& rhs, const DoFHandler<dim, spacedim> & dof_handler, boundary_index tag)
{
    if(A.get_n_rows() != A.get_n_cols())
    {
        throw std::runtime_error("apply_homogeneous_dirichlet(): The matrix must be square");
    }

    if(A.get_n_rows() != dof_handler.get_n_dofs())
    {
        throw std::runtime_error("apply_homogeneous_dirichlet(): The matrix and the dof handler must have the same number of rows");
    }

    for(auto it = dof_handler.boundary_dofs_begin(tag); it != dof_handler.boundary_dofs_end(tag); ++it)
    {
        global_dof_index dof = *it;
        A.set_row_col_to_zero(dof);
        rhs[dof] = 0.0;

        //maybe we should set the diagonal entry to something different than one to keep the matrix well conditioned
        A.set_entry(dof, dof, 1.0);
    }
}

template <unsigned int dim, unsigned int spacedim>
void MatrixTools::apply_homogeneous_dirichlet(CSRMatrix& A, Vector& rhs, const DoFHandler<dim, spacedim> & dof_handler, boundary_index tag)
{
    if(A.get_n_rows() != A.get_n_cols())
    {
        throw std::runtime_error("apply_homogeneous_dirichlet(): The matrix must be square");
    }

    if(A.get_n_rows() != dof_handler.get_n_dofs())
    {
        throw std::runtime_error("apply_homogeneous_dirichlet(): The matrix and the dof handler must have the same number of rows");
    }

    // create a map that associates a column index to the indices of the values in the values array
    std::map<global_dof_index, std::vector<ff_index>> col_to_values;
    std::vector<ff_index> &col_indices = A.base_pattern->col_indices;

    for(ff_index i = 0; i < col_indices.size(); ++i)
    {
        col_to_values[col_indices[i]].push_back(i);
    }

    for(auto it = dof_handler.boundary_dofs_begin(tag); it != dof_handler.boundary_dofs_end(tag); ++it)
    {
        global_dof_index dof = *it;
        A.set_row_col_to_zero(dof, col_to_values);
        rhs[dof] = 0.0;

        //maybe we should set the diagonal entry to something different than one to keep the matrix well conditioned
        A.set_entry(dof, dof, 1.0);
    }
}

void MatrixTools::add_local_matrix_to_global(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index>& local_dofs)
{
    for(size_t i = 0; i < local_dofs.size(); ++i)
    {
        for(size_t j = 0; j < local_dofs.size(); ++j)
        {
            A.accumulate_entry(local_dofs[i], local_dofs[j], local_matrix(i, j));
        }
    }
}

void MatrixTools::add_local_vector_to_global(Vector& global_vector, const Vector& local_vector, const std::vector<global_dof_index>& local_dofs)
{
    for(size_t i = 0; i < local_dofs.size(); ++i)
    {
        global_vector[local_dofs[i]] += local_vector[i];
    }
}

// instantiate the templates

template void MatrixTools::apply_homogeneous_dirichlet<1,1>(SparseMatrix& A, Vector& rhs, const DoFHandler<1, 1>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<1,2>(SparseMatrix& A, Vector& rhs, const DoFHandler<1, 2>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<1,3>(SparseMatrix& A, Vector& rhs, const DoFHandler<1, 3>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<2,2>(SparseMatrix& A, Vector& rhs, const DoFHandler<2, 2>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<2,3>(SparseMatrix& A, Vector& rhs, const DoFHandler<2, 3>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<3,3>(SparseMatrix& A, Vector& rhs, const DoFHandler<3, 3>& dof_handler, boundary_index tag);

template void MatrixTools::apply_homogeneous_dirichlet<1,1>(CSRMatrix& A, Vector& rhs, const DoFHandler<1, 1>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<1,2>(CSRMatrix& A, Vector& rhs, const DoFHandler<1, 2>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<1,3>(CSRMatrix& A, Vector& rhs, const DoFHandler<1, 3>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<2,2>(CSRMatrix& A, Vector& rhs, const DoFHandler<2, 2>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<2,3>(CSRMatrix& A, Vector& rhs, const DoFHandler<2, 3>& dof_handler, boundary_index tag);
template void MatrixTools::apply_homogeneous_dirichlet<3,3>(CSRMatrix& A, Vector& rhs, const DoFHandler<3, 3>& dof_handler, boundary_index tag);

} // namespace linalg
} // namespace FastFem