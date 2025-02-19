#include "FastFem/linalg/sparseMatrices/MatrixTools.hpp"

namespace fastfem{
namespace linalg{
namespace matrixtools{

template <unsigned int dim, unsigned int spacedim>
void apply_homogeneous_dirichlet(SparseMatrix& A, Vector& rhs, const DoFHandler<dim> & dof_handler, size_t tag)
{
    if(A.get_n_rows() != A.get_n_cols())
    {
        throw std::runtime_error("The matrix must be square");
    }

    if(A.get_n_rows() != dof_handler.get_n_dofs())
    {
        throw std::runtime_error("The matrix and the DoF handler must have the same number of DoFs");
    }

    for(auto it = dof_handler.boundary_dofs_begin(tag); it != dof_handler.boundary_dofs_end(tag); ++it)
    {
        global_dof_index_t dof = *it;
        A.set_row_col_to_zero(dof);
        rhs[dof] = 0.0;

        A.set_entry(dof, dof, 1.0);
    }
}

template <unsigned int dim, unsigned int spacedim>
void add_local_matrix_to_global(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index_t>& local_dofs)
{
    for(size_t i = 0; i < local_dofs.size(); ++i)
    {
        for(size_t j = 0; j < local_dofs.size(); ++j)
        {
            A.accumulate_entry(local_dofs[i], local_dofs[j], local_matrix(i, j));
        }
    }
}

// instantiate the templates

template void apply_homogeneous_dirichlet<1,1>(SparseMatrix& A, Vector& rhs, const DoFHandler<1>& dof_handler, size_t tag);
template void apply_homogeneous_dirichlet<1,2>(SparseMatrix& A, Vector& rhs, const DoFHandler<1>& dof_handler, size_t tag);
template void apply_homogeneous_dirichlet<1,3>(SparseMatrix& A, Vector& rhs, const DoFHandler<1>& dof_handler, size_t tag);
template void apply_homogeneous_dirichlet<2,2>(SparseMatrix& A, Vector& rhs, const DoFHandler<2>& dof_handler, size_t tag);
template void apply_homogeneous_dirichlet<2,3>(SparseMatrix& A, Vector& rhs, const DoFHandler<2>& dof_handler, size_t tag);
template void apply_homogeneous_dirichlet<3,3>(SparseMatrix& A, Vector& rhs, const DoFHandler<3>& dof_handler, size_t tag);

template void add_local_matrix_to_global<1,1>(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index_t>& local_dofs);
template void add_local_matrix_to_global<1,2>(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index_t>& local_dofs);
template void add_local_matrix_to_global<1,3>(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index_t>& local_dofs);
template void add_local_matrix_to_global<2,2>(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index_t>& local_dofs);
template void add_local_matrix_to_global<2,3>(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index_t>& local_dofs);
template void add_local_matrix_to_global<3,3>(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index_t>& local_dofs);


} // namespace matrixtools
} // namespace linalg
} // namespace FastFem