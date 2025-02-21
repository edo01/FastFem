#ifndef FASTFEM_MATRIXTOOLS_HPP
#define FASTFEM_MATRIXTOOLS_HPP

#include "FastFem/linalg/sparseMatrices/SparseMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include "FastFem/linalg/Vector.hpp"
#include "FastFem/dof/DofHandler.hpp"
#include "FastFem/types/CommonTypes.hpp"

using fastfem::types::global_dof_index;
using fastfem::types::ff_index;
using fastfem::types::boundary_index;


namespace fastfem{
namespace linalg{

// Forward declarations needed for the friend declaration
class CSRMatrix;

class MatrixTools
{
public:
    MatrixTools() = delete;

    template <unsigned int dim, unsigned int spacedim>
    static void apply_homogeneous_dirichlet(SparseMatrix& A, Vector& rhs, const dof::DoFHandler<dim, spacedim> & dof_handler, boundary_index tag);

    template <unsigned int dim, unsigned int spacedim>
    static void apply_homogeneous_dirichlet(CSRMatrix& A, Vector& rhs, const dof::DoFHandler<dim, spacedim> & dof_handler, boundary_index tag);

    static void add_local_matrix_to_global(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index>& local_dofs);

    static void add_local_vector_to_global(Vector& global_vector, const Vector& local_vector, const std::vector<global_dof_index>& local_dofs);

    static void interpolate(Vector& f_interpolated, const dof::DoFHandler<1, 1>& dof_handler, const std::function<double(double)>& f);

    template <unsigned int dim>
    static void interpolate(Vector& f_interpolated, const dof::DoFHandler<dim, 2>& dof_handler, const std::function<double(double, double)>& f);
    
    template <unsigned int dim>
    static void interpolate(Vector& f_interpolated, const dof::DoFHandler<dim, 3>& dof_handler, const std::function<double(double, double, double)>& f);
};

} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_MATRIXTOOLS_HPP