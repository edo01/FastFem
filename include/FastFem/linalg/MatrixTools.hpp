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

class FullMatrix
{
public:
    FullMatrix(ff_index n_rows, ff_index n_cols) : n_rows(n_rows), n_cols(n_cols), data(n_rows * n_cols) {}
    FullMatrix(ff_index dim) : FullMatrix(dim, dim) {}
    inline const double &operator()(ff_index i, ff_index j) const { return data[i * n_cols + j]; }
    inline double &operator()(ff_index i, ff_index j) { return data[i * n_cols + j]; }

    inline void set_to_zero() { std::fill(data.begin(), data.end(), 0.0); }

    inline ff_index get_n_rows() const { return n_rows; }
    inline ff_index get_n_cols() const { return n_cols; }
private:
    ff_index n_rows;
    ff_index n_cols;
    std::vector<double> data;
};

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