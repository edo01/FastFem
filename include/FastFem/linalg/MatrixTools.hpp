#ifndef FASTFEM_MATRIXTOOLS_HPP
#define FASTFEM_MATRIXTOOLS_HPP

#include "FastFem/linalg/sparseMatrices/SparseMatrix.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include "FastFem/linalg/Vector.hpp"
#include "FastFem/dof/DofHandler.hpp"

using namespace fastfem::dof;
using fastfem::types::global_dof_index;

namespace fastfem{
namespace linalg{
namespace tools{

class FullMatrix
{
public:
    FullMatrix(size_t n_rows, size_t n_cols) : n_rows(n_rows), n_cols(n_cols), data(n_rows * n_cols) {}
    FullMatrix(size_t dim) : FullMatrix(dim, dim) {}
    inline const double &operator()(size_t i, size_t j) const { return data[i * n_cols + j]; }
    inline double &operator()(size_t i, size_t j) { return data[i * n_cols + j]; }

    inline void set_to_zero() { std::fill(data.begin(), data.end(), 0.0); }

private:
    size_t n_rows;
    size_t n_cols;
    std::vector<double> data;
};

template <unsigned int dim, unsigned int spacedim>
void apply_homogeneous_dirichlet(SparseMatrix& A, Vector& rhs, const DoFHandler<dim, spacedim> & dof_handler, size_t tag);

// template <unsigned int dim, unsigned int spacedim>
// void apply_homogeneous_dirichlet(CSRMatrix& A, Vector& rhs, const DoFHandler<dim, spacedim> & dof_handler, size_t tag);

void add_local_matrix_to_global(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index>& local_dofs);

void add_local_vector_to_global(Vector& global_vector, const Vector& local_vector, const std::vector<global_dof_index>& local_dofs);

} // namespace tools
} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_MATRIXTOOLS_HPP