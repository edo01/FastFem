#ifndef FASTFEM_MATRIXTOOLS_HPP
#define FASTFEM_MATRIXTOOLS_HPP

#include "FastFem/linalg/sparseMatrices/SparseMatrix.hpp"
#include "FastFem/linalg/Vector.hpp"
#include "FastFem/dof/DofHandler.hpp"

using namespace fastfem::dof;
using fastfem::types::global_dof_index_t;

namespace fastfem{
namespace linalg{
namespace matrixtools{

class FullMatrix
{
public:
    FullMatrix(size_t n_rows, size_t n_cols) : n_rows(n_rows), n_cols(n_cols), data(n_rows * n_cols) {}
    FullMatrix(size_t dim) : FullMatrix(dim, dim) {}
    inline const double &operator()(size_t i, size_t j) const { return data[i * n_cols + j]; }
    inline double &operator()(size_t i, size_t j) { return data[i * n_cols + j]; }

    inline void set_to_zero() { std::fill(data.begin(), data.end(), 0.0); }

private:
    std::vector<double> data;
    size_t n_rows;
    size_t n_cols;
};

template <unsigned int dim, unsigned int spacedim>
void apply_homogeneous_dirichlet(SparseMatrix& A, Vector& rhs, const DoFHandler<dim> & dof_handler, size_t tag);

template <unsigned int dim, unsigned int spacedim>
void add_local_matrix_to_global(SparseMatrix& A, const FullMatrix& local_matrix, const std::vector<global_dof_index_t>& local_dofs);

} // namespace matrixtools
} // namespace linalg
} // namespace FastFem

#endif // FASTFEM_MATRIXTOOLS_HPP