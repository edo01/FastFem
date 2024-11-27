#ifndef POISSON_HPP
#define POISSON_HPP

#include <stdlib.h>
#include "mesh/mesh.hpp"
#include "linalg/sparse_matrix.hpp"
#include "linalg/system.hpp"
#include "linalg/vector.hpp"

namespace fem::poisson
{
    /******************************************************************************
    * Builds the P1 stiffness and mass matrices of a given mesh.
    * We do not try to assemble different elements together here for simplicity.
    * Both matrices M and S will therefore have 9 * number of triangles.
    */
    void build_fem_matrices(const struct mesh::Mesh *m, struct linalg::SparseMatrix *S,
			struct linalg::SparseMatrix *M);

}

#endif // POISSON_HPP