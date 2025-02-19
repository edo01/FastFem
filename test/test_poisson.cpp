#include <iostream>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"

#include "FastFem/fe/FESimplexP.hpp"
#include "FastFem/dof/DofHandler.hpp"

#include "FastFem/linalg/Vector.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include "FastFem/linalg/MatrixTools.hpp"

#include "FastFem/linalg/iterativeSolvers/CGSolver.hpp"

#include "FastFem/types/CommonTypes.hpp"

#include "FastFem/mesh/MeshIO.hpp"

using namespace fastfem;

double dot(const mesh::Point<2> &v1, const mesh::Point<2> &v2, const mesh::Point<2> &w1, const mesh::Point<2> &w2)
{
    return (v2.coords[0] - v1.coords[0]) * (w2.coords[0] - w1.coords[0]) + (v2.coords[1] - v1.coords[1]) * (w2.coords[1] - w1.coords[1]);
}

#include <cstdlib>

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <N>" << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int N = std::atoi(argv[1]);
    unsigned int dim = 2;

    mesh::SquareMaker mesh_maker(N);
    mesh::Mesh<2> mesh = mesh_maker.make_mesh();

    fe::FESimplexP1<2> fe;
    dof::DoFHandler<2> dof_handler(mesh, std::make_unique<fe::FESimplexP1<2>>());

    dof_handler.distribute_dofs();

    unsigned int n_dofs = dof_handler.get_n_dofs();
    unsigned int n_dofs_per_cell = fe.get_n_dofs_per_element();

    linalg::Vector rhs(n_dofs);

    linalg::CSRPattern csr_pattern = linalg::CSRPattern::create_from_dof_handler(dof_handler);

    linalg::CSRMatrix A(n_dofs, csr_pattern);

    linalg::tools::FullMatrix local_matrix(n_dofs_per_cell);

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it)
    {
        auto &elem = *it;

        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);

        double volume = triangle.volume();

        double edge_01 = mesh::Simplex<1, 2>(triangle.get_vertex(0), triangle.get_vertex(1)).volume();
        double edge_12 = mesh::Simplex<1, 2>(triangle.get_vertex(1), triangle.get_vertex(2)).volume();
        double edge_20 = mesh::Simplex<1, 2>(triangle.get_vertex(2), triangle.get_vertex(0)).volume();

        mesh::Point<2> v0 = triangle.get_vertex(0);
        mesh::Point<2> v1 = triangle.get_vertex(1);
        mesh::Point<2> v2 = triangle.get_vertex(2);

        auto local_dofs = dof_handler.get_ordered_dofs_on_element(elem);

        local_matrix.set_to_zero();

        local_matrix(0, 0) += edge_12 * edge_12 / (4 * volume);
        local_matrix(1, 1) += edge_20 * edge_20 / (4 * volume);
        local_matrix(2, 2) += edge_01 * edge_01 / (4 * volume);

        local_matrix(0, 1) += dot(v1, v2, v2, v0) / (4 * volume);
        local_matrix(0, 2) += dot(v2, v1, v1, v0) / (4 * volume);
        local_matrix(1, 0) += dot(v0, v2, v2, v1) / (4 * volume);
        local_matrix(1, 2) += dot(v2, v0, v0, v1) / (4 * volume);
        local_matrix(2, 0) += dot(v0, v1, v1, v2) / (4 * volume);
        local_matrix(2, 1) += dot(v1, v0, v0, v2) / (4 * volume);

        linalg::tools::add_local_matrix_to_global(A, local_matrix, local_dofs);
    }

    for (unsigned int i = 0; i < n_dofs; ++i)
    {
        rhs[i] = 1.0;
    } 

    linalg::tools::apply_homogeneous_dirichlet(A, rhs, dof_handler, 0);

    linalg::CGSolver solver;
    linalg::Vector sol = solver.solve(A, rhs);

    linalg::Vector Av = A * sol;

    fastfem::mesh::DataIO<2, 2> data_io(mesh, dof_handler, sol);
    data_io.save_vtx("solution.vtk");

    return EXIT_SUCCESS;
}