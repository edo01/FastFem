#include <iostream>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"

#include "FastFem/fe/FESimplexP.hpp"
#include "FastFem/dof/DofHandler.hpp"

#include "FastFem/linalg/Vector.hpp"
#include "FastFem/linalg/sparseMatrices/SkylineMatrix.hpp"
#include "FastFem/linalg/MatrixTools.hpp"

#include "FastFem/linalg/iterativeSolvers/CGSolver.hpp"

#include "FastFem/mesh/MeshIO.hpp"

#include "FastFem/types/CommonTypes.hpp"

#include "FastFem/mesh/MeshIO.hpp"

using namespace fastfem;

double distance(const mesh::Point<2> &v1, const mesh::Point<2> &v2)
{
    return std::sqrt((v2.coords[0] - v1.coords[0]) * (v2.coords[0] - v1.coords[0]) + (v2.coords[1] - v1.coords[1]) * (v2.coords[1] - v1.coords[1]));
}

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
    
    mesh::SquareMaker mesh_maker(N);
    mesh::Mesh<2> mesh = mesh_maker.make_mesh();

    fe::FESimplexP1<2> fe;
    dof::DoFHandler<2> dof_handler(mesh, std::make_unique<fe::FESimplexP1<2>>());

    dof_handler.distribute_dofs();

    unsigned int n_dofs = dof_handler.get_n_dofs();
    unsigned int n_dofs_per_cell = fe.get_n_dofs_per_element();

    linalg::Vector rhs(n_dofs);

    linalg::SkylinePattern skyline_pattern = linalg::SkylinePattern::create_from_dof_handler(dof_handler);

    linalg::SkylineMatrix A_sky(n_dofs, skyline_pattern);

    linalg::FullMatrix local_matrix(n_dofs_per_cell);

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it)
    {
        auto &elem = *it;

        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);

        double volume = triangle.volume();

        mesh::Point<2> v0 = triangle.get_vertex(0);
        mesh::Point<2> v1 = triangle.get_vertex(1);
        mesh::Point<2> v2 = triangle.get_vertex(2);

        double lenght_01 = distance(v0, v1);
        double lenght_12 = distance(v1, v2);
        double lenght_20 = distance(v2, v0);

        local_matrix.set_to_zero();

        local_matrix(0, 0) += lenght_12 * lenght_12 / (4 * volume);
        local_matrix(1, 1) += lenght_20 * lenght_20 / (4 * volume);
        local_matrix(2, 2) += lenght_01 * lenght_01 / (4 * volume);

        local_matrix(0, 1) += dot(v1, v2, v2, v0) / (4 * volume);
        local_matrix(0, 2) += dot(v2, v1, v1, v0) / (4 * volume);
        local_matrix(1, 2) += dot(v2, v0, v0, v1) / (4 * volume);
    
        auto local_dofs = dof_handler.get_ordered_dofs_on_element(elem);
        linalg::MatrixTools::add_local_matrix_to_global(A_sky, local_matrix, local_dofs);
    } 

    linalg::MatrixTools::apply_homogeneous_dirichlet(A_sky, rhs, dof_handler, 0);

    rhs[1443] = 2;
    rhs[367] = 2;
    rhs[1900] = -2;

    linalg::CGSolver solver;
    linalg::Vector sol = solver.solve(A_sky, rhs);

    fastfem::mesh::DataIO<2, 2> data_io(mesh, dof_handler, sol);
    data_io.save_vtx("solution_skyline_cg.vtk");

    A_sky.cholesky_factorize();
    sol = A_sky.cholesky_solve(rhs);

    fastfem::mesh::DataIO<2, 2> data_io_cholesky(mesh, dof_handler, sol);
    data_io_cholesky.save_vtx("solution_skyline_cholesky.vtk");

    return EXIT_SUCCESS;
}