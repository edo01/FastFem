#include <iostream>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"

#include "FastFem/fe/FESimplexP.hpp"
#include "FastFem/dof/DofHandler.hpp"

#include "FastFem/linalg/Vector.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
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

    //auto f_constant = [](double x, double y) { return 1; };
    // auto f = [](double x, double y) { return 4 - 2 * (x * x + y * y); };
    auto f = [](double x, double y) { return 10*10*10*10 * std::exp(- ((x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5))/0.001); };

    linalg::Vector rhs(n_dofs);

    linalg::CSRPattern csr_pattern = linalg::CSRPattern::create_from_dof_handler(dof_handler);

    linalg::CSRMatrix A(n_dofs, csr_pattern);

    linalg::FullMatrix local_matrix(n_dofs_per_cell);
    linalg::Vector local_rhs(n_dofs_per_cell);

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it)
    {
        auto &elem = *it;

        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);

        mesh::Point<2> v0 = triangle.get_vertex(0);
        mesh::Point<2> v1 = triangle.get_vertex(1);
        mesh::Point<2> v2 = triangle.get_vertex(2);

        mesh::Point<2> centroid = triangle.get_centroid();

        double volume = triangle.volume();

        double lenght_01 = distance(v0, v1);
        double lenght_12 = distance(v1, v2);
        double lenght_20 = distance(v2, v0);

        local_matrix.set_to_zero();
        local_rhs.fill(0.0);

        local_matrix(0, 0) += lenght_12 * lenght_12 / (4 * volume);
        local_matrix(1, 1) += lenght_20 * lenght_20 / (4 * volume);
        local_matrix(2, 2) += lenght_01 * lenght_01 / (4 * volume);

        local_matrix(0, 1) += dot(v1, v2, v2, v0) / (4 * volume);
        local_matrix(0, 2) += dot(v2, v1, v1, v0) / (4 * volume);
        local_matrix(1, 0) += dot(v0, v2, v2, v1) / (4 * volume);
        local_matrix(1, 2) += dot(v2, v0, v0, v1) / (4 * volume);
        local_matrix(2, 0) += dot(v0, v1, v1, v2) / (4 * volume);
        local_matrix(2, 1) += dot(v1, v0, v0, v2) / (4 * volume);
    
        local_rhs[0] += f(centroid.coords[0], centroid.coords[1]) * volume / 3;
        local_rhs[1] += f(centroid.coords[0], centroid.coords[1]) * volume / 3;
        local_rhs[2] += f(centroid.coords[0], centroid.coords[1]) * volume / 3;

        auto local_dofs = dof_handler.get_ordered_dofs_on_element(elem);
        linalg::MatrixTools::add_local_matrix_to_global(A, local_matrix, local_dofs);
        linalg::MatrixTools::add_local_vector_to_global(rhs, local_rhs, local_dofs);

    } 

    linalg::MatrixTools::apply_homogeneous_dirichlet(A, rhs, dof_handler, 0);

    linalg::CGSolver solver(1000, 1e-12);
    linalg::Vector sol = solver.solve(A, rhs);

    fastfem::mesh::DataIO<2, 2> data_io(mesh, dof_handler, sol);
    data_io.save_vtx("solution_csr.vtk");

    std::cout << "Max of solution: " << sol.max() << std::endl;

    auto exact_f = [](double x, double y) { return (1 - x * x) * (1 - y * y); };

    linalg::Vector exact_sol(n_dofs);

    // Questa può diventare una funzione di utility
    for(auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it)
    {
        auto &elem = *it;

        std::vector<fastfem::types::global_dof_index> global_dofs = dof_handler.get_ordered_dofs_on_element(elem);

        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);

        for(types::local_dof_index i = 0; i < global_dofs.size(); ++i)
        {
            mesh::Point<2> p_dof = fe.get_dof_coords(triangle, i);
            exact_sol[global_dofs[i]] = exact_f(p_dof.coords[0], p_dof.coords[1]);
        }
    }

    std::cout << "Norm of difference: " << (sol - exact_sol).norm() << std::endl;

    mesh::DataIO<2, 2> data_io_exact(mesh, dof_handler, exact_sol);
    data_io_exact.save_vtx("exact_solution.vtk");

    return EXIT_SUCCESS;
}