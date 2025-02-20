#include <iostream>
#include <chrono>
#include <cstdlib>

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

using namespace fastfem;

double distance(const mesh::Point<2> &v1, const mesh::Point<2> &v2)
{
    return std::sqrt((v2.coords[0] - v1.coords[0]) * (v2.coords[0] - v1.coords[0]) + 
                     (v2.coords[1] - v1.coords[1]) * (v2.coords[1] - v1.coords[1]));
}

double dot(const mesh::Point<2> &v1, const mesh::Point<2> &v2, const mesh::Point<2> &w1, const mesh::Point<2> &w2)
{
    return (v2.coords[0] - v1.coords[0]) * (w2.coords[0] - w1.coords[0]) + 
           (v2.coords[1] - v1.coords[1]) * (w2.coords[1] - w1.coords[1]);
}

auto f = [](double x, double y) { return 4 - 2 * (x * x + y * y); };

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <N>" << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int N = std::atoi(argv[1]);

    auto start_total = std::chrono::high_resolution_clock::now();

    auto start_mesh = std::chrono::high_resolution_clock::now();
    mesh::SquareMaker mesh_maker(N);
    mesh::Mesh<2> mesh = mesh_maker.make_mesh();
    auto end_mesh = std::chrono::high_resolution_clock::now();

    auto start_dof = std::chrono::high_resolution_clock::now();
    fe::FESimplexP1<2> fe;
    dof::DoFHandler<2> dof_handler(mesh, std::make_unique<fe::FESimplexP1<2>>());
    dof_handler.distribute_dofs();
    auto end_dof = std::chrono::high_resolution_clock::now();

    unsigned int n_dofs = dof_handler.get_n_dofs();
    unsigned int n_dofs_per_cell = fe.get_n_dofs_per_element();

    linalg::Vector rhs(n_dofs);
    linalg::SkylinePattern csr_pattern_sym = linalg::SkylinePattern::create_from_dof_handler(dof_handler);
    linalg::SkylineMatrix A_sky(n_dofs, csr_pattern_sym);

    auto start_assembly = std::chrono::high_resolution_clock::now();
    linalg::FullMatrix local_matrix(n_dofs_per_cell);
    linalg::Vector local_rhs(n_dofs_per_cell);

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it)
    {
        auto &elem = *it;
        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);
        double volume = triangle.volume();
        mesh::Point<2> v0 = triangle.get_vertex(0);
        mesh::Point<2> v1 = triangle.get_vertex(1);
        mesh::Point<2> v2 = triangle.get_vertex(2);

        local_matrix.set_to_zero();
        local_rhs.fill(0.0);

        local_matrix(0, 0) += distance(v1, v2) * distance(v1, v2) / (4 * volume);
        local_matrix(1, 1) += distance(v2, v0) * distance(v2, v0) / (4 * volume);
        local_matrix(2, 2) += distance(v0, v1) * distance(v0, v1) / (4 * volume);

        local_matrix(0, 1) += dot(v1, v2, v2, v0) / (4 * volume);
        local_matrix(0, 2) += dot(v2, v1, v1, v0) / (4 * volume);
        local_matrix(1, 2) += dot(v2, v0, v0, v1) / (4 * volume);

        mesh::Point<2> centroid = triangle.get_centroid();
        double f_value = f(centroid.coords[0], centroid.coords[1]);
        for (unsigned int i = 0; i < n_dofs_per_cell; i++)
            local_rhs[i] += f_value * volume / 3;

        auto local_dofs = dof_handler.get_ordered_dofs_on_element(elem);
        linalg::MatrixTools::add_local_matrix_to_global(A_sky, local_matrix, local_dofs);
        linalg::MatrixTools::add_local_vector_to_global(rhs, local_rhs, local_dofs);
    }
    auto end_assembly = std::chrono::high_resolution_clock::now();

    auto start_bc = std::chrono::high_resolution_clock::now();
    linalg::MatrixTools::apply_homogeneous_dirichlet(A_sky, rhs, dof_handler, 0);
    auto end_bc = std::chrono::high_resolution_clock::now();

    auto start_solver = std::chrono::high_resolution_clock::now();
    linalg::CGSolver solver;
    linalg::Vector sol = solver.solve(A_sky, rhs);
    auto end_solver = std::chrono::high_resolution_clock::now();

    auto start_save = std::chrono::high_resolution_clock::now();
    fastfem::mesh::DataIO<2, 2> data_io(mesh, dof_handler, sol);
    data_io.save_vtx("solution_skyline.vtk");
    auto end_save = std::chrono::high_resolution_clock::now();

    auto end_total = std::chrono::high_resolution_clock::now();

    std::cout << "Timings (seconds):\n";
    std::cout << "Mesh Generation: " << std::chrono::duration<double>(end_mesh - start_mesh).count() << "s\n";
    std::cout << "DoF Distribution: " << std::chrono::duration<double>(end_dof - start_dof).count() << "s\n";
    std::cout << "Matrix and RHS Assembly: " << std::chrono::duration<double>(end_assembly - start_assembly).count() << "s\n";
    std::cout << "Apply Boundary Conditions: " << std::chrono::duration<double>(end_bc - start_bc).count() << "s\n";
    std::cout << "CG Solver Execution: " << std::chrono::duration<double>(end_solver - start_solver).count() << "s\n";
    std::cout << "Solution Saving: " << std::chrono::duration<double>(end_save - start_save).count() << "s\n";
    std::cout << "Total Execution Time: " << std::chrono::duration<double>(end_total - start_total).count() << "s\n";

    return EXIT_SUCCESS;
}
