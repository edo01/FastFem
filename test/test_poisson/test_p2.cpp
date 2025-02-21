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

    fe::FESimplexP2<2> fe;
    dof::DoFHandler<2> dof_handler(mesh);

    dof_handler.distribute_dofs(std::make_shared<fe::FESimplexP2<2>>(fe));

    unsigned int n_dofs = dof_handler.get_n_dofs();
    unsigned int n_dofs_per_cell = fe.get_n_dofs_per_element();

    // auto f = [](double x, double y) { return 20; };
    auto f = [](double x, double y) { return 4 - 2 * (x * x + y * y); };
    //auto f = [](double x, double y) { return 10*10*10*10 * std::exp(- ((x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5))/0.001); };

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

        local_matrix.set_to_zero();
        local_rhs.fill(0.0);

        static double stiffness_ref[6][6] = {
            {1.0, -2.0/3.0, 1.0/6.0, 0.0, 1.0/6.0, -2.0/3.0},
            {-2.0/3.0, 8.0/3.0, -2.0/3.0, -4.0/3.0, 0.0, 0.0},
            {1.0/6.0, -2.0/3.0, 1.0/2.0, 0.0, 0.0, 0.0},
            {0.0, -4.0/3.0, 0.0, 8.0/3.0, 0.0, -4.0/3.0},
            {1.0/6.0, 0.0, 0.0, 0.0, 1.0/2.0, -2.0/3.0},
            {-2.0/3.0, 0.0, 0.0, -4.0/3.0, -2.0/3.0, 8.0/3.0}
        };

        //dubious about f term, it has too many zeroes. if the rhs is manually filled at random, the solution obtained doesnt look that bad. maybe we need to interpolate f better. rhs has always norm very close to zero
        static double shape_integral_on_ref[6] = {0, 1/6, 0, 1/6, 0, 1/6};

        for(types::local_dof_index i = 0; i < n_dofs_per_cell; ++i)
        {
            for(types::local_dof_index j = 0; j < n_dofs_per_cell; ++j)
            {
                local_matrix(i, j) += stiffness_ref[i][j] / (2 * volume);
            }
            local_rhs[i] += f(centroid.coords[0], centroid.coords[1]) * shape_integral_on_ref[i] * volume;
        }

        auto local_dofs = dof_handler.get_ordered_dofs_on_element(elem);
        linalg::MatrixTools::add_local_matrix_to_global(A, local_matrix, local_dofs);
        linalg::MatrixTools::add_local_vector_to_global(rhs, local_rhs, local_dofs);

    } 

    //instead doing this manually, we get reasonable results
    for(auto it = dof_handler.elem_begin() + 1000; it < dof_handler.elem_begin() + 1300; it++)
    {
        auto global_dofs = dof_handler.get_ordered_dofs_on_element(*it);
        for(auto dof : global_dofs)
        {
            rhs[dof] = 5;
        }
    }

    linalg::MatrixTools::apply_homogeneous_dirichlet(A, rhs, dof_handler, 0);

    //A.print();

    std::cout << "norm of rhs: " << rhs.norm() << std::endl;

    linalg::CGSolver solver(1000, 1e-7);
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