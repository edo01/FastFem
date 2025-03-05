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

#include <cstdlib>

static inline double distance(const mesh::Point<2> &v1, const mesh::Point<2> &v2)
{
    return std::sqrt((v2.coords[0] - v1.coords[0]) * (v2.coords[0] - v1.coords[0]) + (v2.coords[1] - v1.coords[1]) * (v2.coords[1] - v1.coords[1]));
}

static inline double dot(const mesh::Point<2> &v1, const mesh::Point<2> &v2, const mesh::Point<2> &w1, const mesh::Point<2> &w2)
{
    return (v2.coords[0] - v1.coords[0]) * (w2.coords[0] - w1.coords[0]) + (v2.coords[1] - v1.coords[1]) * (w2.coords[1] - w1.coords[1]);
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <N>" << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int N = std::atoi(argv[1]);

    auto f = [](double var1, double var2) { return 4 - 2 * (var1 * var1 + var2 * var2); };
    auto exact_f = [](double var1, double var2) { return (1 - var1 * var1) * (1 - var2 * var2); };
    
    /**
     * CREATE THE MESH
     */
     mesh::SquareMaker mesh_maker(N);
     mesh::Mesh<2> mesh = mesh_maker.make_mesh();
 
    std::cout << "============ MESH INFO =============" << std::endl;
    std::cout << mesh.vtx_count() << " vertices" << std::endl;
    std::cout << mesh.elem_count() << " elements" << std::endl;
    std::cout << mesh.boundary_elem_count(0) << " boundary elements" << std::endl;
    std::cout << "====================================" << std::endl << std::endl;

    /**
     * CREATE THE FE AND DOF HANDLER
     */
    fe::FESimplexP1<2> fe;
    dof::DoFHandler<2> dof_handler(mesh);

    // create the solver
    linalg::CGSolver solver(1000, 1e-7);

    /**
     * DOFs DISTRIBUTION
     */
    dof_handler.distribute_dofs(std::make_shared<fe::FESimplexP1<2>>(fe));

    unsigned int n_dofs = dof_handler.get_n_dofs();
    unsigned int n_dofs_per_cell = fe.get_n_dofs_per_element();

    std::cout << "========= DOF DISTRIBUTION =========" << std::endl;
    std::cout << n_dofs << " total dofs" << std::endl;
    std::cout << n_dofs_per_cell << " dofs per element" << std::endl;
    std::cout << "====================================" << std::endl << std::endl;

    /**
     * INITIALIZE THE LINEAR SYSTEM
     */
    linalg::Vector rhs(n_dofs);
    linalg::SkylinePattern skyline_pattern = linalg::SkylinePattern::create_from_dof_handler(dof_handler);
    linalg::SkylineMatrix A(n_dofs, skyline_pattern);

    /**
     * ASSEMBLE THE LINEAR SYSTEM
     */
    linalg::FullMatrix local_matrix(n_dofs_per_cell);
    linalg::Vector local_rhs(n_dofs_per_cell);

    // shape integral on the reference triangle
    static double shape_integral_on_ref[3] = {1.0/6, 1.0/6, 1.0/6};

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it)
    {
        auto &elem = *it;

        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);
        mesh::Point<2> v0 = triangle.get_vertex(0);
        mesh::Point<2> v1 = triangle.get_vertex(1);
        mesh::Point<2> v2 = triangle.get_vertex(2);

        double volume = triangle.volume();

        double length01 = distance(v0, v1);
        double length12 = distance(v1, v2);
        double length20 = distance(v2, v0);

        // reset the local matrix and rhs
        local_matrix.set_to_zero();
        local_rhs.fill(0.0);

        local_matrix(0, 0) += length12 * length12 / (4 * volume);
        local_matrix(1, 1) += length20 * length20 / (4 * volume);
        local_matrix(2, 2) += length01 * length01 / (4 * volume);

        local_matrix(0, 1) += dot(v1, v2, v2, v0) / (4 * volume);
        local_matrix(0, 2) += dot(v2, v1, v1, v0) / (4 * volume);
        local_matrix(1, 2) += dot(v2, v0, v0, v1) / (4 * volume);

    
        // average of the function f over the element
        double avg = f(v0[0], v0[1]) + f(v1[0], v1[1]) + f(v2[0], v2[1]);
        avg /= 3.0; 

        for(types::local_dof_index i = 0; i < n_dofs_per_cell; ++i)
        {
            /*
             * Approximation of the integral of f over the element
             */
            local_rhs[i] += avg * shape_integral_on_ref[i] * 2 * volume;
        }

        auto local_dofs = dof_handler.get_ordered_dofs_on_element(elem);
        linalg::MatrixTools::add_local_matrix_to_global(A, local_matrix, local_dofs);
        linalg::MatrixTools::add_local_vector_to_global(rhs, local_rhs, local_dofs);

    } 

    // apply homogeneous dirichlet boundary conditions
    linalg::MatrixTools::apply_homogeneous_dirichlet(A, rhs, dof_handler, 0);

    /**
     * SOLVE THE LINEAR SYSTEM
     */
    std::cout << "Solving the linear system..." << std::endl << std::endl;
    linalg::Vector sol = solver.solve(A, rhs);

    /**
     * SAVE THE SOLUTION
     */

    fastfem::mesh::DataIO<2, 2> data_io(mesh, dof_handler, sol);
    data_io.save_vtx("solution_csr.vtk");

    linalg::Vector exact_sol(n_dofs);

    // interpolate the exact solution
    linalg::MatrixTools::interpolate(exact_sol, dof_handler, exact_f);

    std::cout << "===== SOLUTION INFO =====" << std::endl;
    std::cout << "CG converged in " << solver.get_last_step() << " iterations" << std::endl;
    std::cout << "Last residual: " << solver.get_error() << std::endl;
    std::cout << "Max of solution: " << sol.max() << std::endl;
    std::cout << "Norm of difference: " << (sol - exact_sol).norm() << std::endl;
    std::cout << "=========================" << std::endl << std::endl;

    mesh::DataIO<2, 2> data_io_exact(mesh, dof_handler, exact_sol);
    data_io_exact.save_vtx("exact_solution.vtk");

    return EXIT_SUCCESS;
}