#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

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

#define SQUARE(x) ((x) * (x))

using namespace fastfem;

//-----------------------------
// Utility functions
//-----------------------------
double distance(const mesh::Point<2> &v1, const mesh::Point<2> &v2)
{
    return std::sqrt((v2.coords[0] - v1.coords[0]) * (v2.coords[0] - v1.coords[0]) +
                     (v2.coords[1] - v1.coords[1]) * (v2.coords[1] - v1.coords[1]));
}

double dot(const mesh::Point<2> &v1, const mesh::Point<2> &v2,
           const mesh::Point<2> &w1, const mesh::Point<2> &w2)
{
    return (v2.coords[0] - v1.coords[0]) * (w2.coords[0] - w1.coords[0]) +
           (v2.coords[1] - v1.coords[1]) * (w2.coords[1] - w1.coords[1]);
}

double compute_den(double xa, double ya, double xb, double yb, double xc, double yc) {
    return SQUARE(xc * (-ya + yb) + xb * (ya - yc) + xa * (-yb + yc));
}

void compute_stiffness_loc(mesh::Simplex<2, 2> triangle, double matrix[6][6]) {

    double xa = triangle.get_vertex(0).coords[0];
    double ya = triangle.get_vertex(0).coords[1];
    double xb = triangle.get_vertex(1).coords[0];
    double yb = triangle.get_vertex(1).coords[1];
    double xc = triangle.get_vertex(2).coords[0];
    double yc = triangle.get_vertex(2).coords[1];

    matrix[0][0] = (SQUARE(xb - xc) + SQUARE(yb - yc)) / 2.;
    matrix[0][1] = (-2 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 3.;
    matrix[0][2] = ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc)) / 6.;
    matrix[0][3] = 0;
    matrix[0][4] = (-1 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 6.;
    matrix[0][5] = (2 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 3.;

    matrix[1][0] = (-2 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 3.;
    matrix[1][1] = (4 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - (ya + yb) * yc + SQUARE(yc))) / 3.;
    matrix[1][2] = (-2 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 3.;
    matrix[1][3] = (4 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 3.;
    matrix[1][4] = 0;
    matrix[1][5] = (-4 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 3.;

    matrix[2][0] = ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc)) / 6.;
    matrix[2][1] = (-2 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 3.;
    matrix[2][2] = (SQUARE(xa - xc) + SQUARE(ya - yc)) / 2.;
    matrix[2][3] = (-2 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 3.;
    matrix[2][4] = (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc)) / 6.;
    matrix[2][5] = 0;

    matrix[3][0] = 0;
    matrix[3][1] = (4 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 3.;
    matrix[3][2] = (-2 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 3.;
    matrix[3][3] = (4 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - (ya + yb) * yc + SQUARE(yc))) / 3.;
    matrix[3][4] = (-2 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 3.;
    matrix[3][5] = (-4 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 3.;

    matrix[4][0] = (-1 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 6.;
    matrix[4][1] = 0;
    matrix[4][2] = (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc)) / 6.;
    matrix[4][3] = (-2 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 3.;
    matrix[4][4] = (SQUARE(xa - xb) + SQUARE(ya - yb)) / 2.;
    matrix[4][5] = (2 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 3.;

    matrix[5][0] = (2 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 3.;
    matrix[5][1] = (-4 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 3.;
    matrix[5][2] = 0;
    matrix[5][3] = (-4 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 3.;
    matrix[5][4] = (2 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 3.;
    matrix[5][5] = (4 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - (ya + yb) * yc + SQUARE(yc))) / 3.;

    // multiply by jacobian
    double jacobian = triangle.volume() * 2;
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 6; ++j)
        {
            matrix[i][j] *= jacobian/compute_den(xa, ya, xb, yb, xc, yc);
        }
    }
}

//-----------------------------
// Function: run_p1_simulation
// This function runs the p1 finite element simulation (assumed to use FESimplexP1)
// and returns the error norm.
//-----------------------------
double run_p1_simulation(unsigned int N)
{
    // Create mesh
    mesh::SquareMaker mesh_maker(N);
    mesh::Mesh<2> mesh = mesh_maker.make_mesh();

    // Use P1 finite element
    fe::FESimplexP1<2> fe;
    dof::DoFHandler<2> dof_handler(mesh);
    dof_handler.distribute_dofs(std::make_shared<fe::FESimplexP1<2>>(fe));

    unsigned int n_dofs = dof_handler.get_n_dofs();
    unsigned int n_dofs_per_cell = fe.get_n_dofs_per_element();

    auto f = [](double var1, double var2) { return 4 - 2 * (var1 * var1 + var2 * var2); };

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

    auto exact_f = [](double var1, double var2) { return (1 - var1 * var1) * (1 - var2 * var2); };
    linalg::Vector exact_sol(n_dofs);
    linalg::MatrixTools::interpolate(exact_sol, dof_handler, exact_f);
    double error = (sol - exact_sol).norm()/*/exact_sol.norm()*/;

    return error;
}

//-----------------------------
// Function: run_p2_simulation
// This function runs the p2 finite element simulation (the code you provided)
// and returns the error norm.
//-----------------------------
double run_p2_simulation(unsigned int N) {
    mesh::SquareMaker mesh_maker(N);
    mesh::Mesh<2> mesh = mesh_maker.make_mesh();

    fe::FESimplexP2<2> fe;
    dof::DoFHandler<2> dof_handler(mesh);
    dof_handler.distribute_dofs(std::make_shared<fe::FESimplexP2<2>>(fe));

    unsigned int n_dofs = dof_handler.get_n_dofs();
    unsigned int n_dofs_per_cell = fe.get_n_dofs_per_element();

    auto f = [](double var1, double var2) { return 4 - 2 * (var1 * var1 + var2 * var2); };

    linalg::Vector rhs(n_dofs);
    linalg::CSRPattern csr_pattern = linalg::CSRPattern::create_from_dof_handler(dof_handler);
    linalg::CSRMatrix A(n_dofs, csr_pattern);

    linalg::FullMatrix local_matrix(n_dofs_per_cell);
    linalg::Vector local_rhs(n_dofs_per_cell);

    static double shape_integral_on_ref[6] = {0, 1.0 / 6, 0, 1.0 / 6, 0, 1.0 / 6};

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it) {
        auto &elem = *it;
        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);

        double stiffness[6][6];
        compute_stiffness_loc(triangle, stiffness);

        local_matrix.set_to_zero();
        local_rhs.fill(0.0);

        mesh::Point<2> v0 = triangle.get_vertex(0);
        mesh::Point<2> v1 = triangle.get_vertex(1);
        mesh::Point<2> v2 = triangle.get_vertex(2);

        double volume = triangle.volume();

        double average = f(v0[0], v0[1]) + f(v1[0], v1[1]) + f(v2[0], v2[1]);
        average /= 3;

        for (types::local_dof_index i = 0; i < n_dofs_per_cell; ++i) {
            for (types::local_dof_index j = 0; j < n_dofs_per_cell; ++j) {
                local_matrix(i, j) += stiffness[i][j];
            }
            local_rhs[i] += /*f(centroid.coords[0], centroid.coords[1])*/ average * shape_integral_on_ref[i] * 2 * volume;
        }

        auto local_dofs = dof_handler.get_ordered_dofs_on_element(elem);
        linalg::MatrixTools::add_local_matrix_to_global(A, local_matrix, local_dofs);
        linalg::MatrixTools::add_local_vector_to_global(rhs, local_rhs, local_dofs);
    }

    linalg::MatrixTools::apply_homogeneous_dirichlet(A, rhs, dof_handler, 0);

    linalg::CGSolver solver(2000, 1e-7);
    linalg::Vector sol = solver.solve(A, rhs);

    auto exact_f = [](double var1, double var2) { return (1 - var1 * var1) * (1 - var2 * var2); };
    linalg::Vector exact_sol(n_dofs);
    linalg::MatrixTools::interpolate(exact_sol, dof_handler, exact_f);

    return (sol - exact_sol).norm()/*/exact_sol.norm()*/;
}
    

//-----------------------------
// Main function
//-----------------------------
int main()
{
    // Open CSV file to store the convergence results
    std::ofstream results("../plot/convergence_results.csv");
    if (!results.is_open()) {
        std::cerr << "Error: Could not open CSV file for writing." << std::endl;
        return EXIT_FAILURE;
    }
    results << "N,error_p1,error_p2,theoretical_p1,theoretical_p2\n";

    // Theoretical convergence orders
    double theoretical_p1 = 2.0;
    double theoretical_p2 = 3.0;

    // Loop over mesh sizes from 20 to 320, doubling each time
    for (unsigned int N = 20; N <= 320; N *= 2)
    {
        double error_p1 = run_p1_simulation(N);
        double error_p2 = run_p2_simulation(N);

        // Write one row per mesh size in the CSV file
        results << N << "," << error_p1 << "," << error_p2 << "," << theoretical_p1 << "," << theoretical_p2 << "\n";

        // Also print results to console (optional)
        std::cout << "N = " << N 
                  << " | error_p1 = " << error_p1 
                  << " | error_p2 = " << error_p2 
                  << std::endl;
    }

    results.close();
    return EXIT_SUCCESS;
}
