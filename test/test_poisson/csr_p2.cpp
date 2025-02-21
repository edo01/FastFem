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

#define x 0.1666666666
#define y 0.1666666666
#define x2 0.0833333333
#define y2 0.0833333333
#define xy 0.0416666666
#define area 0.5


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

    //auto f_constant = [](double x, double y) { return 1; };
    auto f = [](double var1, double var2) { return 4 - 2 * (var1 * var1 + var2 * var2); };
    //auto f = [](double x, double y) { return 10*10*10*10 * std::exp(- ((x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5))/0.001); };

    linalg::Vector rhs(n_dofs);

    linalg::CSRPattern csr_pattern = linalg::CSRPattern::create_from_dof_handler(dof_handler);

    linalg::CSRMatrix A(n_dofs, csr_pattern);

    linalg::FullMatrix local_matrix(n_dofs_per_cell);
    linalg::Vector local_rhs(n_dofs_per_cell);

    local_matrix.set_to_zero();

    local_matrix(0, 0) += 32*x2 + 32*y2 + 64*xy - 48*x - 48*y + 18*area;
    local_matrix(0, 1) += -48*x2 - 16*y2  -64*xy -4*x - 28*y - 12*area;
    local_matrix(1, 0) += -48*x2 - 16*y2  -64*xy -4*x - 28*y - 12*area;
    local_matrix(0, 2) += 16*x2 + 16*xy - 16*x -4*y +3*area;
    local_matrix(2, 0) += 16*x2 + 16*xy - 16*x -4*y +3*area;
    local_matrix(0, 3) += 16*x2 + 16*y2 + 32*xy - 12*x - 12*y;
    local_matrix(3, 0) += 16*x2 + 16*y2 + 32*xy - 12*x - 12*y;
    local_matrix(0, 4) += 16*y2 + 16*xy - 16*y -4*x + 3*area;
    local_matrix(4, 0) += 16*y2 + 16*xy - 16*y -4*x + 3*area;
    local_matrix(0, 5) += -16*x2 - 48*y2 - 64*xy + 52*y + 28*x - 12*area;
    local_matrix(5, 0) += -16*x2 - 48*y2 - 64*xy + 52*y + 28*x - 12*area;
    
    local_matrix(1, 1) += 16*x2 + 16*y2 + 48*xy - 16*y + 16*area;
    local_matrix(1, 2) += -32*x2 - 16*xy + 24*x + 4*y - 4*area;
    local_matrix(2, 1) += -32*x2 - 16*xy + 24*x + 4*y - 4*area;
    local_matrix(1, 3) += -32*xy + 16*y - 16*y2 - 16*x2;
    local_matrix(3, 1) += -32*xy + 16*y - 16*y2 - 16*x2;
    local_matrix(1, 4) += 4*x - 16*xy;
    local_matrix(4, 1) += 4*x - 16*xy;
    local_matrix(1, 5) += 64*xy + 16*y2 - 16*x + 16*x2;
    local_matrix(5, 1) += 64*xy + 16*y2 - 16*x + 16*x2;

    local_matrix(2, 2) += 16*x2 + 1 - 8*x;
    local_matrix(2, 3) += 16*xy+ 4*y;
    local_matrix(3, 2) += 16*xy+ 4*y;
    //local_matrix(2, 4) += 0;
    //local_matrix(4, 2) += 0;
    local_matrix(2, 5) += -16*xy + 4*y;
    local_matrix(5, 2) += -16*xy + 4*y;

    local_matrix(3, 3) += 16*y2 + 16*x2;
    local_matrix(3, 4) += 16*xy - 4*x;
    local_matrix(4, 3) += 16*xy - 4*x;
    local_matrix(3, 5) += -16*x2 - 16*y2 -32*xy + 16*x;
    local_matrix(5, 3) += -16*x2 - 16*y2 -32*xy + 16*x;

    local_matrix(4, 4) += 16*y2 - 8*y + area;
    local_matrix(4, 5) += 16*xy + 4*x - 8*y - 4*area;
    local_matrix(5, 4) += 16*xy + 4*x - 8*y - 4*area;

    local_matrix(5, 5) += 80*y2 + 64*xy + 16*x - 64*y + 16*area;

    // apply the jacobian
    for(size_t i = 0; i < n_dofs_per_cell; ++i)
    {
        for(size_t j = 0; j < n_dofs_per_cell; ++j)
        {
            local_matrix(i, j) *= 1;   
            printf("(%d,%d)%f \t", i, j,  local_matrix(i, j));
        }
        printf("\n");
    }

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it)
    {
        auto &elem = *it;

        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);

        mesh::Point<2> v0 = triangle.get_vertex(0);
        mesh::Point<2> v1 = triangle.get_vertex(1);
        mesh::Point<2> v2 = triangle.get_vertex(2);

        mesh::Point<2> centroid = triangle.get_centroid();

        double jacobian = triangle.volume()/2;

        double lenght_01 = distance(v0, v1);
        double lenght_12 = distance(v1, v2);
        double lenght_20 = distance(v2, v0);

        local_matrix.set_to_zero();
        local_rhs.fill(0.0);

        local_matrix(0, 0) += 32*x2 + 32*y2 + 64*xy - 48*x - 48*y + 18*area;
        local_matrix(0, 1) += -48*x2 - 16*y2  -64*xy -4*x - 28*y - 12*area;
        local_matrix(1, 0) += -48*x2 - 16*y2  -64*xy -4*x - 28*y - 12*area;
        local_matrix(0, 2) += 16*x2 + 16*xy - 16*x -4*y +3*area;
        local_matrix(2, 0) += 16*x2 + 16*xy - 16*x -4*y +3*area;
        local_matrix(0, 3) += 16*x2 + 16*y2 + 32*xy - 12*x - 12*y;
        local_matrix(3, 0) += 16*x2 + 16*y2 + 32*xy - 12*x - 12*y;
        local_matrix(0, 4) += 16*y2 + 16*xy - 16*y -4*x + 3*area;
        local_matrix(4, 0) += 16*y2 + 16*xy - 16*y -4*x + 3*area;
        local_matrix(0, 5) += -16*x2 - 48*y2 - 64*xy + 52*y + 28*x - 12*area;
        local_matrix(5, 0) += -16*x2 - 48*y2 - 64*xy + 52*y + 28*x - 12*area;
        
        local_matrix(1, 1) += 16*x2 + 16*y2 + 48*xy - 16*y + 16*area;
        local_matrix(1, 2) += -32*x2 - 16*xy + 24*x + 4*y - 4*area;
        local_matrix(2, 1) += -32*x2 - 16*xy + 24*x + 4*y - 4*area;
        local_matrix(1, 3) += -32*xy + 16*y - 16*y2 - 16*x2;
        local_matrix(3, 1) += -32*xy + 16*y - 16*y2 - 16*x2;
        local_matrix(1, 4) += 4*x - 16*xy;
        local_matrix(4, 1) += 4*x - 16*xy;
        local_matrix(1, 5) += 64*xy + 16*y2 - 16*x + 16*x2;
        local_matrix(5, 1) += 64*xy + 16*y2 - 16*x + 16*x2;

        local_matrix(2, 2) += 16*x2 + 1 - 8*x;
        local_matrix(2, 3) += 16*xy+ 4*y;
        local_matrix(3, 2) += 16*xy+ 4*y;
        //local_matrix(2, 4) += 0;
        //local_matrix(4, 2) += 0;
        local_matrix(2, 5) += -16*xy + 4*y;
        local_matrix(5, 2) += -16*xy + 4*y;

        local_matrix(3, 3) += 16*y2 + 16*x2;
        local_matrix(3, 4) += 16*xy - 4*x;
        local_matrix(4, 3) += 16*xy - 4*x;
        local_matrix(3, 5) += -16*x2 - 16*y2 -32*xy + 16*x;
        local_matrix(5, 3) += -16*x2 - 16*y2 -32*xy + 16*x;

        local_matrix(4, 4) += 16*y2 - 8*y + area;
        local_matrix(4, 5) += 16*xy + 4*x - 8*y - 4*area;
        local_matrix(5, 4) += 16*xy + 4*x - 8*y - 4*area;

        local_matrix(5, 5) += 80*y2 + 64*xy + 16*x - 64*y + 16*area;

        // apply the jacobian
        for(size_t i = 0; i < n_dofs_per_cell; ++i)
        {
            for(size_t j = 0; j < n_dofs_per_cell; ++j)
            {
                local_matrix(i, j) *= jacobian;
            }
        }

        local_rhs[0] += f(centroid.coords[0], centroid.coords[1]) * (2*x2 + 2*y2 + 4*xy - 3*x - 3*y + area)*jacobian;
        local_rhs[1] += f(centroid.coords[0], centroid.coords[1]) * (-4*x2 - 4*xy + 4*x)*jacobian;
        local_rhs[2] += f(centroid.coords[0], centroid.coords[1]) * (2*x2 - x)*jacobian;
        local_rhs[3] += f(centroid.coords[0], centroid.coords[1]) * (4*xy)*jacobian;
        local_rhs[4] += f(centroid.coords[0], centroid.coords[1]) * (2*y2 - y)*jacobian;
        local_rhs[5] += f(centroid.coords[0], centroid.coords[1]) * (-4*y2 - 4*xy + 4*y)*jacobian;

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

    auto exact_f = [](double var1, double var2) { return (1 - var1 * var1) * (1 - var2 * var2); };

    linalg::Vector exact_sol(n_dofs);

    linalg::MatrixTools::interpolate(exact_sol, dof_handler, exact_f);

    std::cout << "Norm of difference: " << (sol - exact_sol).norm() << std::endl;

    mesh::DataIO<2, 2> data_io_exact(mesh, dof_handler, exact_sol);
    data_io_exact.save_vtx("exact_solution.vtk");

    return EXIT_SUCCESS;
}