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

    // auto f = [](double x, double y) { return 20; };
    auto f = [](double var1, double var2) { return 4 - 2 * (var1 * var1 + var2 * var2); };
    //auto f = [](double x, double y) { return 10*10*10*10 * std::exp(- ((x - 0.5) * (x - 0.5) - (y - 0.5) * (y - 0.5))/0.001); };

    linalg::Vector rhs(n_dofs);

    linalg::CSRPattern csr_pattern = linalg::CSRPattern::create_from_dof_handler(dof_handler);

    linalg::CSRMatrix A(n_dofs, csr_pattern);

    linalg::FullMatrix local_matrix(n_dofs_per_cell);
    linalg::Vector local_rhs(n_dofs_per_cell);

    static double stiffness_ref[6][6] = {
        {1.0, -2.0/3.0, 1.0/6.0, 0.0, 1.0/6.0, -2.0/3.0},
        {-2.0/3.0, 8.0/3.0, -2.0/3.0, -4.0/3.0, 0.0, 0.0},
        {1.0/6.0, -2.0/3.0, 1.0/2.0, 0.0, 0.0, 0.0},
        {0.0, -4.0/3.0, 0.0, 8.0/3.0, 0.0, -4.0/3.0},
        {1.0/6.0, 0.0, 0.0, 0.0, 1.0/2.0, -2.0/3.0},
        {-2.0/3.0, 0.0, 0.0, -4.0/3.0, -2.0/3.0, 8.0/3.0}
    };
/*   static double stiffness_ref[6][6] = {
        {1.33333, -0.333333,  0.166667, -0.333333, -0.166667, -0.333333},
        {-0.333333,  2.66667, -0.666667, -1.33333,  0.0,       0.0},
        { 0.166667, -0.666667,  0.5,      0.0,      0.0,       0.0},
        {-0.333333, -1.33333,   0.0,      2.66667,  0.0,      -1.33333},
        {-0.166667,  0.0,       0.0,      0.0,      0.5,      -0.666667},
        {-0.333333,  0.0,       0.0,     -1.33333, -0.666667,  2.66667}
    }; */

/*     static double stiffness_ref[6][6];

    stiffness_ref[0][0] = 32*x2 + 20*y2 + 48*xy - 48*x - 36*y + 18*area;
    stiffness_ref[0][1] = -48*x2 - 56*xy - 16*y2 + 52*x + 28*y - 12*area;
    stiffness_ref[1][0] = -48*x2 - 56*xy - 16*y2 + 52*x + 28*y - 12*area;
    stiffness_ref[0][2] = 16*x2 + 16*xy - 16*x - 4*y + 3*area;
    stiffness_ref[2][0] = 16*x2 + 16*xy - 16*x - 4*y + 3*area;
    stiffness_ref[0][3] = 16*x2 + 16*y2 + 24*xy - 12*x - 12*y;
    stiffness_ref[3][0] = 16*x2 + 16*y2 + 24*xy - 12*x - 12*y;
    stiffness_ref[0][4] = 8*y2 + 16*xy - 4*x - 14*y + 3*area;
    stiffness_ref[4][0] = 8*y2 + 16*xy - 4*x - 14*y + 3*area;
    stiffness_ref[0][5] = -16*x2 - 32*y2 - 56*xy + 28*x + 44*y - 12*area;
    stiffness_ref[5][0] = -16*x2 - 32*y2 - 56*xy + 28*x + 44*y - 12*area;
    stiffness_ref[1][1] = 80*x2 + 16*y2 + 64*xy - 64*x - 32*y + 16*area;
    stiffness_ref[1][2] = -32*x2 - 16*xy + 24*x + 4*y - 4*area;
    stiffness_ref[2][1] = -32*x2 - 16*xy + 24*x + 4*y - 4*area;
    stiffness_ref[1][3] = -16*x2 - 32*xy - 16*y2 + 16*y;
    stiffness_ref[3][1] = -16*x2 - 32*xy - 16*y2 + 16*y;
    stiffness_ref[1][4] = 4*x - 16*xy;
    stiffness_ref[4][1] = 4*x - 16*xy;
    stiffness_ref[1][5] = 16*x2 + 16*y2 + 64*xy - 16*x - 16*y;
    stiffness_ref[5][1] = 16*x2 + 16*y2 + 64*xy - 16*x - 16*y;
    stiffness_ref[2][2] = 16*x2 - 8*x + 1*area;
    stiffness_ref[2][3] = 16*xy - 4*y;
    stiffness_ref[3][2] = 16*xy - 4*y;
    stiffness_ref[2][4] = 0;
    stiffness_ref[4][2] = 0;
    stiffness_ref[2][5] = 4*y - 16*xy;
    stiffness_ref[5][2] = 4*y - 16*xy;
    stiffness_ref[3][3] = 16*x2 + 16*y2;
    stiffness_ref[3][4] = -4*x + 16*xy;
    stiffness_ref[4][3] = -4*x + 16*xy;
    stiffness_ref[3][5] = -16*x2 - 16*y2 - 32*xy + 16*x;
    stiffness_ref[5][3] = -16*x2 - 16*y2 - 32*xy + 16*x;
    stiffness_ref[4][4] = 16*y2 - 8*y + 1*area;
    stiffness_ref[4][5] = -32*y2 - 16*xy + 4*x + 24*y - 4*area;
    stiffness_ref[5][4] = -32*y2 - 16*xy + 4*x + 24*y - 4*area;
    stiffness_ref[5][5] = 16*x2 + 80*y2 + 64*xy - 32*x - 64*y + 16*area; */

    //print stiffness matrix
    for(size_t i = 0; i < n_dofs_per_cell; ++i)
    {
        for(size_t j = 0; j < n_dofs_per_cell; ++j)
        {
            std::cout << "(" << i << "," << j << ")" << stiffness_ref[i][j] << "\t";
        }
        std::cout << std::endl;
    }



    static double mass_matrix[6][6] = {
        0.0611078, -0.0111113,  0.00277811, -0.0444454, -0.019444,  -0.0333353,
        -0.0111111,  0.0888891,  0,  0.0444463, -0.0111112,  0.0444437,
        0.00277821, 0 , 0.0166667, 0, -0.00277767, -0.0111116,
       -0.0444454,  0.0444448, 0,  0.0888916, 0,  0.044444,
       -0.0194451, -0.0111117, -0.00277739,  0,  0.0166668, 0,
       -0.0333338,  0.0444455, -0.0111112,  0.044446, 0,  0.0888872
    };

    //dubious about f term, it has too many zeroes. if the rhs is manually filled at random, the solution obtained doesnt look that bad. maybe we need to interpolate f better. rhs has always norm very close to zero
    static double shape_integral_on_ref[6] = {0, 1.0/6, 0, 1.0/6, 0, 1.0/6};
    //static double shape_integral_on_ref[6] = {1.0/6, 0, 1.0/6, 0, 1.0/6, 0};

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


        for(types::local_dof_index i = 0; i < n_dofs_per_cell; ++i)
        {
            for(types::local_dof_index j = 0; j < n_dofs_per_cell; ++j)
            {
                //local_matrix(i, j) += stiffness_ref[i][j] * 2 / (volume);
                local_matrix(i, j) = stiffness_ref[i][j] * 2 * volume;
            }
            // avarage of the function f over the element
            double avg = f(v0[0], v0[1]) + f(v1[0], v1[1]) + f(v2[0], v2[1]);
            avg /= 3.0;

            local_rhs[i] += avg * shape_integral_on_ref[i] * 2 * volume;
            
            //local_rhs[i] += f(centroid.coords[0], centroid.coords[1]) * shape_integral_on_ref[i] * 2 * volume;
        }

        /* linalg::Vector f_local(n_dofs_per_cell);

        for(types::local_dof_index dof_i = 0; dof_i < n_dofs_per_cell; ++dof_i)
        {
            mesh::Point<2> dof_coords = fe.get_dof_coords(triangle, dof_i);
            f_local[dof_i] = f(dof_coords.coords[0], dof_coords.coords[1]);
        }
        
        
        // multiply local rhs by the mass matrix
        for(types::local_dof_index i = 0; i < n_dofs_per_cell; ++i)
        {
            for(types::local_dof_index j = 0; j < n_dofs_per_cell; ++j)
            {
                local_rhs[i] += mass_matrix[i][j] * f_local[j] * volume *2;
            } 
        }*/

        auto local_dofs = dof_handler.get_ordered_dofs_on_element(elem);
        linalg::MatrixTools::add_local_matrix_to_global(A, local_matrix, local_dofs);
        linalg::MatrixTools::add_local_vector_to_global(rhs, local_rhs, local_dofs);

    } 


    linalg::MatrixTools::apply_homogeneous_dirichlet(A, rhs, dof_handler, 0);

    //A.print();

    std::cout << "norm of rhs: " << rhs.norm() << std::endl;
    //output the rhs to see if it is correct

    fastfem::mesh::DataIO<2, 2> rhs_io(mesh, dof_handler, rhs);
    rhs_io.save_vtx("rhs.vtk");

    linalg::CGSolver solver(1000, 1e-7);
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