#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"
#include "FastFem/mesh/MeshIO.hpp"
#include "FastFem/mesh/MeshAdjacency.hpp"
#include "FastFem/dof/DofHandler.hpp"
#include "FastFem/fe/FESimplexP.hpp"
#include "FastFem/linalg/Vector.hpp"

using namespace fastfem::mesh;
using namespace fastfem::fe;
using namespace fastfem::dof;
using namespace fastfem::types;

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <N>" << std::endl;
        return 1;
    }
    int N = std::atoi(argv[1]);


    SquareMaker square(N);
    Mesh<2, 2> mesh = square.make_mesh();

    std::cout << "SQUARE: Final Number of vertices: " << mesh.vtx_count() << std::endl;
    std::cout << "SQUARE: Final Number of triangles: " << mesh.elem_count() << std::endl;
    std::cout << "SQUARE: Final Number of boundary elements: " << mesh.boundary_elem_count(0) << std::endl;

    assert(mesh.elem_count() == (size_t)(2 * N * N));
    assert(mesh.vtx_count() == (size_t)((N + 1) * (N + 1)));
    assert(mesh.boundary_elem_count(0) == (size_t)(4 * N));


    MeshIO io(mesh);
    io.save_vtu("square.vtu");
    io.save_msh("square.msh");


    FESimplexP1<2, 2> fe(1);
    DoFHandler<2, 2> dof_handler(mesh, std::make_unique<FESimplexP1<2, 2>>(fe));

    unsigned int n_dofs = dof_handler.distribute_dofs();

    std::cout << "Number of DoFs using P1: " << n_dofs << std::endl;

    assert(n_dofs == (size_t)(N + 1) * (N + 1));

    // fill a random solution vector
    fastfem::linalg::Vector solution(n_dofs);

    for (unsigned int i = 0; i < n_dofs; i++)
        solution[i] = i;
    
    fastfem::mesh::DataIO<2,2> data_io(mesh, dof_handler, solution);
    data_io.save_vtx("square.vtk");

    return 0;
}