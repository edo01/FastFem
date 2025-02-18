#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <memory>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"
#include "FastFem/fe/FESimplexP.hpp"
#include "FastFem/dof/DofHandler.hpp"
#include "FastFem/types/CommonTypes.hpp"

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

    assert(mesh.elem_count() == (size_t)(2 * N * N));
    assert(mesh.vtx_count() == (size_t)((N + 1) * (N + 1)));

    
    FESimplexP1<2, 2> fe(1);
    DoFHandler<2, 2> dof_handler(mesh, std::make_unique<FESimplexP1<2, 2>>(fe));

    unsigned int n_dofs = dof_handler.distribute_dofs();

    std::cout << "Number of DoFs using P1: " << n_dofs << std::endl;

    assert(n_dofs == (size_t)(N + 1) * (N + 1));

    int n_components = 2;
    FESimplexP2<2, 2> fe2(n_components);
    DoFHandler<2, 2> dof_handler2(mesh, std::make_unique<FESimplexP2<2, 2>>(fe2));

    unsigned int n_dofs2 = dof_handler2.distribute_dofs();

    std::cout << "Number of DoFs using P2: " << n_dofs2 << std::endl;

    assert(n_dofs2 == (size_t)(n_components*((N + 1) * (N + 1) + N*N + 2*N*(N+1))));
    
    n_components = 1;
    FESimplexP3<2, 2> fe3(n_components);
    DoFHandler<2, 2> dof_handler3(mesh, std::make_unique<FESimplexP3<2, 2>>(fe3));

    unsigned int n_dofs3 = dof_handler3.distribute_dofs();

    std::cout << "Number of DoFs using P3: " << n_dofs3 << std::endl;

    assert(n_dofs3 == (size_t)(n_components*((N + 1) * (N + 1) + 2*(N*N + 2*N*(N+1)) + (2 * N * N))));

    //for each element in the mesh, print the global indices of the dofs
/*     for(auto it = mesh.elem_begin(); it != mesh.elem_end(); ++it){
        MeshSimplex<2, 2> T = *it;
        std::vector<dof_index_t> dofs = dof_handler2.get_dofs_on_cell(T);
        std::cout << "Element: " << T.get_vertex(0) << " " << T.get_vertex(1) << " " << T.get_vertex(2) << std::endl;
        std::cout << "Dofs: ";
        for(auto &d : dofs){
            std::cout << d << " ";
        }
        std::cout << std::endl;
    } */
    //dof_handler.print_dofs();

    return 0;
}