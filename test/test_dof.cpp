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

    SquareMaker square3(2);
    Mesh<2, 2> mesh3 = square3.make_mesh();

    FESimplexP1<2, 2> fe4(2);
    DoFHandler<2, 2> dof_handler4(mesh3, std::make_unique<FESimplexP1<2, 2>>(fe4));

    unsigned int n_dofs4 = dof_handler4.distribute_dofs();

    std::cout << " TEST mappping " << std::endl;
    std::cout << "Number of elements: " << mesh3.elem_count() << std::endl;
    std::cout << "Number of vertices: " << mesh3.vtx_count() << std::endl;
    std::cout << "Number of DoFs using P2: " << n_dofs4 << std::endl;


    // print all elements
    for(int i=0; i<mesh3.elem_count(); ++i){
        MeshSimplex<2, 2> T = mesh3.get_mesh_element(i);
        std::cout << "Element " << i << ": ";
        for(int j=0; j<T.vertex_count(); ++j){
            std::cout << T.get_vertex(j) << " ";
        }
        std::cout << std::endl;
    }

    // print the DoFs on all elements
    for(int i=0; i<mesh3.elem_count(); ++i){
        MeshSimplex<2, 2> T = mesh3.get_mesh_element(i);
        std::vector<fastfem::types::global_dof_index_t> ordered_dofs = dof_handler4.get_ordered_dofs_on_element(T);
        std::cout << "DoFs on element " << i << ": ";
        for(int j=0; j<ordered_dofs.size(); ++j){
            std::cout << ordered_dofs[j] << " ";
        }
        std::cout << std::endl;
    }   


    return 0;
}