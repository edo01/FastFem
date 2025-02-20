#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"
#include "FastFem/mesh/MeshIO.hpp"
#include "FastFem/mesh/MeshAdjacency.hpp"

/**
 * @brief Test the mesh generation for different geometries
 * 
 */
using namespace fastfem::mesh;
int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <N>" << std::endl;
        return 1;
    }
    int N = std::atoi(argv[1]);
    CubeSurfaceMaker cube(N);
    Mesh<2, 3> mesh = cube.make_mesh();

    std::cout << "CUBE: Final Number of vertices: " << mesh.vtx_count() << std::endl;
    std::cout << "CUBE: Final Number of triangles: " << mesh.elem_count() << std::endl;
    
    assert(mesh.elem_count() == (size_t) (12 * N * N));
    assert(mesh.vtx_count() == (size_t)(6 * (N + 1) * (N + 1) - 12 * (N+1) + 8));

    SphereSurfaceMaker sphere(N);
    Mesh<2, 3> mesh2 = sphere.make_mesh();

    std::cout << "SPHERE: Final Number of vertices: " << mesh2.vtx_count() << std::endl;
    std::cout << "SPHERE: Final Number of triangles: " << mesh2.elem_count() << std::endl;

    assert(mesh2.elem_count() == (size_t)(12 * N * N));
    assert(mesh2.vtx_count() == (size_t)(6 * (N + 1) * (N + 1) - 12 * (N+1) + 8));

    SquareMaker square(N);
    Mesh<2, 2> mesh3 = square.make_mesh();

    std::cout << "SQUARE: Final Number of vertices: " << mesh3.vtx_count() << std::endl;
    std::cout << "SQUARE: Final Number of triangles: " << mesh3.elem_count() << std::endl;
    std::cout << "SQUARE: Final Number of boundary elements: " << mesh3.boundary_elem_count(0) << std::endl;

    assert(mesh3.elem_count() == (size_t)(2 * N * N));
    assert(mesh3.vtx_count() == (size_t)((N + 1) * (N + 1)));
    assert(mesh3.boundary_elem_count(0) == (size_t)(4 * N));

    // skip the test if the mesh is too small
    if(N<100) return 0;

    MeshAdjacency<2,2> adj(mesh3);
    //adj.print_adjacency();
    MeshSimplex<2,2> T = mesh3.get_mesh_element(2*(N+1)+(N+1)/2);
    std::vector<MeshSimplex<2,2>> adj_simplices = adj.get_adjacent_simplices(T);
    std::cout << "Number of adjacent simplices to the element "<< 2*(N+1)+(N+1)/2 << ": " << adj_simplices.size() << std::endl;
    assert(adj_simplices.size() == 12);
    adj_simplices.clear();
    T = mesh3.get_mesh_element(0);
    adj_simplices = adj.get_adjacent_simplices(T);
    std::cout << "Number of adjacent simplices to the element 0: " << adj_simplices.size() << std::endl;
    assert(adj_simplices.size() == 6);



    return 0;
}