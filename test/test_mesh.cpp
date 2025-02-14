#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"

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

    std::cout << "Final Number of vertices: " << mesh.vtx_count() << std::endl;
    std::cout << "Final Number of triangles: " << mesh.elem_count() << std::endl;
    
    assert(mesh.elem_count() == 12 * N * N);
    assert(mesh.vtx_count() == 6 * (N + 1) * (N + 1) - 12 * (N+1) + 8);

    SphereSurfaceMaker sphere(N);
    Mesh<2, 3> mesh2 = sphere.make_mesh();

    std::cout << "Final Number of vertices: " << mesh2.vtx_count() << std::endl;
    std::cout << "Final Number of triangles: " << mesh2.elem_count() << std::endl;

    assert(mesh2.elem_count() == 12 * N * N);
    assert(mesh2.vtx_count() == 6 * (N + 1) * (N + 1) - 12 * (N+1) + 8);
    return 0;
}