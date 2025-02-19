#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"
#include "FastFem/mesh/MeshIO.hpp"
#include "FastFem/mesh/MeshAdjacency.hpp"

using namespace fastfem::mesh;
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

    MeshIO io(mesh);
    io.save_vtu("square.vtu");
    io.save_msh("square.msh");

    return 0;
}