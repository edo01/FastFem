#include <cassert>
#include <iostream>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"
#include "MeshTools.hpp"

namespace fastfem {
namespace mesh {

void CubeSurfaceMaker::build_cube_vertices(Mesh<2,3> &mesh) const
{
	int V = N + 1;
	
    // we reserve space for 6 * V^2 vertices
    mesh.reserve_vertices(6 * V * V);

	for (double row = 0; row < V; row++) {
		for (double col = 0; col < V; col++) {
            double Nd = this->N;

            mesh.add_vertex({col, 0, row});
            mesh.add_vertex({Nd, col, row});
            mesh.add_vertex({Nd - col, Nd, row});
            mesh.add_vertex({0, Nd - col, row});
            mesh.add_vertex({col, Nd - row, 0});
            mesh.add_vertex({col, row, Nd});
		}
	}

    assert(mesh.vtx_count() == 6 * V * V);
};

void CubeSurfaceMaker::build_cube_triangles(Mesh<2,3> &mesh) const
{
	int V = N + 1;
    mesh.reserve_elements(12 * N * N);
    
    int indices[3];
	for (int face = 0; face < 6; face++) {
		for (int row = 0; row < N; row++) {
			for (int col = 0; col < N; col++) {
				int v = face * V * V + row * V + col;
                
                indices[0] = v; indices[1] = v + 1; indices[2] = v + 1 + V;
                Simplex<2> t1(indices);
                mesh.add_element(t1);
                indices[0] = v; indices[1] = v + 1 + V; indices[2] = v + V;
                Simplex<2> t2(indices);
                mesh.add_element(t2);
			}
		}
	}
};

Mesh<2,3> CubeSurfaceMaker::make_mesh() const
{
    Mesh<2,3> mesh;

	/* We fill the vertices and then the faces */
	build_cube_vertices(mesh);
	build_cube_triangles(mesh);

	/* We fix-up vertex duplication */
	int vtx_count = dedup_mesh_vertices<2, 3>(mesh);
	
	int V = N+1;
	assert(vtx_count == (6 * V * V - 12 * V + 8));

	// resize the vertices array to the new size
	mesh.resize_vertices(vtx_count);

	// Rescale to unit cube centered at the origin
	 for (int i = 0; i < mesh.vtx_count(); ++i) {
		struct Vertex<3> &v = mesh.get_vertex(i);
		v.coords[0] = 2 * v.coords[0] / N - 1;
		v.coords[1] = 2 * v.coords[1] / N - 1;
		v.coords[2] = 2 * v.coords[2] / N - 1;
	}

	return mesh;
}


} // namespace mesh
} // namespace fastfem