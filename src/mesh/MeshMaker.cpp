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
    
    size_t indices[3];
	for (int face = 0; face < 6; face++) {
		for (int row = 0; row < N; row++) {
			for (int col = 0; col < N; col++) {
				int v = face * V * V + row * V + col;
                
                indices[0] = v; indices[1] = v + 1; indices[2] = v + 1 + V;
                MeshSimplex<2,3> t1(indices);
                mesh.add_element(t1);
                indices[0] = v; indices[1] = v + 1 + V; indices[2] = v + V;
                MeshSimplex<2, 3> t2(indices);
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
	for (auto v = mesh.vtx_begin(); v != mesh.vtx_end(); ++v) {
		v->point.coords[0] = 2 * v->point.coords[0] / N - 1;
		v->point.coords[1] = 2 * v->point.coords[1] / N - 1;
		v->point.coords[2] = 2 * v->point.coords[2] / N - 1;
	}

	return mesh;
};

void SphereSurfaceMaker::sendPointsToSphere(Mesh<2,3> &mesh) const{
	/**
	 * The so-called spherical cube, built by simply normalizing all vertices of
	 * the cube mesh so that they end up in S^2
	 */
	//use foreach on the vertices of the mesh
	for(auto v = mesh.vtx_begin(); v != mesh.vtx_end(); ++v){
		double norm = sqrt(v->point.coords[0] * v->point.coords[0] + v->point.coords[1] * v->point.coords[1] + v->point.coords[2] * v->point.coords[2]);
		v->point.coords[0] /= norm;
		v->point.coords[1] /= norm;
		v->point.coords[2] /= norm;
	}
};

Mesh<2,3> SphereSurfaceMaker::make_mesh() const{
	Mesh<2,3> mesh = cubeSurfaceMaker.make_mesh();
	sendPointsToSphere(mesh);
	return mesh;
};

void SquareMaker::build_square_vertices(Mesh<2,2>& mesh) const
{
	int V = N + 1;
	
	// we reserve space for V^2 vertices
	mesh.reserve_vertices(V * V);

	for (double row = 0; row < V; row++) {
		for (double col = 0; col < V; col++) {
			mesh.add_vertex({col, row});
		}
	}

	assert(mesh.vtx_count() == V * V);
}

void SquareMaker::build_square_triangles(Mesh<2,2>& mesh) const
{
	int V = N + 1;
	mesh.reserve_elements(2 * N * N);
	
	size_t indices[3];
	for (int row = 0; row < N; row++) {
		for (int col = 0; col < N; col++) {
			int v = row * V + col;
			
			indices[0] = v; indices[1] = v + 1; indices[2] = v + 1 + V;
			MeshSimplex<2,2> t1(indices);
			mesh.add_element(t1);
			indices[0] = v; indices[1] = v + 1 + V; indices[2] = v + V;
			MeshSimplex<2,2> t2(indices);
			mesh.add_element(t2);
		}
	}

	assert(mesh.elem_count() == 2 * N * N);
}

void SquareMaker::build_square_boundaries(Mesh<2,2>& mesh) const{

	size_t indices[2]; // indices of the dofs of the edges on the boundary
	int V = N + 1;

	// we reserve space for 4 * N elements
	mesh.reserve_boundary_elements(4 * N);
	
	for(int i = 0; i < N; i++){
		// top edge
		indices[0] = i; indices[1] = i + 1;
		MeshSimplex<1,2> e1(indices);
		mesh.add_boundary_element(e1);
		
		// left edge
		indices[0] = i * V; indices[1] = (i + 1) * V;
		MeshSimplex<1,2> e2(indices);
		mesh.add_boundary_element(e2);

		// bottom edge
		indices[0] = V * (V - 1) + i; indices[1] = V * (V - 1) + i + 1;
		MeshSimplex<1,2> e3(indices);
		mesh.add_boundary_element(e3);

		// right edge
		indices[0] = V * i + V - 1; indices[1] = V * (i+1) + V - 1;
		MeshSimplex<1,2> e4(indices);
		mesh.add_boundary_element(e4);
	}

	assert(mesh.boundary_elem_count() == 4 * N);
}

Mesh<2,2> SquareMaker::make_mesh() const {
	Mesh<2,2> mesh;

	build_square_vertices(mesh);
	build_square_triangles(mesh);
	build_square_boundaries(mesh);

	// Rescale to unit cube centered at the origin
	for (auto v = mesh.vtx_begin(); v != mesh.vtx_end(); ++v) {
		v->point.coords[0] = 2 * v->point.coords[0] / N - 1;
		v->point.coords[1] = 2 * v->point.coords[1] / N - 1;
	}
	
	return mesh;
}

} // namespace mesh
} // namespace fastfem