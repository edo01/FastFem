#include "mesh/cube.hpp"

namespace mesh
{

/******************************************************************************
 * Building a cube surface mesh. N is the number of subdivisions per side.
 * ***************************************************************************/
int build_cube_vertices(struct Vertex *vert, int N)
{
	int V = N + 1;
	int nvf = V * V;
	int v = 0;
	
	for (int row = 0; row < V; row++) {
		for (int col = 0; col < V; col++) {
			vert[0 * nvf + v].x = col;
			vert[0 * nvf + v].y = 0;
			vert[0 * nvf + v].z = row;

			vert[1 * nvf + v].x = N;
			vert[1 * nvf + v].y = col;
			vert[1 * nvf + v].z = row;

			vert[2 * nvf + v].x = N - col;
			vert[2 * nvf + v].y = N;
			vert[2 * nvf + v].z = row;

			vert[3 * nvf + v].x = 0;
			vert[3 * nvf + v].y = N - col;
			vert[3 * nvf + v].z = row;

			vert[4 * nvf + v].x = col;
			vert[4 * nvf + v].y = N - row;
			vert[4 * nvf + v].z = 0;

			vert[5 * nvf + v].x = col;
			vert[5 * nvf + v].y = row;
			vert[5 * nvf + v].z = N;

			v++;
		}
	}

	return 6 * nvf;
}

int build_cube_triangles(struct Triangle *tri, int N)
{
	int V = N + 1;
	int t = 0;
	for (int face = 0; face < 6; face++) {
		for (int row = 0; row < N; row++) {
			for (int col = 0; col < N; col++) {
				int v = face * V * V + row * V + col;
				tri[t++] = {v, v + 1, v + 1 + V};
				tri[t++] = {v, v + 1 + V, v + V};
			}
		}
	}
	assert(t == 12 * N * N);
	return t;
}    

void build_cube_mesh(struct Mesh *m, int N)
{
	/* Number of vertices per side = number of divisions + 1 */
	int V = N + 1;

	/* We allocate for 6 * V^2 vertices */
	int max_vert = 6 * V * V;
	m->vertices = (struct Vertex *)malloc(max_vert * sizeof(struct Vertex));
	m->vtx_count = 0;

	/* We allocate for 12 * N^2 triangles */
	int tri_count = 12 * N * N;
	m->triangles =
	    (struct Triangle *)malloc(tri_count * sizeof(struct Triangle));
	m->tri_count = 0;

	/* We fill the vertices and then the faces */
	m->vtx_count = build_cube_vertices(m->vertices, N);
	m->tri_count = build_cube_triangles(m->triangles, N);

	/* We fix-up vertex duplication */
	m->vtx_count = dedup_mesh_vertices(m);
	assert(m->vtx_count == 6 * V * V - 12 * V + 8);

	/* Rescale to unit cube centered at the origin */
	for (int i = 0; i < m->vtx_count; ++i) {
		struct Vertex *v = &m->vertices[i];
		v->x = 2 * v->x / N - 1;
		v->y = 2 * v->y / N - 1;
		v->z = 2 * v->z / N - 1;
	}
}

/**
 * The so-called spherical cube, built by simply normalizing all vertices of
 * the cube mesh so that they end up in S^2
 */
void send_cube_to_sphere(struct Vertex *vert, int vtx_count)
{
	for (int i = 0; i < vtx_count; i++) {
		struct Vertex *v = &vert[i];
		double norm = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
		v->x /= norm;
		v->y /= norm;
		v->z /= norm;
	}
}

} // mesh