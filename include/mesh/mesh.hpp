#ifndef MESH_HPP
#define MESH_HPP
/**
 * @TODO: add introduction on the type of mesh that has been chosen
 */
#include "linalg/vector.hpp"

namespace mesh
{

/* A Vertex is simply a point in R^3 */
struct Vertex 
{
	double x;
	double y;
	double z;
	// stored globally
	double operator==(const Vertex &v) const 
	{
		return (x == v.x && y == v.y && z == v.z);
	}
};

/* A triangle refers to three vertices by their **index** */
struct Triangle
{
	int a;
	int b;
	int c;
};

/* A mesh is an array of vertices, and an array of triangles
 * built over these vertices
 */
struct Mesh
{
	int vtx_count;
	int tri_count;
	struct Vertex *vertices; // GEOMETRY
	struct Triangle *triangles; // TOPOLOGY
};

/* Returns a vector from its end points */
linalg::Vector vector(struct Vertex A, struct Vertex B);


}// mesh

#endif // MESH_HPP