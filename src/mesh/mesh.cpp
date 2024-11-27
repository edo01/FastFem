#include "mesh/mesh.hpp"

namespace mesh
{

/* Returns a vector from its end points */
linalg::Vector vector(struct Vertex A, struct Vertex B)
{
	// From the affine space to the space of vectors
	linalg::Vector res;
	res.x = B.x - A.x;
	res.y = B.y - A.y;
	res.z = B.z - A.z;
	return res;
}

}