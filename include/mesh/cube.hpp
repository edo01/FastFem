#ifndef CUBE_HPP
#define CUBE_HPP

#include <assert.h>
#include "mesh.hpp"
#include "mesh_tools.hpp"

namespace mesh
{

void build_cube_mesh(struct Mesh *m, int N);

/**
 * The so-called spherical cube, built by simply normalizing all vertices of
 * the cube mesh so that they end up in S^2
 */
void send_cube_to_sphere(struct Vertex *vert, int vtx_count);

} // mesh


#endif //CUBE_HPP