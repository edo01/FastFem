#ifndef MESH_TOOLS_HPP
#define MESH_TOOLS_HPP

#include <stdint.h>
#include "common/hash_table.h"
#include "common/hash.h"
#include "common/vec3.h"
#include "mesh/mesh.hpp"

namespace mesh
{

int dedup_mesh_vertices(struct Mesh *m);

} // mesh

#endif // MESH_TOOLS_HPP