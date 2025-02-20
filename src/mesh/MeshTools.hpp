#ifndef MESHTOOLS_HPP
#define MESHTOOLS_HPP

/**
 * Utility functions for mesh manipulation.
 * No need to explicitly instantiate the functions.
 */

#include <stdint.h>
#include "FastFem/common/hash_table.h"
#include "FastFem/common/hash.h"
#include "FastFem/common/vec3.h"
#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/VertexHasher.hpp"
#include "FastFem/types/CommonTypes.hpp"

namespace fastfem {
namespace mesh {

using global_vertex_index  = fastfem::types::global_vertex_index;
using local_vertex_index   = fastfem::types::local_vertex_index;

template <unsigned int dim, unsigned int spacedim>
int dedup_mesh_vertices(Mesh<dim, spacedim> &mesh)
{
	global_vertex_index vtx_count = 0;
	global_vertex_index V = mesh.vtx_count();
    std::vector<global_vertex_index> remap(V);
    
	/*
	 * We have replaced the linear search (O(n^2)) with an
	 * hash table to deduplicate the vertices. 
	 */
	
	/*
	 * 	To build the hasher we tried to follow the design of the professor
	 *  using a custom Hasher struct. 
	 *
	 *	Our goal is to remap the vertices. To achieve that we build an 
	 *  hash table that given the index of a vertex returns the deduplicated index.
	 *  We use an Hasher struct which is able to calc the hash of the vertex given 
     *  the key (index inside the vertices array) and also being able to recognize
     *  when two different vertices have the same hash.
	 *
	 *  So the hash table will be built using the index of the vertex as key and the
	 *  deduplicated index as value. 
	 *
	 */
	VertexHasher<dim, spacedim> hasher(mesh);
	HashTable<uint64_t, global_vertex_index, VertexHasher<dim, spacedim>> vtx_remap(V, hasher);

	for (global_vertex_index i = 0; i < V; ++i) {
		size_t *p = vtx_remap.get_or_set(i, vtx_count);
		if (p) {
			remap[i] = *p;
		} else {
			remap[i] = vtx_count;
			vtx_count++;
		}
	}

	/* Remap vertices */
	for (global_vertex_index i = 0; i < V; i++) {
		//m->vertices[remap[i]] = m->vertices[i];
        mesh.set_vertex(remap[i], mesh.get_vertex(i));
        
	}
	/* Remap triangle indices */
	for (auto it = mesh.elem_begin(); it != mesh.elem_end(); ++it) {
		MeshSimplex<dim,spacedim> T = *it;
        for(local_vertex_index j = 0; j < dim; j++){
            T.set_vertex(j, remap[T.get_vertex(j)]);
            assert(T.get_vertex(j) < vtx_count);
        }
	}
    
	return vtx_count;
}


} // mesh
} // fastfem

#endif // MESH_TOOLS_HPP