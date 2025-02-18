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

namespace fastfem {
namespace mesh {

template <unsigned int dim, unsigned int spacedim>
int dedup_mesh_vertices(Mesh<dim, spacedim> &mesh)
{
	size_t vtx_count = 0;
	size_t V = mesh.vtx_count();
    std::vector<int> remap(V);
    
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
	HashTable<uint64_t, size_t, VertexHasher<dim, spacedim>> vtx_remap(V, hasher);

	for (size_t i = 0; i < V; ++i) {
		size_t *p = vtx_remap.get_or_set(i, vtx_count);
		if (p) {
			remap[i] = *p;
		} else {
			remap[i] = vtx_count;
			vtx_count++;
		}
	}

	/* Remap vertices */
	for (size_t i = 0; i < V; i++) {
		//m->vertices[remap[i]] = m->vertices[i];
        mesh.set_vertex(remap[i], mesh.get_vertex(i));
        
	}
	/* Remap triangle indices */
	for (auto it = mesh.elem_begin(); it != mesh.elem_end(); ++it) {
		MeshSimplex<dim,spacedim> T = *it;
        for(unsigned int j = 0; j < dim; j++){
            T.set_vertex(j, remap[T.get_vertex(j)]);
            assert(T.get_vertex(j) < vtx_count);
        }
	}
    
	return vtx_count;
}


} // mesh
} // fastfem

#endif // MESH_TOOLS_HPP