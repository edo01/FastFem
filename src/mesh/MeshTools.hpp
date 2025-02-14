#ifndef MESHTOOLS_HPP
#define MESHTOOLS_HPP

#include <stdint.h>
#include "FastFem/common/hash_table.h"
#include "FastFem/common/hash.h"
#include "FastFem/common/vec3.h"
#include "FastFem/mesh/Mesh.hpp"
#include <functional>

namespace fastfem {
namespace mesh {

/*	
 *	The Hasher class must implement 'hash', 'is_empty', 'is_equal'
 *	functions, and reserve a key named empty_key for marking empty 
 *	slots.  
 */
template <unsigned int dim, unsigned int spacedim>
struct VertexHasher 
{
	//const Vertex<spacedim> *pos;
    Mesh<dim, spacedim> &mesh;

	static constexpr size_t empty_key = ~static_cast<size_t>(0); // reserved key to mark empty slots
	size_t hash(size_t key) const
	{
        std::size_t hash = 0;

        std::hash<double> hasher;
        std::hash<int> hasher_int;
        Vertex<spacedim> v = mesh.get_vertex(key);
        
        for(int i=0; i<spacedim; i++){
            std::size_t hvalue = hasher(v.coords[i]);
            hash = hashCombine(hash, hvalue);
            std::size_t hindex = hasher_int(i);
            hash = hashCombine(hash, hindex);
        }

        return hash;
        
	}

	bool is_empty(size_t key) const { return (key == empty_key); }

	bool is_equal(size_t key1, size_t key2) const
	{   
		return mesh.get_vertex(key1) == mesh.get_vertex(key2); // redefinition of the equality operator
	}

    std::size_t hashCombine(std::size_t h1, std::size_t h2) const {
        return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

template <unsigned int dim, unsigned int spacedim>
int dedup_mesh_vertices(Mesh<dim, spacedim> &mesh)
{
	int vtx_count = 0;
	int V = mesh.vtx_count();
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
	VertexHasher<dim, spacedim> hasher{mesh};
	HashTable<uint64_t, size_t, VertexHasher<dim, spacedim>> vtx_remap(V, hasher);

	for (size_t i = 0; i < V; ++i) {
		size_t *p = vtx_remap.get_or_set(i, vtx_count);
		if (p) {
            Vertex<spacedim> v = mesh.get_vertex(i);
            Vertex<spacedim> v2 = mesh.get_vertex(*p);
			remap[i] = *p;
		} else {
			remap[i] = vtx_count;
			vtx_count++;
		}
	}

	/* Remap vertices */
	for (int i = 0; i < V; i++) {
		//m->vertices[remap[i]] = m->vertices[i];
        mesh.set_vertex(remap[i], mesh.get_vertex(i));
        
	}
	/* Remap triangle indices */
	for (int i = 0; i < mesh.elem_count(); i++) {
		Simplex<dim> T = mesh.get_element(i);
        for(int j = 0; j < dim; j++){
            T.set_vertex(j, remap[T.get_vertex(j)]);
            assert(T.get_vertex(j) < vtx_count);
        }
	}
    
	return vtx_count;
}


} // mesh
} // fastfem

#endif // MESH_TOOLS_HPP