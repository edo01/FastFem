#include "mesh/mesh_tools.hpp"

namespace mesh
{

/*	
 *	The Hasher class must implement 'hash', 'is_empty', 'is_equal'
 *	functions, and reserve a key named empty_key for marking empty 
 *	slots.  
 */
struct VertexHasher 
{
	const struct Vertex *pos;
	static constexpr uint32_t empty_key = ~static_cast<uint32_t>(0); // reserved key to mark empty slots
	size_t hash(uint32_t key) const
	{
		uint32_t hash = 0;
		/* reinterpret the pointer to the vertex as a pointer to an array of 3 uint32_t values
		   in order to pass it to the murmur2_32 function */
		const uint32_t *p =
		    reinterpret_cast<const uint32_t *>(pos + key);
		/*
		 * hash the three coordinates
		 */
		hash = murmur2_32(hash, p[0]);
		hash = murmur2_32(hash, p[1]);
		hash = murmur2_32(hash, p[2]);
		return hash;
	}

	bool is_empty(uint32_t key) const { return (key == empty_key); }

	bool is_equal(uint32_t key1, uint32_t key2) const
	{
		return pos[key1] == pos[key2]; // redefinition of the equality operator
	}
};

int dedup_mesh_vertices(struct Mesh *m)
{
	int vtx_count = 0;
	int V = m->vtx_count;
	int *remap = (int *)malloc(V * sizeof(int));

	/*
	 * We have replaced the linear search (O(n^2)) with an
	 * hash table to deduplicate the vertices. 
	 */
	
	/*
	 * 	To build the hasher we tried to follow the design of the professor
	 *  using a custom Hasher struct. 
	 *
	 *	Our goal is to remap the vertices. To achieve that we have to build an 
	 *  hash table that given the index of a vertex returns the deduplicated index.
	 *  To achieve that we have to build an Hasher struct which is able to calc 
	 *  the hash of the vertex given the key (index inside the vertices array) and 
	 *  also being able to recognize when two different vertices have the same hash.
	 *
	 *  So the hash table will be built using the index of the vertex as key and the
	 *  deduplicated index as value. 
	 *
	 */
	VertexHasher hasher{m->vertices};
	HashTable<uint32_t, uint32_t, VertexHasher> vtx_remap(V, hasher);

	for (size_t i = 0; i < V; ++i) {
		uint32_t *p = vtx_remap.get_or_set(i, vtx_count);
		if (p) {
			remap[i] = *p;
		} else {
			remap[i] = vtx_count;
			vtx_count++;
		}
	}

	/* Remap vertices */
	for (int i = 0; i < V; i++) {
		m->vertices[remap[i]] = m->vertices[i];
	}
	/* Remap triangle indices */
	for (int i = 0; i < m->tri_count; i++) {
		struct Triangle *T = &m->triangles[i];
		T->a = remap[T->a];
		assert(T->a < vtx_count);
		T->b = remap[T->b];
		assert(T->b < vtx_count);
		T->c = remap[T->c];
		assert(T->c < vtx_count);
	}
	free(remap);
	return vtx_count;
}

} // mesh