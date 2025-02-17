#ifndef VERTEXHASHER_HPP
#define VERTEXHASHER_HPP

#include <functional>
//#include <FastFem/mesh/Mesh.hpp>

namespace fastfem{
namespace mesh{

// Forward declaration of Mesh class template
template <unsigned int dim, unsigned int spacedim>
class Mesh;

// Forward declaration of Vertex struct template
template <unsigned int spacedim>
struct Vertex;

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

    VertexHasher(Mesh<dim, spacedim> &mesh) : mesh(mesh) {}

	static constexpr size_t empty_key = ~static_cast<size_t>(0); // reserved key to mark empty slots
	size_t hash(size_t key) const
	{
        std::size_t hash = 0;

        std::hash<double> hasher;
        std::hash<int> hasher_int;
        Vertex<spacedim> v = mesh.get_vertex(key);
        
        for(unsigned int i=0; i<spacedim; i++){
            std::size_t hvalue = hasher(v.point.coords[i]);
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

} // namespace mesh
} // namespace fastfem

#endif // VERTEXHASHER_HPP