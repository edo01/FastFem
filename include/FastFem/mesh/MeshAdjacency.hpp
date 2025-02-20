/**
 * Utility class to handle the adjacency of the mesh elements.
 * 
 */

#ifndef MESHADJACENCY_HPP
#define MESHADJACENCY_HPP

#include <FastFem/common/hash_table.h>
#include <FastFem/types/CommonTypes.hpp>
#include <FastFem/mesh/Mesh.hpp>
#include <FastFem/mesh/Geometry.hpp>
#include <FastFem/mesh/VertexHasher.hpp>
#include <stdint.h>


namespace fastfem {
namespace mesh {
    template <unsigned int dim, unsigned int spacedim=dim>
    class MeshAdjacency
    {
    public:
        MeshAdjacency(Mesh<dim, spacedim> &mesh);

        std::vector<MeshSimplex<dim, spacedim>> get_adjacent_simplices(MeshSimplex<dim, spacedim> &T);

    private:
        Mesh<dim, spacedim> &mesh;
        VertexHasher<dim, spacedim> hasher;
        HashTable<uint64_t, std::vector<types::global_element_index>, VertexHasher<dim, spacedim>> vtx_adjacency;
    };

} // namespace mesh
} // namespace fastfem

#endif // MESHADJACENCY_HPP