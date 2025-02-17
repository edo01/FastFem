/**
 * Utility class to handle the adjacency of the mesh elements.
 * 
 */

#ifndef MESHADJACENCY_HPP
#define MESHADJACENCY_HPP

#include <FastFem/common/hash_table.h>
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
    MeshAdjacency(Mesh<dim, spacedim> &mesh) : mesh(mesh), hasher(mesh), vtx_adjacency{mesh.vtx_count(), hasher} {

        for(size_t i=0; i<mesh.elem_count(); i++){
            MeshSimplex<dim, spacedim> T = mesh.get_mesh_element(i);
            for(unsigned int j=0; j<dim+1; j++){
                size_t vtx = T.get_vertex(j);
                std::vector<size_t> new_list;
                new_list.push_back(i);
                std::vector<size_t> *p = vtx_adjacency.get_or_set(vtx, new_list);
                if(p != nullptr){
                    p->push_back(i);
                }
            }
        }
    }

    std::vector<MeshSimplex<dim, spacedim>> get_adjacent_simplices(MeshSimplex<dim, spacedim> &T) {
        std::vector<size_t> p_indices;
        for (int j = 0; j < dim + 1; j++) {
            std::vector<size_t> *adj_list = vtx_adjacency.get(T.get_vertex(j));
            // add only the elements that are not already in the list
            for (auto it = adj_list->begin(); it != adj_list->end(); ++it) {
                if (std::find(p_indices.begin(), p_indices.end(), *it) == p_indices.end()) {
                    p_indices.push_back(*it);
                }
            }
        }
        std::vector<MeshSimplex<dim, spacedim>> p;
        for (auto it = p_indices.begin(); it != p_indices.end(); ++it) {
            if(!(mesh.get_mesh_element(*it) == T))
                p.push_back(mesh.get_mesh_element(*it));
        }
        return p;
    }

    void print_adjacency() {
        for (size_t i = 0; i < mesh.vtx_count(); ++i) {
            const std::vector<size_t> *p = vtx_adjacency.get(i);
            std::cout << "Vertex " << i << " is adjacent to: ";
            if (p != nullptr) {
                for (size_t j = 0; j < p->size(); ++j) {
                    std::cout << (*p)[j] << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    ~MeshAdjacency() {}

private:
    Mesh<dim, spacedim> &mesh;
    VertexHasher<dim, spacedim> hasher;
    HashTable<uint64_t, std::vector<size_t>, VertexHasher<dim, spacedim>> vtx_adjacency;
};

} // namespace mesh
} // namespace fastfem

#endif // MESHADJACENCY_HPP