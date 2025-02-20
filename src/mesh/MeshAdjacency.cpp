#include "FastFem/mesh/MeshAdjacency.hpp"

namespace fastfem {
namespace mesh {

using global_element_index = fastfem::types::global_element_index;
using global_vertex_index  = fastfem::types::global_vertex_index;
using local_vertex_index   = fastfem::types::local_vertex_index;


template <unsigned int dim, unsigned int spacedim>
MeshAdjacency<dim, spacedim>::MeshAdjacency(Mesh<dim, spacedim> &mesh) : mesh(mesh), hasher(mesh), vtx_adjacency{mesh.vtx_count(), hasher} {
    
    for(global_element_index i=0; i<mesh.elem_count(); i++){
        MeshSimplex<dim, spacedim> T = mesh.get_mesh_element(i);
        
        for(local_vertex_index j=0; j<dim+1; j++){
            global_vertex_index vtx = T.get_vertex(j);
            std::vector<global_element_index> new_list;
            new_list.push_back(i);
            std::vector<size_t> *p = vtx_adjacency.get_or_set(vtx, new_list);
            if(p != nullptr){
                p->push_back(i);
            }
        }
    }
}

template <unsigned int dim, unsigned int spacedim>
std::vector<MeshSimplex<dim, spacedim>> MeshAdjacency<dim, spacedim>::get_adjacent_simplices(MeshSimplex<dim, spacedim> &T) {
    std::vector<global_element_index> p_indices;
    for (local_vertex_index j = 0; j < dim + 1; j++) {
        std::vector<global_element_index> *adj_list = vtx_adjacency.get(T.get_vertex(j));
        
        for (auto it = adj_list->begin(); it != adj_list->end(); ++it) {
            // avoid duplicates
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

// Explicit template instantiation

template class MeshAdjacency<1, 1>;
template class MeshAdjacency<1, 2>;
template class MeshAdjacency<1, 3>;
template class MeshAdjacency<2, 2>;
template class MeshAdjacency<2, 3>;
template class MeshAdjacency<3, 3>;

} // namespace mesh
} // namespace fastfem