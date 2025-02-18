#include "FastFem/mesh/MeshAdjacency.hpp"

namespace fastfem {
namespace mesh {

template <unsigned int dim, unsigned int spacedim>
MeshAdjacency<dim, spacedim>::MeshAdjacency(Mesh<dim, spacedim> &mesh) : mesh(mesh), hasher(mesh), vtx_adjacency{mesh.vtx_count(), hasher} {
    
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

template <unsigned int dim, unsigned int spacedim>
std::vector<MeshSimplex<dim, spacedim>> MeshAdjacency<dim, spacedim>::get_adjacent_simplices(MeshSimplex<dim, spacedim> &T) {
    std::vector<size_t> p_indices;
    for (unsigned int j = 0; j < dim + 1; j++) {
        std::vector<size_t> *adj_list = vtx_adjacency.get(T.get_vertex(j));
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

// Explicit template instantiation

template class MeshAdjacency<1, 1>;
template class MeshAdjacency<1, 2>;
template class MeshAdjacency<1, 3>;
template class MeshAdjacency<2, 2>;
template class MeshAdjacency<2, 3>;
template class MeshAdjacency<3, 3>;

} // namespace mesh
} // namespace fastfem