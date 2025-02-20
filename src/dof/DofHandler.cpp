#include "FastFem/dof/DofHandler.hpp"

#include <iostream>

namespace fastfem {
namespace dof {

using global_vertex_index = fastfem::types::global_vertex_index;
using global_edge_index = fastfem::types::global_edge_index;
using global_face_index = fastfem::types::global_face_index;
using global_cell_index = fastfem::types::global_cell_index;

using global_dof_index = fastfem::types::global_dof_index;
using local_dof_index  = fastfem::types::local_dof_index;
using global_element_index = fastfem::types::global_element_index;

template <unsigned int dim, unsigned int spacedim>
DoFHandler<dim, spacedim>::DoFHandler(const mesh::Mesh<dim, spacedim> &mesh, std::unique_ptr<fe::FESimplexP<dim, spacedim>> fe)
: mesh(mesh), fe(std::move(fe)) {}

template <unsigned int dim, unsigned int spacedim>
unsigned int DoFHandler<dim, spacedim>::distribute_dofs() {
    n_dofs = 0;
    
    /**
     * --------------------- VERTEX DOFS ---------------------
     */
    unsigned int dofs_per_vertex = (*fe).get_n_dofs_per_vertex();
    for(auto it = mesh.elem_begin(); it != mesh.elem_end(); ++it){
        mesh::MeshSimplex<dim, spacedim> T = *it;

        for(const global_vertex_index &v : T.get_vertex_indices()){
            // if the vertex is not already in the map, we add it
            auto [it, inserted] = vertex_dofs.emplace(v, std::vector<global_dof_index>(dofs_per_vertex, -1));
            // if the vertex was not already in the map, we assign the dofs to it
            if(inserted){
                for(long unsigned int i = 0; i < dofs_per_vertex; ++i){
                    it->second[i] = n_dofs++;
                }
            }  
        }   
    }
    
    /**
     * --------------------- EDGE DOFS ---------------------
     */
    unsigned int dofs_per_edge = (*fe).get_n_dofs_per_edge();
    if(dofs_per_edge > 0){
        for(auto it = mesh.elem_begin(); it != mesh.elem_end(); ++it){
            mesh::MeshSimplex<dim, spacedim> T = *it;

            for(const global_edge_index &e : T.get_edges_indices()){
                auto [it, inserted] = edge_dofs.emplace(e, std::vector<global_dof_index>(dofs_per_edge, -1));
                if(inserted){
                    for(long unsigned int i = 0; i < dofs_per_edge; ++i){
                        it->second[i] = n_dofs++;
                    }
                }
            }
        }
    }

    /**
     * --------------------- FACE DOFS ---------------------
     */
    unsigned int dofs_per_face = (*fe).get_n_dofs_per_face();
    if(dofs_per_face > 0){
        for(auto it = mesh.elem_begin(); it != mesh.elem_end(); ++it){
            mesh::MeshSimplex<dim, spacedim> T = *it;
            
            for(const global_face_index &f : T.get_faces_indices()){
                auto [it, inserted] = face_dofs.emplace(f, std::vector<global_dof_index>(dofs_per_face, -1));
                if(inserted){
                    for(long unsigned int i = 0; i < dofs_per_face; ++i){
                        it->second[i] = n_dofs++;
                    }
                }
            }
        }
    }

    /**
     * --------------------- CELL DOFS ---------------------
     */
    unsigned int dofs_per_cell = (*fe).get_n_dofs_per_cell();
    if(dofs_per_cell>0){
        for(auto it = mesh.elem_begin(); it != mesh.elem_end(); ++it){
            mesh::MeshSimplex<dim, spacedim> T = *it;
            
            for(const global_cell_index &c : T.get_cell_indices()){
                auto [it, inserted] = cell_dofs.emplace(c, std::vector<global_dof_index>(dofs_per_cell, -1));
                if(inserted){
                    for(long unsigned int i = 0; i < dofs_per_cell; ++i){
                        it->second[i] = n_dofs++;
                    }
                }
            }
        }
    }

    /*
     * We now store the dofs of the boundary. The boundary is partitioned in different tags, 
     * each tag is associated to a set of simplices. For each tag, we store the dofs of the boundary
     * in a set in order to avoid duplicates.
     */
    for(auto it = mesh.boundary_begin(); it != mesh.boundary_end(); ++it){
        types::boundary_index tag = it->first;
        // for each boundary tag, we store the dofs of the boundary in a set
        for(auto jt = it->second.begin(); jt != it->second.end(); ++jt){
            mesh::MeshSimplex<dim-1, spacedim> T = *jt;
            std::vector<global_dof_index> dofs = get_unordered_dofs_on_boundary(T);
            map_boundary_dofs[tag].insert(dofs.begin(), dofs.end());
        }
    }

    return n_dofs;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<global_dof_index> DoFHandler<dim, spacedim>::get_unordered_dofs_on_element(const mesh::MeshSimplex<dim, spacedim> &T) const {
    
    unsigned int dofs_per_element = (*fe).get_n_dofs_per_element();
    std::vector<global_dof_index> dofs;
    dofs.reserve(dofs_per_element);

    for(const global_vertex_index &v : T.get_vertex_indices()){
        const std::vector<global_dof_index> &d = vertex_dofs.at(v);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    if((*fe).get_n_dofs_per_edge() == 0) return dofs;
    for(const global_edge_index &e : T.get_edges_indices()){
        const std::vector<global_dof_index> &d = edge_dofs.at(e);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    if((*fe).get_n_dofs_per_face() == 0) return dofs;
    for(const global_face_index &f : T.get_faces_indices()){
        const std::vector<global_dof_index> &d = face_dofs.at(f);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    if((*fe).get_n_dofs_per_cell() == 0) return dofs;
    for(const global_cell_index &t : T.get_cell_indices()){
        const std::vector<global_dof_index> &d = cell_dofs.at(t);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    assert(dofs.size() == dofs_per_element);

    return dofs;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<global_dof_index> DoFHandler<dim, spacedim>::get_unordered_dofs_on_boundary(const mesh::MeshSimplex<dim-1, spacedim> &T) const {
    std::vector<global_dof_index> dofs;

    for(const global_vertex_index &v : T.get_vertex_indices()){
        const std::vector<global_dof_index> &d = vertex_dofs.at(v);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }
    if((*fe).get_n_dofs_per_edge() == 0) return dofs;
    
    for(const global_edge_index &e : T.get_edges_indices()){
        const std::vector<global_dof_index> &d = edge_dofs.at(e);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    if((*fe).get_n_dofs_per_face() == 0) return dofs;

    for(const global_face_index &f : T.get_faces_indices()){
        const std::vector<global_dof_index> &d = face_dofs.at(f);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    return dofs;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<global_dof_index> DoFHandler<dim, spacedim>::get_ordered_dofs_on_element(const mesh::MeshSimplex<dim, spacedim> &T) const {
    unsigned int dofs_per_element = (*fe).get_n_dofs_per_element();
    std::vector<global_dof_index> dofs(dofs_per_element);

    for(const global_vertex_index &v : T.get_vertex_indices()){
        std::vector<local_dof_index> local_dofs = (*fe).get_local_dofs_on_subsimplex(T, v);
        for(unsigned int i=0; i < local_dofs.size(); ++i){
            dofs[local_dofs[i]] = vertex_dofs.at(v)[i];
        }
    }

    if((*fe).get_n_dofs_per_edge() == 0) return dofs;

    for(const global_edge_index &e : T.get_edges_indices()){
        std::vector<local_dof_index> local_dofs = (*fe).get_local_dofs_on_subsimplex(T, e);
        for(unsigned int i=0; i < local_dofs.size(); ++i){
            dofs[local_dofs[i]] = edge_dofs.at(e)[i];
        }
    }

    if((*fe).get_n_dofs_per_face() == 0) return dofs;
    for(const global_face_index &f : T.get_faces_indices()){
        std::vector<local_dof_index> local_dofs = (*fe).get_local_dofs_on_subsimplex(T, f);
        for(unsigned int i=0; i < local_dofs.size(); ++i){
            dofs[local_dofs[i]] = face_dofs.at(f)[i];
        }
    }
    
    if((*fe).get_n_dofs_per_cell() == 0) return dofs;
    for(const global_cell_index &t : T.get_cell_indices()){
        std::vector<local_dof_index> local_dofs = (*fe).get_local_dofs_on_subsimplex(T, t);
        for(unsigned int i=0; i < local_dofs.size(); ++i){
            dofs[local_dofs[i]] = cell_dofs.at(t)[i];
        }
    }

    return dofs;
}

template <unsigned int dim, unsigned int spacedim>
void DoFHandler<dim, spacedim>::print_dofs() const {
    for(auto &v : vertex_dofs){
        std::cout << "Vertex: ";
        for(auto &d : v.second){
            std::cout << d << " ";
        }
        std::cout << std::endl;
    }
}


// Explicit instantiation
template class DoFHandler<1, 1>;
template class DoFHandler<1, 2>;
template class DoFHandler<1, 3>;
template class DoFHandler<2, 2>;
template class DoFHandler<2, 3>;
template class DoFHandler<3, 3>;

} // namespace dof
} // namespace fastfem