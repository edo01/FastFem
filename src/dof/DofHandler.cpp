#include "FastFem/dof/DofHandler.hpp"

#include <iostream>

namespace fastfem {
namespace dof {

using global_vertex_id = fastfem::types::global_vertex_id;
using global_edge_id = fastfem::types::global_edge_id;
using global_face_id = fastfem::types::global_face_id;
using global_cell_id = fastfem::types::global_cell_id;

using global_dof_index_t = fastfem::types::global_dof_index_t;
using local_dof_index_t  = fastfem::types::local_dof_index_t;

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
        for(const global_vertex_id &v : T.get_vertex_indices()){
            if(vertex_dofs.find(v) == vertex_dofs.end()){
                vertex_dofs[v] = std::vector<long unsigned int>(dofs_per_vertex, -1);
                for(long unsigned int i = 0; i < dofs_per_vertex; ++i){
                    vertex_dofs[v][i] = n_dofs++;
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
            for(const global_edge_id &e : T.get_edges_indices()){
                if(edge_dofs.find(e) == edge_dofs.end()){
                    edge_dofs[e] = std::vector<long unsigned int>(dofs_per_edge, -1);
                    for(long unsigned int i = 0; i < dofs_per_edge; ++i){
                        edge_dofs[e][i] = n_dofs++;
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
            for(const global_face_id &f : T.get_faces_indices()){
                if(face_dofs.find(f) == face_dofs.end()){
                    face_dofs[f] = std::vector<long unsigned int>(dofs_per_face, -1);
                    for(long unsigned int i = 0; i < dofs_per_face; ++i){
                        face_dofs[f][i] = n_dofs++;
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
            for(const global_cell_id &c : T.get_cell_indices()){
                if(cell_dofs.find(c) == cell_dofs.end()){
                    cell_dofs[c] = std::vector<long unsigned int>(dofs_per_cell, -1);
                    for(long unsigned int i = 0; i < dofs_per_cell; ++i){
                        cell_dofs[c][i] = n_dofs++;
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
        size_t tag = it->first;
        // for each boundary tag, we store the dofs of the boundary in a set
        for(auto jt = it->second.begin(); jt != it->second.end(); ++jt){
            mesh::MeshSimplex<dim-1, spacedim> T = *jt;
            std::vector<global_dof_index_t> dofs = get_unordered_dofs_on_boundary(T);
            map_boundary_dofs[tag].insert(dofs.begin(), dofs.end());
        }
    }



    return n_dofs;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<global_dof_index_t> DoFHandler<dim, spacedim>::get_unordered_dofs_on_element(const mesh::MeshSimplex<dim, spacedim> &T) const {
    unsigned int dofs_per_element = (*fe).get_n_dofs_per_element();
    std::vector<global_dof_index_t> dofs;
    dofs.reserve(dofs_per_element);

    for(const global_vertex_id &v : T.get_vertex_indices()){
        const std::vector<global_dof_index_t> &d = vertex_dofs.at(v);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }
    if((*fe).get_n_dofs_per_edge() == 0) return dofs;
    
    for(const global_edge_id &e : T.get_edges_indices()){
        const std::vector<global_dof_index_t> &d = edge_dofs.at(e);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    if((*fe).get_n_dofs_per_face() == 0) return dofs;

    for(const global_face_id &f : T.get_faces_indices()){
        const std::vector<global_dof_index_t> &d = face_dofs.at(f);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    if((*fe).get_n_dofs_per_cell() == 0) return dofs;

    for(const global_cell_id &t : T.get_cell_indices()){
        const std::vector<global_dof_index_t> &d = cell_dofs.at(t);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    assert(dofs.size() == dofs_per_element);

    return dofs;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<global_dof_index_t> DoFHandler<dim, spacedim>::get_unordered_dofs_on_boundary(const mesh::MeshSimplex<dim-1, spacedim> &T) const {
    std::vector<global_dof_index_t> dofs;

    for(const global_vertex_id &v : T.get_vertex_indices()){
        const std::vector<global_dof_index_t> &d = vertex_dofs.at(v);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }
    if((*fe).get_n_dofs_per_edge() == 0) return dofs;
    
    for(const global_edge_id &e : T.get_edges_indices()){
        const std::vector<global_dof_index_t> &d = edge_dofs.at(e);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    if((*fe).get_n_dofs_per_face() == 0) return dofs;

    for(const global_face_id &f : T.get_faces_indices()){
        const std::vector<global_dof_index_t> &d = face_dofs.at(f);
        dofs.insert(dofs.end(), d.begin(), d.end());
    }

    return dofs;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<global_dof_index_t> DoFHandler<dim, spacedim>::get_ordered_dofs_on_element(const mesh::MeshSimplex<dim, spacedim> &T) const {
    unsigned int dofs_per_element = (*fe).get_n_dofs_per_element();
    std::vector<global_dof_index_t> dofs(dofs_per_element);

    for(const global_vertex_id &v : T.get_vertex_indices()){
        std::vector<local_dof_index_t> local_dofs = (*fe).get_local_dofs_on_subsimplex(T, v);
        for(int i = 0; i < local_dofs.size(); ++i){
            dofs[local_dofs[i]] = vertex_dofs.at(v)[i];
        }
    }

    if((*fe).get_n_dofs_per_edge() == 0) return dofs;

    for(const global_edge_id &e : T.get_edges_indices()){
        std::vector<local_dof_index_t> local_dofs = (*fe).get_local_dofs_on_subsimplex(T, e);
        for(int i = 0; i < local_dofs.size(); ++i){
            dofs[local_dofs[i]] = edge_dofs.at(e)[i];
        }
    }

    if((*fe).get_n_dofs_per_face() == 0) return dofs;

    for(const global_face_id &f : T.get_faces_indices()){
        std::vector<local_dof_index_t> local_dofs = (*fe).get_local_dofs_on_subsimplex(T, f);
        for(int i = 0; i < local_dofs.size(); ++i){
            dofs[local_dofs[i]] = face_dofs.at(f)[i];
        }
    }
    
    if((*fe).get_n_dofs_per_cell() == 0) return dofs;

    for(const global_cell_id &t : T.get_cell_indices()){
        std::vector<local_dof_index_t> local_dofs = (*fe).get_local_dofs_on_subsimplex(T, t);
        for(int i = 0; i < local_dofs.size(); ++i){
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