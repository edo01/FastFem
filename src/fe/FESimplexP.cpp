#include "FastFem/fe/FESimplexP.hpp"
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>

namespace fastfem {
namespace fe {

using global_vertex_index  = fastfem::types::global_vertex_index;
using global_edge_index    = fastfem::types::global_edge_index;
using global_face_index    = fastfem::types::global_face_index;
using global_cell_index    = fastfem::types::global_cell_index;

using local_vertex_index   = fastfem::types::local_vertex_index;
using local_edge_index     = fastfem::types::local_edge_index;
using local_face_index     = fastfem::types::local_face_index;
using local_cell_index     = fastfem::types::local_cell_index;

using local_dof_index = fastfem::types::local_dof_index;

template <unsigned int dim, unsigned int spacedim>
void FESimplexP<dim, spacedim>::
compute_stiffness_loc(const mesh::Simplex<dim, spacedim> &/*elem*/, linalg::FullMatrix &/*matrix*/) const {
    throw std::runtime_error("compute_stiffness_loc not implemented");
}

template <unsigned int dim, unsigned int spacedim>
std::vector<local_dof_index> FESimplexP<dim, spacedim>::
get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_vertex_index v) const {
    std::vector<local_dof_index> local_dofs_on_vertex(n_dofs_per_vertex);

    local_vertex_index local_vertex = map_global_simplex_to_local(T, v);
    const auto &dofs = vertex_dofs.at(local_vertex);
    std::copy(dofs.begin(), dofs.end(), local_dofs_on_vertex.begin());

    return local_dofs_on_vertex;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<local_dof_index> FESimplexP<dim, spacedim>::
get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_edge_index e) const {
    std::vector<local_dof_index> local_dofs_on_edge(n_dofs_per_edge);

    local_edge_index local_edge = map_global_simplex_to_local(T, e);
    const auto &dofs = edge_dofs.at(local_edge);
    std::copy(dofs.begin(), dofs.end(), local_dofs_on_edge.begin());

    return local_dofs_on_edge;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<local_dof_index> FESimplexP<dim, spacedim>::
get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_face_index f) const {
    std::vector<local_dof_index> local_dofs_on_face(n_dofs_per_face);

    local_face_index local_face = map_global_simplex_to_local(T, f);
    const auto &dofs = face_dofs.at(local_face);
    std::copy(dofs.begin(), dofs.end(), local_dofs_on_face.begin());

    return local_dofs_on_face;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<local_dof_index> FESimplexP<dim, spacedim>::
get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_cell_index c) const {
    std::vector<local_dof_index> local_dofs_on_cell(n_dofs_per_cell);

    local_cell_index local_cell = map_global_simplex_to_local(T, c);
    const auto &dofs = cell_dofs.at(local_cell);
    std::copy(dofs.begin(), dofs.end(), local_dofs_on_cell.begin());

    return local_dofs_on_cell;
}

template <unsigned int dim, unsigned int spacedim>
local_vertex_index FESimplexP<dim, spacedim>::
map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_vertex_index v) const {
    
    for (local_vertex_index i = 0; i < mesh::MeshSimplex<dim, spacedim>::n_vertices; ++i) {
        if (T.get_vertex(i) == v) return i;
    }
    assert(false);
}

template <unsigned int dim, unsigned int spacedim>
local_edge_index FESimplexP<dim, spacedim>::
map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_edge_index e) const {
    constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;
   
    for (local_vertex_index i = 0; i < n_vertices; ++i) {
        for (local_vertex_index j = i + 1; j < n_vertices; ++j) {

            global_edge_index e_simplex = {T.get_vertex(i), T.get_vertex(j)};
            std::sort(e_simplex.begin(), e_simplex.end());

            if (e_simplex == e) {
                return {i, j};
            }

        }
    }
    assert(false);
}

template <unsigned int dim, unsigned int spacedim>
local_face_index FESimplexP<dim, spacedim>::
map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_face_index f) const{
    constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;
    
    for (local_vertex_index i = 0; i < n_vertices; ++i) {
        for (local_vertex_index j = i + 1; j < n_vertices; ++j) {
            for (local_vertex_index k = j + 1; k < n_vertices; ++k) {
                global_face_index f_simplex = {T.get_vertex(i), T.get_vertex(j), T.get_vertex(k)};
                std::sort(f_simplex.begin(), f_simplex.end());
                
                if (f_simplex == f) {
                    return {i, j, k};
                }
            }
        }
    }
    assert(false);
}

template <unsigned int dim, unsigned int spacedim>
local_cell_index FESimplexP<dim, spacedim>::
map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_cell_index c) const {
    constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;
    for (local_vertex_index i = 0; i < n_vertices; ++i) {
        for (local_vertex_index j = i + 1; j < n_vertices; ++j) {
            for (local_vertex_index k = j + 1; k < n_vertices; ++k) {
                for (local_vertex_index l = k + 1; l < n_vertices; ++l) {
                    global_cell_index c_simplex = {T.get_vertex(i), T.get_vertex(j), T.get_vertex(k), T.get_vertex(l)};
                    std::sort(c_simplex.begin(), c_simplex.end());
                    if (c_simplex == c) {
                        return {i, j, k, l};
                    }
                }
            }
        }
    }
    assert(false);
}

template <unsigned int dim, unsigned int spacedim>
mesh::Point<spacedim> FESimplexP<dim, spacedim>::get_dof_coords(mesh::Simplex<dim, spacedim> T, local_dof_index dof) const 
{

    mesh::Point<spacedim> v0 = T.get_vertex(0);
    mesh::Point<spacedim> coords_ref = this->dofs[dof];
    mesh::Point<spacedim> mapped_coords;  // Initialize with v0 directly

    std::array<mesh::Point<spacedim>, dim + 1> vertices;
    for (local_vertex_index j = 0; j < dim + 1; ++j) {
        vertices[j] = T.get_vertex(j);  // Store vertices to avoid multiple calls
    }

    //this->

    for(unsigned int i = 0; i < spacedim; ++i) {
        mapped_coords[i] = v0[i];
        for(unsigned int j = 0; j < spacedim; ++j) {
            //std::cout << "(i, j): (" << i << ", " << j << ") = " << vertices[i + 1][j] - v0[j] << std::endl;
            mapped_coords[i] += (vertices[j + 1][i] - v0[i]) * coords_ref[j];
        }
    }

    return mapped_coords;

}

// Explicit instantiation
template class FESimplexP<1, 1>;
template class FESimplexP<1, 2>;
template class FESimplexP<1, 3>;
template class FESimplexP<2, 2>;
template class FESimplexP<2, 3>;
template class FESimplexP<3, 3>;

} // namespace fe
} // namespace fastfem