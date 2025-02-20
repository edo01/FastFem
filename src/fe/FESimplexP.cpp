#include "FastFem/fe/FESimplexP.hpp"
#include <vector>
#include <algorithm>
#include <cassert>

namespace fastfem {
namespace fe {


using global_vertex_id = fastfem::types::global_vertex_id;
using global_edge_id = fastfem::types::global_edge_id;
using global_face_id = fastfem::types::global_face_id;
using global_cell_id = fastfem::types::global_cell_id;

using local_vertex_id = fastfem::types::local_vertex_id;
using local_edge_id   = fastfem::types::local_edge_id;
using local_face_id   = fastfem::types::local_face_id;
using local_cell_id   = fastfem::types::local_cell_id; 

using local_dof_index_t = fastfem::types::local_dof_index_t;


template <unsigned int dim, unsigned int spacedim>
std::vector<local_dof_index_t> FESimplexP<dim, spacedim>::
get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_vertex_id v) const {
    std::vector<local_dof_index_t> local_dofs_on_vertex(n_dofs_per_vertex);

    local_vertex_id local_vertex = map_global_simplex_to_local(T, v);
    const auto &dofs = vertex_dofs.at(local_vertex);
    std::copy(dofs.begin(), dofs.end(), local_dofs_on_vertex.begin());

    return local_dofs_on_vertex;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<local_dof_index_t> FESimplexP<dim, spacedim>::
get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_edge_id e) const {
    std::vector<local_dof_index_t> local_dofs_on_edge(n_dofs_per_edge);

    local_edge_id local_edge = map_global_simplex_to_local(T, e);
    const auto &dofs = edge_dofs.at(local_edge);
    std::copy(dofs.begin(), dofs.end(), local_dofs_on_edge.begin());

    return local_dofs_on_edge;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<local_dof_index_t> FESimplexP<dim, spacedim>::
get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_face_id f) const {
    std::vector<local_dof_index_t> local_dofs_on_face(n_dofs_per_face);

    local_face_id local_face = map_global_simplex_to_local(T, f);
    const auto &dofs = face_dofs.at(local_face);
    std::copy(dofs.begin(), dofs.end(), local_dofs_on_face.begin());

    return local_dofs_on_face;
}

template <unsigned int dim, unsigned int spacedim>
std::vector<local_dof_index_t> FESimplexP<dim, spacedim>::
get_local_dofs_on_subsimplex(const mesh::MeshSimplex<dim, spacedim> &T, global_cell_id c) const {
    std::vector<local_dof_index_t> local_dofs_on_cell(n_dofs_per_cell);

    local_cell_id local_cell = map_global_simplex_to_local(T, c);
    const auto &dofs = cell_dofs.at(local_cell);
    std::copy(dofs.begin(), dofs.end(), local_dofs_on_cell.begin());

    return local_dofs_on_cell;
}

template <unsigned int dim, unsigned int spacedim>
local_vertex_id FESimplexP<dim, spacedim>::
map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_vertex_id v) const {
    constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;
    for (unsigned int i = 0; i < n_vertices; ++i) {
        if (T.get_vertex(i) == v[0]) return {i};
    }
    assert(false);
}

template <unsigned int dim, unsigned int spacedim>
local_edge_id FESimplexP<dim, spacedim>::
map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_edge_id e) const {
    constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;
    for (int i = 0; i < n_vertices; ++i) {
        for (int j = i + 1; j < n_vertices; ++j) {
            global_edge_id e_simplex = {T.get_vertex(i), T.get_vertex(j)};
            std::sort(e_simplex.begin(), e_simplex.end());
            if (e_simplex == e) {
                return {i, j};
            }
        }
    }
    assert(false);
}

template <unsigned int dim, unsigned int spacedim>
local_face_id FESimplexP<dim, spacedim>::
map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_face_id f) const{
    constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;
    for (int i = 0; i < n_vertices; ++i) {
        for (int j = i + 1; j < n_vertices; ++j) {
            for (int k = j + 1; k < n_vertices; ++k) {
                global_face_id f_simplex = {T.get_vertex(i), T.get_vertex(j), T.get_vertex(k)};
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
local_cell_id FESimplexP<dim, spacedim>::
map_global_simplex_to_local(const mesh::MeshSimplex<dim, spacedim> &T, global_cell_id c) const {
    constexpr unsigned int n_vertices = mesh::MeshSimplex<dim, spacedim>::n_vertices;
    for (int i = 0; i < n_vertices; ++i) {
        for (int j = i + 1; j < n_vertices; ++j) {
            for (int k = j + 1; k < n_vertices; ++k) {
                for (int l = k + 1; l < n_vertices; ++l) {
                    global_cell_id c_simplex = {T.get_vertex(i), T.get_vertex(j), T.get_vertex(k), T.get_vertex(l)};
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


// Explicit instantiation
template class FESimplexP<1, 1>;
template class FESimplexP<1, 2>;
template class FESimplexP<1, 3>;
template class FESimplexP<2, 2>;
template class FESimplexP<2, 3>;
template class FESimplexP<3, 3>;

} // namespace fe
} // namespace fastfem