#ifndef COMMONTYPES_HPP
#define COMMONTYPES_HPP

#include <array>
#include <unordered_map>
#include <functional>
#include <vector>

#define USE_LONG_INDEX

namespace fastfem{
namespace types{

/**
 * The type of the index used to identify the degrees of freedom. We can set it to be
 * either an unsigned int or an unsigned long. In this way, we can control the memory
 * usage of the application.
 */
#ifdef USE_LONG_INDEX
/**
 * Uniquely identifies a degree of freedom in the mesh
 */
typedef unsigned long global_dof_index;
/**
 * Uniquely identifies an element in the mesh
 */
typedef unsigned long global_element_index;
/**
 * Uniquely identifies a vertex in the mesh
 */
typedef unsigned long global_vertex_index;
#else
/**
 * Uniquely identifies a degree of freedom in the mesh
 */
typedef unsigned int global_dof_index;
/**
 * Uniquely identifies an element in the mesh
 */
typedef unsigned int global_element_index;
/**
 * Uniquely identifies a vertex in the mesh
 */
typedef unsigned int global_vertex_index;
#endif
    
/**
 * Uniquely identifies a degree of freedom on the reference simplex
 */
typedef unsigned int local_dof_index;
/**
 * Uniquely identifies a vertex on the reference simplex
 */
typedef unsigned int local_vertex_index;


template <unsigned short D>
using global_simplex_id = std::array<global_vertex_index, D+1>;

template <unsigned short D>
using local_simplex_id = std::array<local_vertex_index, D+1>;

/**
 * A couple of ORDERED global_vertex_id that uniquely identifies an edge in the mesh
 */
typedef global_simplex_id<1> global_edge_id;
/**
 * A triplet of ORDERED global_vertex_id that uniquely identifies a face in the mesh
 */
typedef global_simplex_id<2> global_face_id;
/**
 * A quadruplet of ORDERED global_vertex_id that uniquely identifies a cell in the mesh
 */
typedef global_simplex_id<3> global_cell_id; 

/**
 * A couple of ORDERED local_vertex_id that uniquely identifies an edge in the reference simplex
 */
typedef local_simplex_id<1> local_edge_id;
/**
 * A triplet of ORDERED local_vertex_id that uniquely identifies a face in the reference simplex
 */
typedef local_simplex_id<2> local_face_id;
/**
 * A quadruplet of ORDERED local_vertex_id that uniquely identifies a cell in the reference simplex
 */
typedef local_simplex_id<3> local_cell_id; 


template <unsigned short D>
struct GlobalArrayHasher {
    std::size_t operator()(const global_simplex_id<D>& arr) const {
        std::size_t h = 0;
        for (const auto& elem : arr) {
            h ^= std::hash<std::size_t>()(elem) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

template <unsigned short D>
struct LocalArrayHasher {
    std::size_t operator()(const local_simplex_id<D>& arr) const {
        std::size_t h = 0;
        for (const auto& elem : arr) {
            h ^= std::hash<std::size_t>()(elem) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

template <unsigned short D>
using global_dof_table = std::unordered_map<global_simplex_id<D>, std::vector<global_dof_index>, GlobalArrayHasher<D>>;

template <unsigned short D>
using local_dof_table = std::unordered_map<local_simplex_id<D>, std::vector<local_dof_index>, LocalArrayHasher<D>>;

typedef std::unordered_map<global_vertex_index, std::vector<global_dof_index>, std::hash<global_vertex_index>> global_vertex_dof_table;
typedef global_dof_table<1> global_edge_dof_table;
typedef global_dof_table<2> global_face_dof_table;
typedef global_dof_table<3> global_cell_dof_table;

typedef std::unordered_map<local_vertex_index, std::vector<local_dof_index>, std::hash<local_vertex_index>> local_vertex_dof_table;
typedef local_dof_table<1> local_edge_dof_table;
typedef local_dof_table<2> local_face_dof_table;
typedef local_dof_table<3> local_cell_dof_table;

} // namespace types
} // namespace fastfem

#endif