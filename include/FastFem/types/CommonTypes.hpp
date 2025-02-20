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
typedef unsigned long fastfem_index_t;
#else
typedef unsigned int fastfem_index_t;
#endif
    

typedef fastfem_index_t global_dof_index_t;
typedef unsigned short local_dof_index_t;

template <unsigned short D>
using global_simplex_id = std::array<fastfem_index_t, D+1>;

template <unsigned short D>
using local_simplex_id = std::array<local_dof_index_t, D+1>;

typedef global_simplex_id<0> global_vertex_id;
typedef global_simplex_id<1> global_edge_id;
typedef global_simplex_id<2> global_face_id;
typedef global_simplex_id<3> global_cell_id; 

typedef local_simplex_id<0> local_vertex_id;
typedef local_simplex_id<1> local_edge_id;
typedef local_simplex_id<2> local_face_id;
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
using global_dof_table = std::unordered_map<global_simplex_id<D>, std::vector<global_dof_index_t>, GlobalArrayHasher<D>>;

template <unsigned short D>
using local_dof_table = std::unordered_map<local_simplex_id<D>, std::vector<local_dof_index_t>, LocalArrayHasher<D>>;

} // namespace types
} // namespace fastfem

#endif