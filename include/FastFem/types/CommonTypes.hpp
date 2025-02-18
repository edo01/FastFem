#ifndef COMMONTYPES_HPP
#define COMMONTYPES_HPP

#include <array>
#include <unordered_map>
#include <functional>
#include <vector>


namespace fastfem{
namespace types{
    
typedef unsigned long global_dof_index_t;
typedef unsigned long local_dof_index_t;

template <unsigned short D>
using simplex_id = std::array<std::size_t, D+1>;

typedef simplex_id<0> global_vertex_id;
typedef simplex_id<1> global_edge_id;
typedef simplex_id<2> global_face_id;
typedef simplex_id<3> global_cell_id; 

typedef simplex_id<0> local_vertex_id;
typedef simplex_id<1> local_edge_id;
typedef simplex_id<2> local_face_id;
typedef simplex_id<3> local_cell_id; 


template <unsigned short D>
struct ArrayHasher {
    std::size_t operator()(const simplex_id<D>& arr) const {
        std::size_t h = 0;
        for (const auto& elem : arr) {
            h ^= std::hash<std::size_t>()(elem) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};
    
template <unsigned short D>
using global_dof_table = std::unordered_map<simplex_id<D>, std::vector<global_dof_index_t>, ArrayHasher<D>>;

template <unsigned short D>
using local_dof_table = std::unordered_map<simplex_id<D>, std::vector<local_dof_index_t>, ArrayHasher<D>>;

} // namespace types
} // namespace fastfem

#endif