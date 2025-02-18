#ifndef FASTFEM_DOFHANDLER_HPP
#define FASTFEM_DOFHANDLER_HPP

#include <vector>
#include <random>


namespace fastfem{
namespace dof{

template <unsigned int dim, unsigned int spacedim = dim>
class DofHandler
{
private:
    std::vector<std::vector<unsigned int>> dofs_per_cell;

public:
    DofHandler(unsigned int n_elems)
    {
        //dofs_per_cell.resize(n_elems);
        //rng
        // std::random_device rd;
        // std::mt19937 gen(rd());
        // std::uniform_int_distribution<unsigned int> dis(2, n_elems);
        // unsigned int n_dofs = dis(gen);
        // std::uniform_int_distribution<unsigned int> dis2(0, n_dofs - 1);
        // for(unsigned int i = 0; i < n_elems; ++i){
        //     dofs_per_cell[i].resize(n_dofs);
        //     for(unsigned int j = 0; j < n_dofs; ++j){
        //         dofs_per_cell[i][j] = dis2(gen);
        //     }
        // }

        // dofs_per_cell.resize(3);
        // dofs_per_cell[0] = {0,1};
        // dofs_per_cell[1] = {2,1};
        // dofs_per_cell[2] = {2,0};

        dofs_per_cell.resize(2);
        dofs_per_cell[0] = {0,1,3};
        dofs_per_cell[1] = {0,2,3};
    }
    
    inline unsigned int n_elements() const { return dofs_per_cell.size(); }
    inline const std::vector<unsigned int>& get_element_dofs(unsigned int i) const { return dofs_per_cell[i]; }
    
};

} // namespace dof
} // namespace fastfem

#endif // FASTFEM_DOFHANDLER_HPP