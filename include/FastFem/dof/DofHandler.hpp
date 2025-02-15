/**
 * @file DofHandler.hpp
 * 
 * The idea is the following: we need a class that is able to handle the degrees of freedom (DoFs) of a finite element space. 
 * This class should be able to distribute the DoFs on the mesh, to extract the DoFs on the boundary of the domain.
 * 
 * Before going into the details of the implementation, let's see how the DoFHandler class is used in the MinSurFEM code:
 * We want avoid to have a double loop in the assembly of the linear system O(n_dofs^2), since most of these couples will be zero.
 * 
 * 
 * So the dof handler must expose a method that allows to iterate over the cells of the mesh, compute the local matrix and vector, 
 * and then add them to the global matrix and vector. Something like deal.II does:
 * for (auto &cell : dof_handler.cell_iterators())
 * {
 *  //  1. fill a full matrix of size dofs_per_cell x dofs_per_cell and a vector of size dofs_per_cell
 *  //    The difference with the deal.II code is that we do not provide a class that computes the values of the basis functions
 *  //    and their gradients at the quadrature points. The user should compute them and build the local matrix. 
 * 
 *  //  2. get the global indices of the DoFs of the current cell
 *  //  3. add the local matrix and vector to the global matrix and vector
 * }
 * 
 * In order to do this, the DoFHandler class must store for each cell the global indices of the DoFs of that cell.
 * To do this, we must remember that each DOF may be associated with a vertex, an edge, a face, or a cell of the mesh. A DOF on the vertex has support on the cells that contain that vertex, and so on. 
 * 
 * Where the information about the DoFs in each single cell is stored? Given a cell, we must be able to retrieve the points of the mesh that are associated with the DoFs of that cell. 
 * 
 */