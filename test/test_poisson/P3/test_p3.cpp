#include <iostream>

#include "FastFem/mesh/Mesh.hpp"
#include "FastFem/mesh/MeshMaker.hpp"

#include "FastFem/fe/FESimplexP.hpp"
#include "FastFem/dof/DofHandler.hpp"

#include "FastFem/linalg/Vector.hpp"
#include "FastFem/linalg/sparseMatrices/CSRMatrix.hpp"
#include "FastFem/linalg/MatrixTools.hpp"

#include "FastFem/linalg/iterativeSolvers/CGSolver.hpp"

#include "FastFem/mesh/MeshIO.hpp"

#include "FastFem/types/CommonTypes.hpp"

#include "FastFem/mesh/MeshIO.hpp"

using namespace fastfem;

#define SQUARE(x) ((x) * (x))

double compute_den(double xa, double ya, double xb, double yb, double xc, double yc) {
    return SQUARE(xc * (-ya + yb) + xb * (ya - yc) + xa * (-yb + yc));
}

void compute_stiffness_loc(mesh::Simplex<2, 2> triangle, double matrix[10][10]) {

    double xa = triangle.get_vertex(0).coords[0];
    double ya = triangle.get_vertex(0).coords[1];
    double xb = triangle.get_vertex(1).coords[0];
    double yb = triangle.get_vertex(1).coords[1];
    double xc = triangle.get_vertex(2).coords[0];
    double yc = triangle.get_vertex(2).coords[1];

/*     matrix[0][0] = (17 * (SQUARE(xb - xc) + SQUARE(yb - yc))) / 40.0;
    matrix[0][1] = (3 * (SQUARE(xb) - 19 * xa * (xb - xc) + 17 * xb * xc - 18 * SQUARE(xc) - 19 * ya * yb + SQUARE(yb) + 19 * ya * yc + 17 * yb * yc - 18 * SQUARE(yc))) / 80.0;
    matrix[0][2] = (3 * (SQUARE(xb) + 8 * xa * (xb - xc) - 10 * xb * xc + 9 * SQUARE(xc) + 8 * ya * yb + SQUARE(yb) - 8 * ya * yc - 10 * yb * yc + 9 * SQUARE(yc))) / 80.0;
    matrix[0][3] = (-7 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[0][4] = (-3 * (SQUARE(xb - xc) + SQUARE(yb - yc))) / 80.0;
    matrix[0][5] = (-3 * (SQUARE(xb) - 2 * xb * xc + SQUARE(xc) + SQUARE(yb - yc))) / 80.0;
    matrix[0][6] = (7 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix[0][7] = (3 * (-8 * xa * xb + 9 * SQUARE(xb) + 8 * xa * xc - 10 * xb * xc + SQUARE(xc) - 8 * ya * yb + 9 * SQUARE(yb) + 8 * ya * yc - 10 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[0][8] = (3 * (-18 * SQUARE(xb) + 19 * xa * (xb - xc) + 17 * xb * xc + SQUARE(xc) + 19 * ya * yb - 18 * SQUARE(yb) - 19 * ya * yc + 17 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[0][9] = 0;

    matrix[1][0] = (3 * (SQUARE(xb) - 19 * xa * (xb - xc) + 17 * xb * xc - 18 * SQUARE(xc) - 19 * ya * yb + SQUARE(yb) + 19 * ya * yc + 17 * yb * yc - 18 * SQUARE(yc))) / 80.0;
    matrix[1][1] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[1][2] = (-27 * (SQUARE(xa) + SQUARE(xb) + 2 * xa * (xb - 2 * xc) - 4 * xb * xc + 4 * SQUARE(xc) + SQUARE(ya) + 2 * ya * yb + SQUARE(yb) - 4 * ya * yc - 4 * yb * yc + 4 * SQUARE(yc))) / 80.0;
    matrix[1][3] = (3 * (SQUARE(xa) + 8 * xa * xb - 10 * xa * xc - 8 * xb * xc + 9 * SQUARE(xc) + SQUARE(ya) + 8 * ya * yb - 10 * ya * yc - 8 * yb * yc + 9 * SQUARE(yc))) / 80.0;
    matrix[1][4] = (-27 * (-SQUARE(xb) + xa * (xb - xc) + xb * xc + (ya - yb) * (yb - yc))) / 80.0;
    matrix[1][5] = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix[1][6] = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80.0;
    matrix[1][7] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[1][8] = (-27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 16.0;
    matrix[1][9] = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40.0;

    matrix[2][0] = (3 * (SQUARE(xb) + 8 * xa * (xb - xc) - 10 * xb * xc + 9 * SQUARE(xc) + 8 * ya * yb + SQUARE(yb) - 8 * ya * yc - 10 * yb * yc + 9 * SQUARE(yc))) / 80.0;
    matrix[2][1] = (-27 * (SQUARE(xa) + SQUARE(xb) + 2 * xa * (xb - 2 * xc) - 4 * xb * xc + 4 * SQUARE(xc) + SQUARE(ya) + 2 * ya * yb + SQUARE(yb) - 4 * ya * yc - 4 * yb * yc + 4 * SQUARE(yc))) / 80.0;
    matrix[2][2] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[2][3] = (3 * (SQUARE(xa) - 19 * xa * xb + 17 * xa * xc + 19 * xb * xc - 18 * SQUARE(xc) + SQUARE(ya) - 19 * ya * yb + 17 * ya * yc + 19 * yb * yc - 18 * SQUARE(yc))) / 80.0;
    matrix[2][4] = (27 * (-SQUARE(xb) + xa * (xb - xc) + xb * xc + (ya - yb) * (yb - yc))) / 16.0;
    matrix[2][5] = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix[2][6] = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80.0;
    matrix[2][7] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[2][8] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[2][9] = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40.0;

    matrix[3][0] = (-7 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[3][1] = (3 * (SQUARE(xa) + 8 * xa * xb - 10 * xa * xc - 8 * xb * xc + 9 * SQUARE(xc) + SQUARE(ya) + 8 * ya * yb - 10 * ya * yc - 8 * yb * yc + 9 * SQUARE(yc))) / 80.0;
    matrix[3][2] = (3 * (SQUARE(xa) - 19 * xa * xb + 17 * xa * xc + 19 * xb * xc - 18 * SQUARE(xc) + SQUARE(ya) - 19 * ya * yb + 17 * ya * yc + 19 * yb * yc - 18 * SQUARE(yc))) / 80.0;
    matrix[3][3] = (17 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 40.0;
    matrix[3][4] = (-3 * (18 * SQUARE(xa) - 19 * xa * xb - 17 * xa * xc + 19 * xb * xc - SQUARE(xc) + 18 * SQUARE(ya) - 19 * ya * yb - 17 * ya * yc + 19 * yb * yc - SQUARE(yc))) / 80.0;
    matrix[3][5] = (3 * (9 * SQUARE(xa) + 8 * xb * xc + SQUARE(xc) - 2 * xa * (4 * xb + 5 * xc) + 9 * SQUARE(ya) - 8 * ya * yb - 10 * ya * yc + 8 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[3][6] = (-7 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[3][7] = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80.0;
    matrix[3][8] = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80.0;
    matrix[3][9] = 0;

    matrix[4][0] = (-3 * (SQUARE(xb - xc) + SQUARE(yb - yc))) / 80.0;
    matrix[4][1] = (-27 * (-SQUARE(xb) + xa * (xb - xc) + xb * xc + (ya - yb) * (yb - yc))) / 80.0;
    matrix[4][2] = (27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 16.0;
    matrix[4][3] = (-3 * (18 * SQUARE(xa) - 19 * xa * xb - 17 * xa * xc + 19 * xb * xc - SQUARE(xc) + 18 * SQUARE(ya) - 19 * ya * yb - 17 * ya * yc + 19 * yb * yc - SQUARE(yc))) / 80.0;
    matrix[4][4] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[4][5] = (-27 * (4 * SQUARE(xa) + SQUARE(xb) + 2 * xb * xc + SQUARE(xc) - 4 * xa * (xb + xc) + 4 * SQUARE(ya) - 4 * ya * yb + SQUARE(yb) - 4 * ya * yc + 2 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[4][6] = (3 * (9 * SQUARE(xa) + SQUARE(xb) + 8 * xb * xc - 2 * xa * (5 * xb + 4 * xc) + (ya - yb) * (9 * ya - yb - 8 * yc))) / 80.0;
    matrix[4][7] = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[4][8] = (27 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 80.0;
    matrix[4][9] = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40.0;

    matrix[5][0] = (-3 * (SQUARE(xb) - 2 * xb * xc + SQUARE(xc) + SQUARE(yb - yc))) / 80.0;
    matrix[5][1] = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix[5][2] = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix[5][3] = (3 * (9 * SQUARE(xa) + 8 * xb * xc + SQUARE(xc) - 2 * xa * (4 * xb + 5 * xc) + 9 * SQUARE(ya) - 8 * ya * yb - 10 * ya * yc + 8 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[5][4] = (-27 * (4 * SQUARE(xa) + SQUARE(xb) + 2 * xb * xc + SQUARE(xc) - 4 * xa * (xb + xc) + 4 * SQUARE(ya) - 4 * ya * yb + SQUARE(yb) - 4 * ya * yc + 2 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[5][5] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[5][6] = (-3 * (xa - xb) * (18 * xa + xb - 19 * xc) - 3 * (ya - yb) * (18 * ya + yb - 19 * yc)) / 80.0;
    matrix[5][7] = (-27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 16.0;
    matrix[5][8] = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[5][9] = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40.0;

    matrix[6][0] = (7 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix[6][1] = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80.0;
    matrix[6][2] = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80.0;
    matrix[6][3] = (-7 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[6][4] = (3 * (9 * SQUARE(xa) + SQUARE(xb) + 8 * xb * xc - 2 * xa * (5 * xb + 4 * xc) + (ya - yb) * (9 * ya - yb - 8 * yc))) / 80.0;
    matrix[6][5] = (-3 * (xa - xb) * (18 * xa + xb - 19 * xc) - 3 * (ya - yb) * (18 * ya + yb - 19 * yc)) / 80.0;
    matrix[6][6] = (17 * (SQUARE(xa - xb) + SQUARE(ya - yb))) / 40.0;
    matrix[6][7] = (3 * (SQUARE(xa) - 18 * SQUARE(xb) + xa * (17 * xb - 19 * xc) + 19 * xb * xc + (ya - yb) * (ya + 18 * yb - 19 * yc))) / 80.0;
    matrix[6][8] = (3 * (SQUARE(xa) + 9 * SQUARE(xb) - 8 * xb * xc + xa * (-10 * xb + 8 * xc) + (ya - yb) * (ya - 9 * yb + 8 * yc))) / 80.0;
    matrix[6][9] = 0;

    matrix[7][0] = (3 * (-8 * xa * xb + 9 * SQUARE(xb) + 8 * xa * xc - 10 * xb * xc + SQUARE(xc) - 8 * ya * yb + 9 * SQUARE(yb) + 8 * ya * yc - 10 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[7][1] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[7][2] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[7][3] = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80.0;
    matrix[7][4] = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[7][5] = (-27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 16.0;
    matrix[7][6] = (3 * (SQUARE(xa) - 18 * SQUARE(xb) + xa * (17 * xb - 19 * xc) + 19 * xb * xc + (ya - yb) * (ya + 18 * yb - 19 * yc))) / 80.0;
    matrix[7][7] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[7][8] = (-27 * (SQUARE(xa - 2 * xb + xc) + SQUARE(ya - 2 * yb + yc))) / 80.0;
    matrix[7][9] = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40.0;

    matrix[8][0] = (3 * (-18 * SQUARE(xb) + 19 * xa * (xb - xc) + 17 * xb * xc + SQUARE(xc) + 19 * ya * yb - 18 * SQUARE(yb) - 19 * ya * yc + 17 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[8][1] = (-27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 16.0;
    matrix[8][2] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[8][3] = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80.0;
    matrix[8][4] = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[8][5] = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[8][6] = (3 * (SQUARE(xa) + 9 * SQUARE(xb) - 8 * xb * xc + xa * (-10 * xb + 8 * xc) + (ya - yb) * (ya - 9 * yb + 8 * yc))) / 80.0;
    matrix[8][7] = (-27 * (SQUARE(xa - 2 * xb + xc) + SQUARE(ya - 2 * yb + yc))) / 80.0;
    matrix[8][8] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[8][9] = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40.0;

    matrix[9][0] = 0;
    matrix[9][1] = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40.0;
    matrix[9][2] = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40.0;
    matrix[9][3] = 0;
    matrix[9][4] = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40.0;
    matrix[9][5] = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40.0;
    matrix[9][6] = 0;
    matrix[9][7] = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40.0;
    matrix[9][8] = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40.0;
    matrix[9][9] = (81 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 20.0;
 */

    matrix[0][0] = (17 * (SQUARE(xb - xc) + SQUARE(yb - yc))) / 40.0;
    matrix[0][1] = (3 * (SQUARE(xb) - 19 * xa * (xb - xc) + 17 * xb * xc - 18 * SQUARE(xc) - 19 * ya * yb + SQUARE(yb) + 19 * ya * yc + 17 * yb * yc - 18 * SQUARE(yc))) / 80.0;
    matrix[0][2] = (3 * (SQUARE(xb) + 8 * xa * (xb - xc) - 10 * xb * xc + 9 * SQUARE(xc) + 8 * ya * yb + SQUARE(yb) - 8 * ya * yc - 10 * yb * yc + 9 * SQUARE(yc))) / 80.0;
    matrix[0][3] = (-7 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[0][4] = (-3 * (SQUARE(xb - xc) + SQUARE(yb - yc))) / 80.0;
    matrix[0][5] = (-3 * (SQUARE(xb) - 2 * xb * xc + SQUARE(xc) + SQUARE(yb - yc))) / 80.0;
    matrix[0][6] = (7 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix[0][7] = (3 * (-8 * xa * xb + 9 * SQUARE(xb) + 8 * xa * xc - 10 * xb * xc + SQUARE(xc) - 8 * ya * yb + 9 * SQUARE(yb) + 8 * ya * yc - 10 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[0][8] = (3 * (-18 * SQUARE(xb) + 19 * xa * (xb - xc) + 17 * xb * xc + SQUARE(xc) + 19 * ya * yb - 18 * SQUARE(yb) - 19 * ya * yc + 17 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[0][9] = 0;

    matrix[1][1] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[1][2] = (-27 * (SQUARE(xa) + SQUARE(xb) + 2 * xa * (xb - 2 * xc) - 4 * xb * xc + 4 * SQUARE(xc) + SQUARE(ya) + 2 * ya * yb + SQUARE(yb) - 4 * ya * yc - 4 * yb * yc + 4 * SQUARE(yc))) / 80.0;
    matrix[1][3] = (3 * (SQUARE(xa) + 8 * xa * xb - 10 * xa * xc - 8 * xb * xc + 9 * SQUARE(xc) + SQUARE(ya) + 8 * ya * yb - 10 * ya * yc - 8 * yb * yc + 9 * SQUARE(yc))) / 80.0;
    matrix[1][4] = (-27 * (-SQUARE(xb) + xa * (xb - xc) + xb * xc + (ya - yb) * (yb - yc))) / 80.0;
    matrix[1][5] = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix[1][6] = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80.0;
    matrix[1][7] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[1][8] = (-27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 16.0;
    matrix[1][9] = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40.0;

    matrix[2][2] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[2][3] = (3 * (SQUARE(xa) - 19 * xa * xb + 17 * xa * xc + 19 * xb * xc - 18 * SQUARE(xc) + SQUARE(ya) - 19 * ya * yb + 17 * ya * yc + 19 * yb * yc - 18 * SQUARE(yc))) / 80.0;
    matrix[2][4] = (27 * (-SQUARE(xb) + xa * (xb - xc) + xb * xc + (ya - yb) * (yb - yc))) / 16.0;
    matrix[2][5] = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80.0;
    matrix[2][6] = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80.0;
    matrix[2][7] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[2][8] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[2][9] = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40.0;

    matrix[3][3] = (17 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 40.0;
    matrix[3][4] = (-3 * (18 * SQUARE(xa) - 19 * xa * xb - 17 * xa * xc + 19 * xb * xc - SQUARE(xc) + 18 * SQUARE(ya) - 19 * ya * yb - 17 * ya * yc + 19 * yb * yc - SQUARE(yc))) / 80.0;
    matrix[3][5] = (3 * (9 * SQUARE(xa) + 8 * xb * xc + SQUARE(xc) - 2 * xa * (4 * xb + 5 * xc) + 9 * SQUARE(ya) - 8 * ya * yb - 10 * ya * yc + 8 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[3][6] = (-7 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80.0;
    matrix[3][7] = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80.0;
    matrix[3][8] = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80.0;
    matrix[3][9] = 0;

    matrix[4][4] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[4][5] = (-27 * (4 * SQUARE(xa) + SQUARE(xb) + 2 * xb * xc + SQUARE(xc) - 4 * xa * (xb + xc) + 4 * SQUARE(ya) - 4 * ya * yb + SQUARE(yb) - 4 * ya * yc + 2 * yb * yc + SQUARE(yc))) / 80.0;
    matrix[4][6] = (3 * (9 * SQUARE(xa) + SQUARE(xb) + 8 * xb * xc - 2 * xa * (5 * xb + 4 * xc) + (ya - yb) * (9 * ya - yb - 8 * yc))) / 80.0;
    matrix[4][7] = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[4][8] = (27 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 80.0;
    matrix[4][9] = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40.0;

    matrix[5][5] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[5][6] = (-3 * (xa - xb) * (18 * xa + xb - 19 * xc) - 3 * (ya - yb) * (18 * ya + yb - 19 * yc)) / 80.0;
    matrix[5][7] = (-27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 16.0;
    matrix[5][8] = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80.0;
    matrix[5][9] = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40.0;

    matrix[6][6] = (17 * (SQUARE(xa - xb) + SQUARE(ya - yb))) / 40.0;
    matrix[6][7] = (3 * (SQUARE(xa) - 18 * SQUARE(xb) + xa * (17 * xb - 19 * xc) + 19 * xb * xc + (ya - yb) * (ya + 18 * yb - 19 * yc))) / 80.0;
    matrix[6][8] = (3 * (SQUARE(xa) + 9 * SQUARE(xb) - 8 * xb * xc + xa * (-10 * xb + 8 * xc) + (ya - yb) * (ya - 9 * yb + 8 * yc))) / 80.0;
    matrix[6][9] = 0;

    matrix[7][7] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[7][8] = (-27 * (SQUARE(xa - 2 * xb + xc) + SQUARE(ya - 2 * yb + yc))) / 80.0;
    matrix[7][9] = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40.0;

    matrix[8][8] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16.0;
    matrix[8][9] = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40.0;

    matrix[9][9] = (81 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 20.0;

    // Fill the matrix with the given expressions
/*     matrix[0][0] = (17 * (SQUARE(xb - xc) + SQUARE(yb - yc))) / 40;
    matrix[0][1] = (3 * (SQUARE(xb) - 19 * xa * (xb - xc) + 17 * xb * xc - 18 * SQUARE(xc) - 19 * ya * yb + SQUARE(yb) + 19 * ya * yc + 17 * yb * yc - 18 * SQUARE(yc))) / 80;
    matrix[0][2] = (3 * (SQUARE(xb) + 8 * xa * (xb - xc) - 10 * xb * xc + 9 * SQUARE(xc) + 8 * ya * yb + SQUARE(yb) - 8 * ya * yc - 10 * yb * yc + 9 * SQUARE(yc))) / 80;
    matrix[0][3] = (-7 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80;
    matrix[0][4] = (-3 * (SQUARE(xb - xc) + SQUARE(yb - yc))) / 80;
    matrix[0][5] = (-3 * (SQUARE(xb) - 2 * xb * xc + SQUARE(xc) + SQUARE(yb - yc))) / 80;
    matrix[0][6] = (7 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80;
    matrix[0][7] = (3 * (-8 * xa * xb + 9 * SQUARE(xb) + 8 * xa * xc - 10 * xb * xc + SQUARE(xc) - 8 * ya * yb + 9 * SQUARE(yb) + 8 * ya * yc - 10 * yb * yc + SQUARE(yc))) / 80;
    matrix[0][8] = (3 * (-18 * SQUARE(xb) + 19 * xa * (xb - xc) + 17 * xb * xc + SQUARE(xc) + 19 * ya * yb - 18 * SQUARE(yb) - 19 * ya * yc + 17 * yb * yc + SQUARE(yc))) / 80;
    matrix[0][9] = 0;

    matrix[1][1] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16;
    matrix[1][2] = (-27 * (SQUARE(xa) + SQUARE(xb) + 2 * xa * (xb - 2 * xc) - 4 * xb * xc + 4 * SQUARE(xc) + SQUARE(ya) + 2 * ya * yb + SQUARE(yb) - 4 * ya * yc - 4 * yb * yc + 4 * SQUARE(yc))) / 80;
    matrix[1][3] = (3 * (SQUARE(xa) + 8 * xa * xb - 10 * xa * xc - 8 * xb * xc + 9 * SQUARE(xc) + SQUARE(ya) + 8 * ya * yb - 10 * ya * yc - 8 * yb * yc + 9 * SQUARE(yc))) / 80;
    matrix[1][4] = (-27 * (-SQUARE(xb) + xa * (xb - xc) + xb * xc + (ya - yb) * (yb - yc))) / 80;
    matrix[1][5] = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80;
    matrix[1][6] = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80;
    matrix[1][7] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80;
    matrix[1][8] = (-27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 16;
    matrix[1][9] = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40;

    matrix[2][2] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16;
    matrix[2][3] = (3 * (SQUARE(xa) - 19 * xa * xb + 17 * xa * xc + 19 * xb * xc - 18 * SQUARE(xc) + SQUARE(ya) - 19 * ya * yb + 17 * ya * yc + 19 * yb * yc - 18 * SQUARE(yc))) / 80;
    matrix[2][4] = (27 * (-SQUARE(xb) + xa * (xb - xc) + xb * xc + (ya - yb) * (yb - yc))) / 16;
    matrix[2][5] = (-27 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 80;
    matrix[2][6] = (-3 * (SQUARE(xa) - 2 * xa * xb + SQUARE(xb) + SQUARE(ya - yb))) / 80;
    matrix[2][7] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80;
    matrix[2][8] = (27 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80;
    matrix[2][9] = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40;

    matrix[3][3] = (17 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 40;
    matrix[3][4] = (-3 * (18 * SQUARE(xa) - 19 * xa * xb - 17 * xa * xc + 19 * xb * xc - SQUARE(xc) + 18 * SQUARE(ya) - 19 * ya * yb - 17 * ya * yc + 19 * yb * yc - SQUARE(yc))) / 80;
    matrix[3][5] = (3 * (9 * SQUARE(xa) + 8 * xb * xc + SQUARE(xc) - 2 * xa * (4 * xb + 5 * xc) + 9 * SQUARE(ya) - 8 * ya * yb - 10 * ya * yc + 8 * yb * yc + SQUARE(yc))) / 80;
    matrix[3][6] = (-7 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 80;
    matrix[3][7] = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80;
    matrix[3][8] = (-3 * (SQUARE(xa - xc) + SQUARE(ya - yc))) / 80;
    matrix[3][9] = 0;

    matrix[4][4] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16;
    matrix[4][5] = (-27 * (4 * SQUARE(xa) + SQUARE(xb) + 2 * xb * xc + SQUARE(xc) - 4 * xa * (xb + xc) + 4 * SQUARE(ya) - 4 * ya * yb + SQUARE(yb) - 4 * ya * yc + 2 * yb * yc + SQUARE(yc))) / 80;
    matrix[4][6] = (3 * (9 * SQUARE(xa) + SQUARE(xb) + 8 * xb * xc - 2 * xa * (5 * xb + 4 * xc) + (ya - yb) * (9 * ya - yb - 8 * yc))) / 80;
    matrix[4][7] = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80;
    matrix[4][8] = (27 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 80;
    matrix[4][9] = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40;

    matrix[5][5] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16;
    matrix[5][6] = (-3 * (xa - xb) * (18 * xa + xb - 19 * xc) - 3 * (ya - yb) * (18 * ya + yb - 19 * yc)) / 80;
    matrix[5][7] = (-27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 16;
    matrix[5][8] = (27 * ((xa - xc) * (xb - xc) + (ya - yc) * (yb - yc))) / 80;
    matrix[5][9] = (81 * ((xa - xb) * (xb - xc) + (ya - yb) * (yb - yc))) / 40;

    matrix[6][6] = (17 * (SQUARE(xa - xb) + SQUARE(ya - yb))) / 40;
    matrix[6][7] = (3 * (SQUARE(xa) - 18 * SQUARE(xb) + xa * (17 * xb - 19 * xc) + 19 * xb * xc + (ya - yb) * (ya + 18 * yb - 19 * yc))) / 80;
    matrix[6][8] = (3 * (SQUARE(xa) + 9 * SQUARE(xb) - 8 * xb * xc + xa * (-10 * xb + 8 * xc) + (ya - yb) * (ya - 9 * yb + 8 * yc))) / 80;
    matrix[6][9] = 0;

    matrix[7][7] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16;
    matrix[7][8] = (-27 * ((xa - 2 * xb + xc) * (xa - 2 * xb + xc) + (ya - 2 * yb + yc) * (ya - 2 * yb + yc))) / 80;
    matrix[7][9] = (-81 * (SQUARE(xa) + xb * xc - xa * (xb + xc) + (ya - yb) * (ya - yc))) / 40;

    matrix[8][8] = (27 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 16;
    matrix[8][9] = (-81 * (xa * (xb - xc) - xb * xc + SQUARE(xc) + ya * yb - ya * yc - yb * yc + SQUARE(yc))) / 40;

    matrix[9][9] = (81 * (SQUARE(xa) + SQUARE(xb) - xb * xc + SQUARE(xc) - xa * (xb + xc) + SQUARE(ya) - ya * yb + SQUARE(yb) - ya * yc - yb * yc + SQUARE(yc))) / 20;
 */
    // multiply by jacobian
    double jacobian = triangle.volume() * 2;
    for (size_t i = 0; i < 10; ++i)
    {
        for (size_t j = i; j < 10; ++j)
        {
            matrix[i][j] *= jacobian/compute_den(xa, ya, xb, yb, xc, yc);
            if (i != j)
            {
                matrix[j][i] = matrix[i][j];
            }
        }
    }

    //
}

double distance(const mesh::Point<2> &v1, const mesh::Point<2> &v2)
{
    return std::sqrt((v2.coords[0] - v1.coords[0]) * (v2.coords[0] - v1.coords[0]) + (v2.coords[1] - v1.coords[1]) * (v2.coords[1] - v1.coords[1]));
}

double dot(const mesh::Point<2> &v1, const mesh::Point<2> &v2, const mesh::Point<2> &w1, const mesh::Point<2> &w2)
{
    return (v2.coords[0] - v1.coords[0]) * (w2.coords[0] - w1.coords[0]) + (v2.coords[1] - v1.coords[1]) * (w2.coords[1] - w1.coords[1]);
}

#include <cstdlib>

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <N>" << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int N = std::atoi(argv[1]);

    mesh::SquareMaker mesh_maker(N);
    mesh::Mesh<2> mesh = mesh_maker.make_mesh();

    fe::FESimplexP3<2> fe;
    dof::DoFHandler<2> dof_handler(mesh);

    dof_handler.distribute_dofs(std::make_shared<fe::FESimplexP3<2>>(fe));

    unsigned int n_dofs = dof_handler.get_n_dofs();
    unsigned int n_dofs_per_cell = fe.get_n_dofs_per_element();

    auto f = [](double var1, double var2) { return 4 - 2 * (var1 * var1 + var2 * var2); };

    linalg::Vector rhs(n_dofs);

    linalg::CSRPattern csr_pattern = linalg::CSRPattern::create_from_dof_handler(dof_handler);

    linalg::CSRMatrix A(n_dofs, csr_pattern);

    linalg::FullMatrix local_matrix(n_dofs_per_cell);
    linalg::Vector local_rhs(n_dofs_per_cell);

    //dubious about f term, it has too many zeroes. if the rhs is manually filled at random, the solution obtained doesnt look that bad. maybe we need to interpolate f better. rhs has always norm very close to zero
    static double shape_integral_on_ref[10] = {1.0/60, 3.0/80, 3.0/80, 1.0/60, 3.0/80, 3.0/80, 1.0/60, 3.0/80, 3.0/80, 9.0/40};


    int token = 1;

    for (auto it = dof_handler.elem_begin(); it != dof_handler.elem_end(); ++it)
    {
        auto &elem = *it;

        mesh::Simplex<2, 2> triangle = mesh.get_Simplex(elem);


        if(token){
            token =1;
            std::cout << "Element vertices: " << elem.get_vertex(0) << " " << elem.get_vertex(1) << " " << elem.get_vertex(2);
            std::cout << "coords: " << "(" <<triangle.get_vertex(0).coords[0] << "," << triangle.get_vertex(0).coords[1] << ") ";
            std::cout << "(" <<triangle.get_vertex(1).coords[0] << "," << triangle.get_vertex(1).coords[1] << ") ";
            std::cout << "(" <<triangle.get_vertex(2).coords[0] << "," << triangle.get_vertex(2).coords[1] << ") " << std::endl;
            for(unsigned int i = 0; i < fe.get_n_dofs_per_element(); ++i)
            {
                mesh::Point<2> p = fe.get_dof_coords(triangle, i);
                std::cout << "position of dof " << i << " in element : " << p[0] << " " << p[1] << std::endl;
            }
        }


        mesh::Point<2> v0 = triangle.get_vertex(0);
        mesh::Point<2> v1 = triangle.get_vertex(1);
        mesh::Point<2> v2 = triangle.get_vertex(2);

        double volume = triangle.volume();

        local_matrix.set_to_zero();
        local_rhs.fill(0.0);

        double stiffness[10][10];
        compute_stiffness_loc(triangle, stiffness);

        // avarage of the function f over the element
        double avg = f(v0[0], v0[1]) + f(v1[0], v1[1]) + f(v2[0], v2[1]);
        avg /= 3.0;

        for(types::local_dof_index i = 0; i < n_dofs_per_cell; ++i)
        {
            for(types::local_dof_index j = 0; j < n_dofs_per_cell; ++j)
            {
                local_matrix(i, j) = stiffness[i][j];
            }

            local_rhs[i] = avg * shape_integral_on_ref[i] * 2 * volume;

            /**
             * TO CHECK
             */
            //local_rhs[i] += f(centroid.coords[0], centroid.coords[1]) * shape_integral_on_ref[i] * 2 * volume;
        }

        auto local_dofs = dof_handler.get_ordered_dofs_on_element(elem);
        std::cout << "dofs on element: " << std::endl;
        for (auto dof : local_dofs)
        {
            std::cout << dof << " ";
        }
        std::cout << std::endl;
        linalg::MatrixTools::add_local_matrix_to_global(A, local_matrix, local_dofs);
        linalg::MatrixTools::add_local_vector_to_global(rhs, local_rhs, local_dofs);

    } 


    linalg::MatrixTools::apply_homogeneous_dirichlet(A, rhs, dof_handler, 0);

    /* A.print();
    rhs.print(); */

    std::cout << "norm of rhs: " << rhs.norm() << std::endl;
    //output the rhs to see if it is correct

    fastfem::mesh::DataIO<2, 2> rhs_io(mesh, dof_handler, rhs);
    rhs_io.save_vtx("rhs.vtk");

    linalg::CGSolver solver(3000, 1e-12);
    linalg::Vector sol = solver.solve(A, rhs);

    fastfem::mesh::DataIO<2, 2> data_io(mesh, dof_handler, sol);
    data_io.save_vtx("solution_csr.vtk");

    std::cout << "Max of solution: " << sol.max() << std::endl;

    auto exact_f = [](double var1, double var2) { return (1 - var1 * var1) * (1 - var2 * var2); };

    linalg::Vector exact_sol(n_dofs);

    linalg::MatrixTools::interpolate(exact_sol, dof_handler, exact_f);

    std::cout << "Norm of difference: " << (sol - exact_sol).norm() << std::endl;

    mesh::DataIO<2, 2> data_io_exact(mesh, dof_handler, exact_sol);
    data_io_exact.save_vtx("exact_solution.vtk");

    return EXIT_SUCCESS;
}