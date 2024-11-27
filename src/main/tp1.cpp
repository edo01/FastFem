/**
 * @TODO:
 * - One possible optimization that can be carried is the optimization of the deduplication. 
 *   The current implementation search duplicates verteces on the entire mesh, when
 *   duplicates can be found only on the borders. On possible solution is to limit the search only on 
 *   the borders of the square. One way could be to mark each vertex to signal borders and then search
 *   only on the marked verteces.
 *
 * - Implement compensated dot product and see how the accuracy of the solver is affected
 * 
 * - 
 */

#include "mesh/mesh.hpp"
#include "mesh/mesh_tools.hpp"
#include "mesh/cube.hpp"
#include "fem/poisson.hpp"
#include "linalg/system.hpp"
#include "linalg/sparse_matrix.hpp"
#include "omp.h"
using namespace mesh;
using namespace linalg;
using namespace fem::poisson;

/******************************************************************************
 * Let's choose our right hand side f of -\Delta u + u = f
 *****************************************************************************/
double f(double x, double y, double z)
{
	(void)z; /* Avoids compiler warning about unused variable */
	return x * x - y * y;
}

/******************************************************************************
 * Main routine
 *****************************************************************************/
int main(int argc, char **argv)
{
	if (argc < 2){
		printf("Usage: %s N\n", argv[0]);
		return EXIT_FAILURE;
	}

	//print omp num threads
	printf("Number of threads: %d\n", omp_get_max_threads());

	struct Mesh m;

	// Build a cube mesh with N subdivisions per side
	build_cube_mesh(&m, atoi(argv[1]));
	send_cube_to_sphere(m.vertices, m.vtx_count);
	int N = m.vtx_count;
	printf("Number of DOF : %d\n", N);

	struct SparseMatrix M;
	struct SparseMatrix S;
	build_fem_matrices(&m, &S, &M);

	/* Fill F */
	double *F = array(N);
	for (int i = 0; i < N; i++) {
		struct Vertex v = m.vertices[i];
		F[i] = f(v.x, v.y, v.z);
	}
	/* Fill B = MF */
	double *B = array(N);
	matrix_vector_product(&M, F, B);

	/* Solve (S + M)U = B */
	double *U = array(N);
	// int iter = gradient_system_solve(&S, &M, B, U, N);
 	int iter = CG_system_solve(&S, &M, B, U, N);
	printf("System solved in %d iterations.\n", iter);

	printf("Integrity check :\n");
	printf("-----------------\n");
	for (int i = 0; i < 8; i++) {
		if (abs(F[i]) > 1e-12) {
			printf("Ratio U/F : %f\n", U[i] / F[i]);
		}
	}

	return (EXIT_SUCCESS);
}





