#include "fem/poisson.hpp"

namespace fem::poisson
{
    using namespace mesh;
    using namespace linalg;

/**
 * Given a simplician conformal mesh defined over S^2, with N vertices,
 * the stiffness matrix S and the mass matrix M are defined as follows:
 * - S=(s_ij) where s_ij = \int_{\Omega} \nabla \phi_j \cdot \nabla \phi_i + \phi_j \phi_i
 * 						 = a(\phi_j, \phi_i) with i,j=0..N-1
 * - M=(m_ij) where m_ij = \int_{\Omega} \phi_i \phi_j with i,j=0..N-1
 * 
 * We can decompose the integral into multiple integrals defined on the elements
 * of the triangulation. For example the mass matrix:
 * m_ij = \sum_k^{N_elements}\int_{\T_k} \phi_i \phi_j,
 * where T_k are the elements of the triangulation.
 * 
 * Computing a(\phi_j, \phi_i) for every i and j is not necessary and it is 
 * extremely inefficient since we are dealing with sparse systems. In fact,
 * a(\phi_j, \phi_i) = 0 when the two supports of \phi_i \phi_j does not 
 * intersect.
 * 
 * A smarter solution leverages this property and computes only non-zero entries.
 * In this way we can compute, for each element of the mesh, 3*3 integrals: one 
 * for each ordered pair of element's vertices.  
 * 
 * An important clarification must be made: even if we perform 9*#elements integrals,
 * most of them will be aggregated in the same entry of the matrix. This is due to the
 * fact that two adjacent triangles share a two common vertices. 
 * 
 * Another important consideration is that both matrices are symmetric, so we can
 * compute only the upper triangular part of the matrix and then copy the values
 * in the lower triangular part, but for simplicity we will compute all the entries.
 *  
 */
void build_fem_matrices(const struct Mesh *m, struct SparseMatrix *S,
			struct SparseMatrix *M)
{

    int N = m->vtx_count; // There is a DOF per vertex
	S->rows = S->cols = M->rows = M->cols = m->vtx_count;

	// Simplified version of the calculation, we are using more memory
	// than necessary but it is easier to understand
	S->nnz = M->nnz = 9 * m->tri_count;
	S->coeffs = (struct Coeff *)malloc(S->nnz * sizeof(struct Coeff));
	M->coeffs = (struct Coeff *)malloc(M->nnz * sizeof(struct Coeff));

	// For each triangle
	for(int t=0; t< m->tri_count; t++){
		// get the vertices of the current triangle
		int a = m->triangles[t].a;
		int b = m->triangles[t].b;
		int c = m->triangles[t].c;
		// get the coordinates of the vertices
		struct Vertex A = m->vertices[a];
		struct Vertex B = m->vertices[b];
		struct Vertex C = m->vertices[c];

		struct Vector AB = vector(A, B);
		struct Vector BC = vector(B, C);
		struct Vector CA = vector(C, A);
		
		// compute the area of the triangle
		double area = norm(cross(AB, CA))* 0.5;

		// get the pointer to the next 9 coefficients
		struct Coeff *mass = &M->coeffs[9*t];
		// diagonal terms
		mass[0] = {a, a, area/6};
		mass[1] = {b, b, area/6};
		mass[2] = {c, c, area/6};
		// off-diagonal terms
		mass[3] = {a, b, area/12};
		mass[4] = {b, a, area/12};
		mass[5] = {a, c, area/12};
		mass[6] = {c, a, area/12};
		mass[7] = {b, c, area/12};
		mass[8] = {c, b, area/12};

		// get the pointer to the next 9 coefficients
		struct Coeff *stiffness = &S->coeffs[9*t];
		double r = 1. / (4 * area);
		// diagonal terms
		stiffness[0] = {a, a, dot(BC, BC) * r};
		stiffness[1] = {b, b, dot(CA, CA) * r};
		stiffness[2] = {c, c, dot(AB, AB) * r};

		// off-diagonal terms
		stiffness[3] = {a, b, dot(BC, CA) * r};
		stiffness[4] = {b, a, dot(BC, CA) * r};
		stiffness[5] = {a, c, dot(AB, BC) * r};
		stiffness[6] = {c, a, dot(AB, BC) * r};
		stiffness[7] = {b, c, dot(AB, CA) * r};
		stiffness[8] = {c, b, dot(AB, CA) * r};
	}
}

} // fem::poisson