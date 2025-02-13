#include "fem/navier_stokes.hpp"

namespace fem::navier_stokes
{   

    class Navier_stokes_solver
    {   
        Navier_stokes_solver(const Mesh& mesh, double dt, double T):
        init_omega(), build_fem_matrices(), mesh(mesh), dt(dt), T(T), N_DOFs(mesh.vtx_count), omega(N_DOFs), Momega(N_DOFs), psi(N_DOFs), M_col_sum(N_DOFs)
        {
            // at the beginning we calculate the volume of the mesh and the column
            volume = 0;
            for(int i=0; i<M.nnz; i++){
                volume += M.coeffs[i].val;
            }

            // we compute M_col_sum once and for all at the beginning of the code
            for(int i=0; i<M.nnz; i++){
                M_col_sum[M.coeffs[i].i] += M.coeffs[i].val;
            }

            // set the initial phi to zero
            std::fill(myVector.begin(), myVector.end(), 0);
        };
         
        void init_omega()
        {
            for(int i = 0; i < N_DOFs; i++)
            {
                const auto &v = mesh.verteces[i];
                // TODO: implement omega_initializer
                omega[i] = omega_initializer(v.x, v.y, v.z); // set the initial value of omega for each vertex
            }
        };

        // set the mean omega to zero
        void set_zero_mean()
        {
            linalg::matrix_vector_product(&M, omega.data(), Momega.data()); 
            double sum = 0;
            // we compute the 
            for(int i=0; i<N_DOFs; i++){
                sum += Momega[i];
            }
            double alpha = sum/volume;

            for(int i=0; i<N_DOFs; i++){
                omega[i] -= alpha;
            }

            linalg::blas_axpby(-alpha, M_col_sum.data(), +1, Momega.data(), N_DOFs);
        }

        void build_fem_matrices()
        {
            S->rows = S->cols = M->rows = M->cols = N_DOFs;

            // Simplified version of the calculation, we are using more memory
            // than necessary but it is easier to understand
            S->nnz = M->nnz = 9 * mesh.tri_count;
            S->coeffs = (struct Coeff *)malloc(S->nnz * sizeof(struct Coeff));
            M->coeffs = (struct Coeff *)malloc(M->nnz * sizeof(struct Coeff));

            // For each triangle
            for(int t=0; t< mesh.tri_count; t++){
                // get the vertices of the current triangle
                int a = mesh.triangles[t].a;
                int b = mesh.triangles[t].b;
                int c = mesh.triangles[t].c;
                // get the coordinates of the vertices
                struct Vertex A = mesh.vertices[a];
                struct Vertex B = mesh.vertices[b];
                struct Vertex C = mesh.vertices[c];

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
        };

        void compute_stream_function(){
            // set the mean of omega to zero
            set_zero_mean();

            // solve the system S phi ^ t+dt = -M omega ^ t+dt

            // conjugate gradient
            
            /* Set phi to the initial guess = 0 */
            //8std::fill(myVector.begin(), myVector.end(), 0); // maybe we can avoid this, the initial guess is the last solution

            Momega = linalg::matrix_vector_product(M, omega.data(), Momega.data()); // Momega = M * omega

            // r now is Momega

            double *p = array(N);
            memcpy(p, r, N * sizeof(double));


            //double *Ap = array(N);
            //double *Mp = array(N);
            
            //double error2 = blas_dot(r, r, N);

            double tol2 = 1e-6;
            int max_iter = 1000;
            int iterate = 0;

            while (error2 > tol2 && iterate < max_iter)
            {
                matrix_vector_product(S, p, Ap);
                matrix_vector_product(M, p, Mp);
                blas_axpby(1, Mp, 1, Ap, N);

                double alpha = blas_dot(r, r, N) / blas_dot(p, Ap, N); // alpha_k = (r_k, r_k) / (p_k, A * p_k)

                blas_axpby(alpha, p, 1, U, N); // u_k+1 = u_k + alpha_k * p_k

                blas_axpby(-alpha, Ap, 1, r, N); // r_k+1 = r_k - alpha_k * A * p_k

                // if(blas_dot(r, r, N) < tol2) break;

                double beta = blas_dot(r, r, N) / error2; // beta_k = (r_k+1, r_k+1) / (r_k, r_k)

                blas_axpby(1, r, beta, p, N); // p_k+1 = r_k+1 + beta_k * p_k

                error2 = blas_dot(r, r, N);

                iterate++;
            }             
        }

        void step()
        {
            // 
            // compute stream function (set_zero_mean)
            // assemble transport
            // solve backward euler

        }


        void solve()
        {
            for(double t=0; t<T; t+=dt){
                step();
            }
        }




    }

} // fem::navier_stokes