#include <stdlib.h>
#include <math.h>
#include "../headers/matrix.h"
#include "../headers/eigen.h"

// Remove lines and columns corresponding to boundary nodes from K and M
void remove_bnd_lines (Matrix *K, Matrix *M, size_t *bnd_nodes, size_t n_bnd_nodes, Matrix **K_new, Matrix **M_new, int* invperm){
  size_t n = K->n;
  size_t n_new = K->n - 2*n_bnd_nodes;
  *K_new = allocate_matrix(n_new, n_new);
  *M_new = allocate_matrix(n_new, n_new);

  size_t *new_lines = malloc(n_new * sizeof(size_t));
  size_t i_bnd = 0;
  size_t i_new = 0;
  for (size_t i=0; i<(int)(n/2); i++){
    if (i == bnd_nodes[i_bnd]) {
      i_bnd++;
      continue;
    }
    new_lines[2*i_new  ] = 2*i;
    new_lines[2*i_new+1] = 2*i+1;
    i_new += 1;
  }
  for (size_t i=0; i<n_new; i++){
    for (size_t j=0; j<n_new; j++){
      (*K_new)->a[i][j] = K->a[new_lines[i]][new_lines[j]];
      (*M_new)->a[i][j] = M->a[new_lines[i]][new_lines[j]];
    }
  }
  free(new_lines);
}

// Matrix times vector product
int matrix_times_vector(Matrix *A, double* b, double* res){
  int n = A->n;
  for (int i=0; i<n; i++){
    double b_i = 0.;
    for (int j=0; j<n; j++){
      b_i += A->a[i][j]*b[j];
    }
    res[i] = b_i;
  }
  return 0;
}

// inner product of 2 vectors
double inner(double* a, double* b, int n){
  double res = 0.;
  for (int i=0; i<n; i++) res += a[i]*b[i];
  return res;
}

// normalize vector 
double normalize(double* a, int n){
  double norm = 0.;
  for (int i=0; i<n; i++) norm += a[i]*a[i];
  norm = sqrt(norm);
  for (int i=0; i<n; i++) a[i] /= norm;
  return norm;
}

double power_iteration(Matrix *A, double *v) {
  int n = A->n;

  // initial guess
  double Av[n];
  for (int i=0; i<n; i++) v[i] = 0;
  v[0] = v[1] = 1;
  normalize(v,n);

  double lambda = 1, lambda_prev, diff;
  double rtol = 1e-9; // relative tolerance on lambda

  for(int it = 0; it < 1e4; it++) {
    lambda_prev = lambda;
    // printf("\n===== Iteration %d =====\n", it+1);
    matrix_times_vector(A, v, Av);
    lambda = inner(v, Av, n);
    // printf("λ = %.9e\n", lambda);
    diff = fabs((lambda-lambda_prev) / lambda_prev);
    for(int i = 0; i < n; i++) v[i] = Av[i];
    normalize(v, n);
    if(diff < rtol) break;
  }

  return lambda;
}


double compute_freq(double r1, double r2, double e, double l, double meshSizeFactor, char * filename,int k,char *argv[],int argc){

  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);
  
  designTuningFork(6e-3, 11e-3, 38e-3, 82e-3, 0.3, NULL);
  
  // Number of vibration modes to find
  int k = atoi(argv[1]);

  // Assemble the 2 matrices of the linear elasticity problem: 
  // M is the mass matrix
  // K is the stiffness matrix
  Matrix *K, *M;
  size_t* boundary_nodes;
  size_t n_boundary_nodes;
  double * coord;
  assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);

  // Remove lines from matrix that are boundary
  Matrix *K_new;
  Matrix *M_new;
  remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);

  // Compute M := K^{-1} M 
  int n = K_new->n;
  lu(K_new);
  double Mj[n]; // column j
  for (int j=0; j<n; j++){ // column by column
    for (int i=0; i<n; i++) Mj[i] = M_new->a[i][j];
    solve(K_new, Mj);
    for (int i=0; i<n; i++){
      M_new->a[i][j] = Mj[i];
    }
  } // -> Now, in M_new we have stored K^-1M

  // Power iteration + deflation to find k largest eigenvalues
  Matrix * A = M_new; // alias for conciseness
  double * v = malloc(A->m * sizeof(double));
  double lambda, freq;
  FILE * file = fopen(argv[2], "w"); // open file to write frequencies
  for(int ki = 0; ki < k; ki++) {
    lambda = power_iteration(A, v);
    freq = 1./(2*M_PI*sqrt(lambda));

    fprintf(file, "%.9lf ", freq);

    printf("lambda = %.9e, f = %.3lf\n", lambda, freq);
    // Deflate matrix
    for(int i = 0; i < A->m; i++)
      for(int j = 0; j < A->n; j++)
        A->a[i][j] -= lambda * v[i] * v[j];

    // Put appropriate BC and plot
    double * vall = calloc(K->m, sizeof(double));
    int iv = 0, i_bnd = 0; 
    for(int i = 0; i < K->m/2; i++) {
      if(i_bnd < n_boundary_nodes && i == boundary_nodes[i_bnd]) {
        i_bnd++;
        continue;
      }
      vall[2*(i)]   = v[2*iv];
      vall[2*(i)+1] = v[2*iv+1];
      iv++;
    }
    visualize_in_gmsh(vall, K->m/2);
  }
  fclose(file);
  gmshFltkRun(&ierr);

  free_matrix (K);
  free_matrix (M);
  free_matrix (K_new);
  free_matrix (M_new);
  free(boundary_nodes);
}