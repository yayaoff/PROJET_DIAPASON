#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "eigen.h"

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
