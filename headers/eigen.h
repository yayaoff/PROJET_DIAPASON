#include <stdio.h>

void remove_bnd_lines (Matrix *K, Matrix *M, size_t *bnd_nodes, size_t n_bnd_nodes, Matrix **K_new, Matrix **M_new, int* invperm);

// inner product of 2 vectors
double inner(double* a, double* b, int n);

// normalize vector 
double normalize(double* a, int n);

double power_iteration(Matrix *A, double *v);

double compute_freq(double r1, double r2, double e, double l, double meshSizeFactor, char * filename,int k,char *argv[],int argc);