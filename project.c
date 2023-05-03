#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../gmsh-sdk/include/gmshc.h"
#include "headers/matrix.h"
#include "headers/elasticity.h"
#include "math.h"
#include "headers/lu.h"
#include "headers/design.h"
#include "headers/eigen.h"
#include "headers/opti.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

int k;
double r1; double r2; double e; double l; double meshSizeFactor; char * filename;

// tag_vis == 1 si on souhaite visualiser le résultat , 0 sinon
double *compute_freq( char *argv[],int argc, int tag_vis){
  
  // Define physical constants
  double E = 0.7e11;  // Young's modulus for Aluminum
  double nu = 0.3;    // Poisson coefficient
  double rho = 3000;  // Density of Aluminum

  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);
  
  designTuningFork(r1, r2, l, e, meshSizeFactor, filename);
  
  // Number of vibration modes to find
  k = atoi(argv[1]);
  double *frequencies = malloc(sizeof(double)*k);

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

    frequencies[ki] = freq;

    if(tag_vis==1) fprintf(file, "%.9lf ", freq);

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
    if(tag_vis ==1) visualize_in_gmsh(vall, K->m/2);
  }
  fclose(file);
  if(tag_vis ==1) gmshFltkRun(&ierr);

  free_matrix (K);
  free_matrix (M);
  free_matrix (K_new);
  free_matrix (M_new);
  free(boundary_nodes);
  return frequencies;
}

int main (int argc, char *argv[]) {

  if (argc < 2){
    printf("Usage: \n"
			"./project <k> <out>\n" 
			"---------------------------- \n\n"
			"- k is the number of frequencies to compute. \n "
			"- out is the output file to write the frequencies. \n "
      "\n");
		return -1;
  } 

  double f_target =  1244.51 ; //frequence for D5# [Hz] 
  // We want fabs(freq-f_target) <= EPS

  // Parameters to optimize : r1, r2, e, l
  // Ici, on optimise le paramètre r1 et on fixe les autres arbitrairement-> Tester pour les 4 et voir
  // pour quel paramètre on converge le plus vite ?
  // Rmq : attention aux contraintes : r1 < r2
  r1=6e-3; r2=11e-3; e=38e-3; l=82e-3; meshSizeFactor=0.3;filename=NULL;
  double* frequencies = compute_freq(argv,argc,0);
  double lower_bound = 0; double upper_bound = 1;
  for(int mode_k=0; mode_k < k; k++){
    double freq_act = frequencies[mode_k];
    printf("Frequence before optimisation (k=%d) = %f\n",mode_k,freq_act);
    double param_biss = bissection_method(lower_bound,upper_bound);
    printf("Optimisation ended : found new paramater %f for previous paramater %f\n",param_biss,r1);
    // Tester pour le nouveau paramètre
    r1 = param_biss;
  }
  double* new_frequencies = compute_freq(argv,argc,1);
  // for(int mode_k=0; mode_k < k; k++) printf("Frequence for mode %d = %f\n",mode_k, new_frequencies[mode_k]);
  return 0;
}