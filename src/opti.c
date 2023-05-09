#include "../headers/opti.h"

// Méthode basique : méthode de bissection 
// Défaut : Optimisation univariée, or plusieurs paramètres !!
// On fixe 3 paramètres sur 4 et on effectue la méthode sur le dernier paramètre
// Amélioration potentielle : Utiliser une méthode d'optimisation multivariée (descente de gradient?) 

// tag_vis == 1 si on souhaite visualiser le résultat , 0 sinon
double *compute_freq(int tag_vis, param_t *param){
  
  // Define physical constants
  double E = 0.7e11;  // Young's modulus for Aluminum
  double nu = 0.3;    // Poisson coefficient
  double rho = 3000;  // Density of Aluminum
  // Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(0,NULL, 0, 0, &ierr);

  designTuningFork(param->r1, param->r2, param->l, param->e, param->meshSizeFactor, NULL);
  
  // Number of vibration modes to find
  int k = param->k;
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
  // FILE * file = fopen(argv[2], "w"); // open file to write frequencies
  for(int ki = 0; ki < k; ki++) {
    lambda = power_iteration(A, v);
    freq = 1./(2*M_PI*sqrt(lambda));

    frequencies[ki] = freq;

    // if(tag_vis==1) fprintf(file, "%.9lf ", freq);

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
    // if(tag_vis ==1) visualize_in_gmsh(vall, K->m/2);
  }
  gmshFltkRun(&ierr);
  free_matrix (K);
  free_matrix (M);
  free_matrix (K_new);
  free_matrix (M_new);
  free(boundary_nodes);
  return frequencies;
}

double f(double param_to_opt,param_t* param_act){
    // // Fonction dont on cherche la racine
    param_act->r1 = param_to_opt;
    printf("Testing for param : %lf\n",param_to_opt);
    double f_actual = compute_freq(0, param_act)[0];
    return f_actual - f_target;
}

double bissection_method(double low_bound, double up_bound,param_t* param){
    
    double mid = (low_bound + up_bound) / 2;
    while (fabs(low_bound-up_bound) > TOL) {
        if (f(mid,param) == 0) return mid;
        else if (f(mid,param)*f(low_bound,param) < 0) up_bound = mid;
        else low_bound = mid;
        mid = (up_bound + low_bound) / 2.0;
    }
    return mid;
    
}