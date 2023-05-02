#ifndef MATRIX_H // le header guard
#define MATRIX_H 

typedef struct Triplet {
	int i,j;	// index of unknowns
    double val;	// value of entry to add to matrix
} Triplet;

typedef struct Matrix {
	int m, n;       // dimensions de la matrice
	double * data;  // tableau 1D de taille m*n contenant les entrées de la matrice
	double ** a;    // tableau 1D de m pointeurs vers chaque ligne, pour pouvoir appeler a[i][j]
} Matrix;

typedef struct BandMatrix {
	int m, k;		// dimension de la matrice et largeur de bande
	double * data;	// tableau 1D de taille m*(2*k+1) contenant les entrées de la matrice
	double ** a;	// tableau 1D de m pointeurs vers chaque ligne, pour pouvoir appeler a[i][j]
} BandMatrix;


Matrix * allocate_matrix(int m, int n); // alloue une matrice de dimensions données
BandMatrix * allocate_band_matrix(int m, int k); // alloue une matrice bande de dimensions données
void free_matrix(Matrix * mat);         // libère la mémoire de la structure Matrix donnée
void free_band_matrix(BandMatrix * mat);         // libère la mémoire de la structure BandMatrix donnée

// Calcule une permutation pour réduire le fill-in.
// On peut utiliser les coordonnées des noeuds ainsi que les entrées de la matrice (éventuellement)
int compute_permutation(int * perm, double * coord, int n_nodes); 

void print_vector(double * v, int n);   // imprime le contenu d'un vecteur (tableau) de taille n
void print_matrix(Matrix * A);          // imprime le contenu d'une matrice

int is_symmetric(Matrix * K);

#endif