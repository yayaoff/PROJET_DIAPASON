#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lu.h"

int lu(Matrix * A) {
	int m = A->m;
	Matrix * LU = A; // factorization is done in-place
	for(int k = 0; k < m-1; k++) {
		if(fabs(LU->a[k][k]) < EPS) return -1;
		for(int j = k+1; j < m; j++) {
			LU->a[j][k] /= LU->a[k][k];
			for(int l = k+1; l < m; l++)
				LU->a[j][l] -= LU->a[j][k] * LU->a[k][l];
		}
	}
	return 0;
}

int solve(Matrix * LU, double * y) {
	int m = LU->m;
	double * x = y;

	// Résolution de L*x = y par substitution avant
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < i; j++)
			x[i] -= LU->a[i][j] * x[j];
	}
	// Résolution de U*x = L^{-1} y par substitution arrière
	for(int i = m-1; i >= 0; i--) {
		for(int j = i+1; j < m; j++)
			x[i] -= LU->a[i][j] * x[j];
		x[i] /= LU->a[i][i];
	}
	
	return 0;
}