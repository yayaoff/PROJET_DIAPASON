#include "../headers/opti.h"

// Méthode basique : méthode de bissection 
// Défaut : Optimisation univariée, or plusieurs paramètres !!
// On fixe 3 paramètres sur 4 et on effectue la méthode sur le dernier paramètre
// Amélioration potentielle : Utiliser une méthode d'optimisation multivariée (descente de gradient?) 

#define f_target 1244.51

double f(double param){
    // // Fonction dont on cherche la racine
    double f_actual;
    return f_actual - f_target;
}

double bissection_method(double low_bound, double up_bound){

    double mid = (low_bound + up_bound) / 2;
    while (fabs(low_bound-up_bound) > TOL) {
        if (f(mid) == 0) return mid;
        else if (f(mid)*f(low_bound) < 0) up_bound = mid;
        else low_bound = mid;
        mid = (up_bound + low_bound) / 2.0;
    }
    return mid;

}