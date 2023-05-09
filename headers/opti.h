#include "eigen.h"
#include <math.h>
#include <stdlib.h>
#include "lu.h"
#include "elasticity.h"
#include "../../gmsh-sdk/include/gmshc.h"
#include "design.h"

#define TOL 1e-9
#define f_target 1244.51

typedef struct{
    int k;
    double r1;
    double r2;
    double e;
    double l; 
    double meshSizeFactor;
}param_t;

double *compute_freq(int tag_vis, param_t *param);
double f(double param_to_opt, param_t* param_act);
double bissection_method(double low_bound, double up_bound, param_t *param);