#include "elasticity.h"

static void p1_stifness_matrix (double dphi[6], double h[3][3] , double det, double S[6][6]){

  double d1dx = dphi[0];
  double d1dy = dphi[1];
  double d2dx = dphi[2];
  double d2dy = dphi[3];
  double d3dx = dphi[4];
  double d3dy = dphi[5];
  double BT[3][6] = {{d1dx, d2dx, d3dx, 0, 0, 0},
		     {0, 0, 0, d1dy, d2dy, d3dy},
		     {d1dy, d2dy, d3dy, d1dx, d2dx, d3dx}};
  double Bh[6][3];
  for (int i=0;i < 6; i++){
    for (int j=0;j < 3; j++){
      Bh[i][j] = 0;
      for (int k=0;k < 3; k++){
	Bh[i][j] += BT[k][i]*h[k][j];
      }
    }
  }

  for (int i=0;i < 6; i++){
    for (int j=0;j < 6; j++){
      S[i][j] = 0;
      for (int k=0;k < 3; k++){
	S[i][j] += det*0.5*Bh[i][k]*BT[k][j];
      }
    }
  }
}

static void hooke_plain_stress (double E, double nu, double h[3][3]){
  double C = E/(1-nu*nu);
  h[0][0] = h[1][1] = C;
  h[2][2] = C*(1-nu);
  h[0][1] = h[1][0] = C*nu;
  h[0][2] = h[2][0] = h[1][2] = h[2][1] = 0;
}

void p1_stifness_matrix_plane_stress (double E, double nu, double dphi[6], double det, double S[6][6]){
  double HookeMatrix[3][3];
  hooke_plain_stress (E,nu,HookeMatrix);    
  p1_stifness_matrix (dphi,HookeMatrix,det,S);  
}

void p1_geometry(const double *x, double *detptr, double dxidx[2][2], double *dphi)
{
  double dxdxi[2][2] = {
    {x[1*2+0]-x[0*2+0],x[2*2+0]-x[0*2+0]},
    {x[1*2+1]-x[0*2+1],x[2*2+1]-x[0*2+1]}
  };
  double det = dxdxi[0][0]*dxdxi[1][1]-dxdxi[0][1]*dxdxi[1][0];
  dxidx[0][0] = dxdxi[1][1]/det;
  dxidx[0][1] = -dxdxi[0][1]/det;
  dxidx[1][0] = -dxdxi[1][0]/det;
  dxidx[1][1] = dxdxi[0][0]/det;
  for (int i = 0; i < 2; ++i) {
    dphi[0*2+i] = -dxidx[0][i] - dxidx[1][i];
    dphi[1*2+i] = dxidx[0][i];
    dphi[2*2+i] = dxidx[1][i];
  }
  *detptr = det;
}


void p1_laplace_matrix (double alpha, double dphi[6], double det, double S[3][3]){

  double d1dx = dphi[0];
  double d1dy = dphi[1];
  double d2dx = dphi[2];
  double d2dy = dphi[3];
  double d3dx = dphi[4];
  double d3dy = dphi[5];
  double BT[3][2] = {{d1dx, d1dy},
		     {d2dx, d2dy},
		     {d3dx, d3dy}};
  for (int i=0;i < 3; i++){
    for (int j=0;j < 3; j++){
      S[i][j] = 0;
      for (int k=0;k < 2; k++){
	S[i][j] += alpha*det*0.5*BT[i][k]*BT[k][j];
      }
    }
  }
}

int get_triangles_p1 (size_t ** elementTags, size_t * elementTags_n,
		       size_t ** nodeTags, size_t * nodeTags_n) {
  int ierr;
  gmshModelMeshGetElementsByType(2, elementTags, elementTags_n,
				 nodeTags, nodeTags_n, -1, 0, 1, &ierr);
  return ierr;
}

int get_nodes (size_t ** nodeTags, size_t * nodeTags_n,
		double ** coord, size_t * coord_n) {
  int ierr;
  gmshModelMeshGetNodes(nodeTags, nodeTags_n, coord, coord_n, 0, 0, 2, -1, 1, 0, &ierr);
  return ierr;
}

Matrix * compute_stiffnes (int nNodes,
				       double *coord,
				       size_t ntriangles,
				       size_t *triangleNodes, 
               double E, 
               double nu){
  Matrix *K = allocate_matrix (2*nNodes, 2*nNodes);
  for (size_t i = 0 ; i < ntriangles ; i++){
    size_t *n = & triangleNodes [3*i];
    double x[6] = {coord[3*(n[0]-1)],coord[3*(n[0]-1)+1],
		   coord[3*(n[1]-1)],coord[3*(n[1]-1)+1],
		   coord[3*(n[2]-1)],coord[3*(n[2]-1)+1]};
    
    double det,dxidx[2][2],dphi[6];
    p1_geometry(x,&det, dxidx, dphi);
    double StiffnessMatrix[6][6];
    p1_stifness_matrix_plane_stress (E,nu,dphi,det,StiffnessMatrix);
    
    size_t dofs[6] = {2*(n[0]-1),2*(n[1]-1),2*(n[2]-1),
		      2*(n[0]-1)+1,2*(n[1]-1)+1,2*(n[2]-1)+1};
    for (size_t j = 0;j<6;j++){
      for (size_t k = 0;k<6;k++){
        K->a[dofs[j]][dofs[k]] += StiffnessMatrix[j][k];
      }
    }
  }
  return K;
}

Matrix * compute_mass (int nNodes,
				   double *coord,
				   size_t ntriangles,
				   size_t *triangleNodes, double rho){
  Matrix *M =allocate_matrix(2*nNodes, 2*nNodes);
  for (size_t i = 0 ; i < ntriangles ; i++){
    size_t *n = & triangleNodes [3*i];
    double x[6] = {coord[3*(n[0]-1)],coord[3*(n[0]-1)+1],
		   coord[3*(n[1]-1)],coord[3*(n[1]-1)+1],
		   coord[3*(n[2]-1)],coord[3*(n[2]-1)+1]};
    
    double det,dxidx[2][2],dphi[6];
    p1_geometry(x,&det, dxidx, dphi);
    const double A = rho * det / 12;
    const double B = A/2.0;
    double MassMatrix[6][6] ={
      {A,B,B,0,0,0},
      {B,A,B,0,0,0},
      {B,B,A,0,0,0},
      {0,0,0,A,B,B},
      {0,0,0,B,A,B},
      {0,0,0,B,B,A}
    };
    size_t dofs[6] = {2*(n[0]-1),2*(n[1]-1),2*(n[2]-1),
		      2*(n[0]-1)+1,2*(n[1]-1)+1,2*(n[2]-1)+1};
    for (size_t j = 0;j<6;j++){
      for (size_t k = 0;k<6;k++){
        M->a[dofs[j]][dofs[k]] += MassMatrix[j][k];
      }
    }
  }
  return M;
}

void get_boundary_nodes (size_t** bnd_nTags, size_t* bnd_nTags_n){
  // boundary conditions -- physical names allowed are --> for generalised eigen value problem, only "0" boundary conditions are allowed
  // This function returns node numbers that should be removed from matrix
  int *dimTags = NULL, ierr;
  size_t dimTags_n;
  gmshModelGetPhysicalGroups(&dimTags, &dimTags_n, -1, &ierr);
  for (size_t i=0;i<dimTags_n;i+=2){  
    char *name =malloc(525*sizeof(char));
    gmshModelGetPhysicalName(dimTags[i],dimTags[i+1],&name, &ierr);
    // "clamped" (force x and y = 0)
    double *cc = NULL;
    size_t  cc_n;
    if (strcmp(name,"clamped") == 0){
      gmshModelMeshGetNodesForPhysicalGroup(dimTags[i],dimTags[i+1], bnd_nTags, bnd_nTags_n, &cc, &cc_n, &ierr);
      printf("clamped found physical group %d\n",dimTags[i+1]);
    }
    if (cc) free(cc);
    free (name);
  }
  for (size_t i=0; i<*bnd_nTags_n; i++) (*bnd_nTags)[i]--;
  if(dimTags)free(dimTags);
}

int assemble_system(Matrix** K, Matrix** M, double** coord, size_t** boundary_nodes, size_t* n_boundary_nodes, double E, double nu, double rho){
  int ierr;
  
  gmshModelMeshRebuildNodeCache(1,&ierr);

  size_t *triangleTags=0, ntriangles, *triangleNodes=NULL, temp;
  get_triangles_p1 (&triangleTags, &ntriangles,&triangleNodes, &temp);

  size_t *nodeTags=0, nNodes;
  double *coord_bad = 0;
  get_nodes (&nodeTags, &nNodes, &coord_bad, &temp);

  get_boundary_nodes(boundary_nodes, n_boundary_nodes); // the position in the matrix is at bnd_nodes[i]-1

  *coord = (double*)malloc(3*nNodes*sizeof(double));
  // Gmsh does not provide nodes in the right order.
  for (size_t i=0;i<nNodes;i++){
    (*coord)[3*(nodeTags[i]-1)] = coord_bad[3*i];
    (*coord)[3*(nodeTags[i]-1)+1] = coord_bad[3*i+1];
    (*coord)[3*(nodeTags[i]-1)+2] = coord_bad[3*i+2];
  }
  
  *K = compute_stiffnes (nNodes,*coord,ntriangles,triangleNodes, E, nu);
  *M = compute_mass (nNodes,*coord,ntriangles,triangleNodes, rho);
  
  free (coord_bad);
  free (nodeTags);
  free (triangleTags);
  free (triangleNodes);
  
  return 0;
}

void visualize_in_gmsh(double* SOL, int n_nodes){
  int ierr;

  size_t *nodeTags=malloc(n_nodes*sizeof(size_t));
  // double *coord_bad = 0;
  // get_nodes (&nodeTags, &nNodes, &coord_bad, &temp);


  double* sol_3D = malloc(3*n_nodes*sizeof(double)); // in gmsh, vectors are 3D
  for (int i=0; i<n_nodes; i++){
    nodeTags[i] = i+1;
    sol_3D[3*i+0] = SOL[2*i+0];
    sol_3D[3*i+1] = SOL[2*i+1];
    sol_3D[3*i+2] = 0.;
  }

  // for (int i=0; i<3*n_nodes; i++) printf("%f \n", sol_3D[i]);
  // printf("\n");
  char *model_name;
  gmshModelGetCurrent(&model_name, &ierr);
  int view_tag = gmshViewAdd("displacement", -1, &ierr);
  gmshViewAddHomogeneousModelData(view_tag, 0, model_name,"NodeData", nodeTags, n_nodes, sol_3D, 3*n_nodes, 0., 3, -1, &ierr);
  gmshViewOptionSetNumber(view_tag, "VectorType", 5, &ierr);
  gmshViewOptionSetNumber(view_tag, "DisplacementFactor", 0.1, &ierr);
}