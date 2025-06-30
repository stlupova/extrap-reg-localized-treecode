
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <time.h>
#include <omp.h>

#include "KITC.h"

using namespace std;


//**********************************************************//
vector<panel> tree;
vector<size_t> leaf;
int max_level = 0;


static size_t node_count = 0;

//**********************************************************************//
long getTickCount() {
  return clock();
}

//**********************************************************************//
void build_tree_init(int N_cube) {
  panel temp_panel;
  
  // indices of particles belonging to panel
  temp_panel.members[0] = 0;
  temp_panel.members[1] = N_cube - 1;
  
  temp_panel.xinterval[0] = -1.0; // interval defining the panel
  temp_panel.xinterval[1] = 1.0;
  temp_panel.yinterval[0] = -1.0;
  temp_panel.yinterval[1] = 1.0;
  temp_panel.zinterval[0] = -1.0;
  temp_panel.zinterval[1] = 1.0;
  temp_panel.xc = 0.5 * (temp_panel.xinterval[0] + temp_panel.xinterval[1]);
  temp_panel.yc = 0.5 * (temp_panel.yinterval[0] + temp_panel.yinterval[1]);
  temp_panel.zc = 0.5 * (temp_panel.zinterval[0] + temp_panel.zinterval[1]);
  temp_panel.MAC = (3 * 10 * 10/ 4) / sq_theta; // MAC = r^2 / theta^2
  
  tree.push_back(temp_panel);
  node_count = 1;
}

//**********************************************************************//
void Swap(size_t i, size_t j, double *lambda[3], struct xyz &s) {
  if (i == j)
    return;
  
  double x = s.x[i];
  double y = s.y[i];
  double z = s.z[i];
  size_t index = s.index[i];
  size_t old_index = s.old_index[i];
  double nx = s.nx[i];
  double ny = s.ny[i];
  double nz = s.nz[i];
  double area = s.area[i];
  double lam0 = lambda[0][i];
  double lam1 = lambda[1][i];
  double lam2 = lambda[2][i];
  
  s.x[i] = s.x[j];
  s.y[i] = s.y[j];
  s.z[i] = s.z[j];
  s.index[i] = s.index[j];
  s.old_index[i] = s.old_index[j];
  s.nx[i] = s.nx[j];
  s.ny[i] = s.ny[j];
  s.nz[i] = s.nz[j];
  s.area[i] = s.area[j];
  lambda[0][i] = lambda[0][j];
  lambda[1][i] = lambda[1][j];
  lambda[2][i] = lambda[2][j];
  
  s.x[j] = x;
  s.y[j] = y;
  s.z[j] = z;
  s.index[j] = index;
  s.old_index[j] = old_index;
  s.nx[j] = nx;
  s.ny[j] = ny;
  s.nz[j] = nz;
  s.area[j] = area;
  lambda[0][j] = lam0;
  lambda[1][j] = lam1;
  lambda[2][j] = lam2;
}
//***********************************************************************//
void split_tree_node(size_t panel_index,
		     double *lambda[3],
		     struct xyz &particles) {

  panel child[8];
  
  double tp_x0 = tree[panel_index].xinterval[0];
  double tp_x1 = tree[panel_index].xinterval[1];
  double tp_y0 = tree[panel_index].yinterval[0];
  double tp_y1 = tree[panel_index].yinterval[1];
  double tp_z0 = tree[panel_index].zinterval[0];
  double tp_z1 = tree[panel_index].zinterval[1];
  
  double midpointx = (tp_x0 + tp_x1) / 2.0;
  double midpointy = (tp_y0 + tp_y1) / 2.0;
  double midpointz = (tp_z0 + tp_z1) / 2.0;
  
  double xc0 = (tp_x0 + midpointx) / 2.0;
  double xc1 = (tp_x1 + midpointx) / 2.0;
  double yc0 = (tp_y0 + midpointy) / 2.0;
  double yc1 = (tp_y1 + midpointy) / 2.0;
  double zc0 = (tp_z0 + midpointz) / 2.0;
  double zc1 = (tp_z1 + midpointz) / 2.0;
  
  child[0].xinterval[0] = tp_x0;
  child[0].xinterval[1] = midpointx;
  child[0].yinterval[0] = tp_y0;
  child[0].yinterval[1] = midpointy;
  child[0].zinterval[0] = tp_z0;
  child[0].zinterval[1] = midpointz;
  child[0].xc = xc0;
  child[0].yc = yc0;
  child[0].zc = zc0;
  child[0].MAC = ((midpointx - xc0) * (midpointx - xc0)
		  + (midpointy - yc0) * (midpointy - yc0)
		  + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
  
  child[1].xinterval[0] = midpointx;
  child[1].xinterval[1] = tp_x1;
  child[1].yinterval[0] = tp_y0;
  child[1].yinterval[1] = midpointy;
  child[1].zinterval[0] = tp_z0;
  child[1].zinterval[1] = midpointz;
  child[1].xc = xc1;
  child[1].yc = yc0;
  child[1].zc = zc0;
  child[1].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1)
		  + (midpointy - yc0) * (midpointy - yc0)
		  + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
  
  child[2].xinterval[0] = tp_x0;
  child[2].xinterval[1] = midpointx;
  child[2].yinterval[0] = midpointy;
  child[2].yinterval[1] = tp_y1;
  child[2].zinterval[0] = tp_z0;
  child[2].zinterval[1] = midpointz;
  child[2].xc = xc0;
  child[2].yc = yc1;
  child[2].zc = zc0;
  child[2].MAC = ((midpointx - xc0) * (midpointx - xc0)
		  + (tp_y1 - yc1) * (tp_y1 - yc1)
		  + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
  
  child[3].xinterval[0] = midpointx;
  child[3].xinterval[1] = tp_x1;
  child[3].yinterval[0] = midpointy;
  child[3].yinterval[1] = tp_y1;
  child[3].zinterval[0] = tp_z0;
  child[3].zinterval[1] = midpointz;
  child[3].xc = xc1;
  child[3].yc = yc1;
  child[3].zc = zc0;
  child[3].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1)
		  + (tp_y1 - yc1) * (tp_y1 - yc1)
		  + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
  
  child[4].xinterval[0] = tp_x0;
  child[4].xinterval[1] = midpointx;
  child[4].yinterval[0] = tp_y0;
  child[4].yinterval[1] = midpointy;
  child[4].zinterval[0] = midpointz;
  child[4].zinterval[1] = tp_z1;
  child[4].xc = xc0;
  child[4].yc = yc0;
  child[4].zc = zc1;
  child[4].MAC = ((midpointx - xc0) * (midpointx - xc0)
		  + (midpointy - yc0) * (midpointy - yc0)
		  + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
  
  child[5].xinterval[0] = midpointx;
  child[5].xinterval[1] = tp_x1;
  child[5].yinterval[0] = tp_y0;
  child[5].yinterval[1] = midpointy;
  child[5].zinterval[0] = midpointz;
  child[5].zinterval[1] = tp_z1;
  child[5].xc = xc1;
  child[5].yc = yc0;
  child[5].zc = zc1;
  child[5].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1)
		  + (midpointy - yc0) * (midpointy - yc0)
		  + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
  
  child[6].xinterval[0] = tp_x0;
  child[6].xinterval[1] = midpointx;
  child[6].yinterval[0] = midpointy;
  child[6].yinterval[1] = tp_y1;
  child[6].zinterval[0] = midpointz;
  child[6].zinterval[1] = tp_z1;
  child[6].xc = xc0;
  child[6].yc = yc1;
  child[6].zc = zc1;
  child[6].MAC = ((midpointx - xc0) * (midpointx - xc0)
		  + (tp_y1 - yc1) * (tp_y1 - yc1)
		  + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
  
  child[7].xinterval[0] = midpointx;
  child[7].xinterval[1] = tp_x1;
  child[7].yinterval[0] = midpointy;
  child[7].yinterval[1] = tp_y1;
  child[7].zinterval[0] = midpointz;
  child[7].zinterval[1] = tp_z1;
  child[7].xc = xc1;
  child[7].yc = yc1;
  child[7].zc = zc1;
  child[7].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1)
		  + (tp_y1 - yc1) * (tp_y1 - yc1)
		  + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
  
  vector<size_t> v[8];
  size_t start = tree[panel_index].members[0];
  size_t end = tree[panel_index].members[1];
  size_t* addr_table = new size_t[end - start + 1];
  
  size_t index;
  for (index = start; index <= end; index++)
    {
      particles.index[index] = index;
      addr_table[index - start] = index;
      
      if (particles.x[index] <= midpointx
	  && particles.y[index] <= midpointy &&
	  particles.z[index] <= midpointz)
	v[0].push_back(index);
      else if (particles.x[index] > midpointx
	       && particles.y[index] <= midpointy &&
	       particles.z[index] <= midpointz )
	v[1].push_back(index);
      else if (particles.x[index] <= midpointx
	       && particles.y[index] > midpointy &&
	       particles.z[index]<= midpointz)
	v[2].push_back(index);
      else if (particles.x[index] > midpointx
	       && particles.y[index] > midpointy &&
	       particles.z[index] <= midpointz)
	v[3].push_back(index);
      else if(particles.x[index] <= midpointx
	      && particles.y[index] <= midpointy &&
	      particles.z[index] > midpointz )
	v[4].push_back(index);
      else if (particles.x[index] > midpointx
	       && particles.y[index] <= midpointy &&
	       particles.z[index] > midpointz)
	v[5].push_back(index);
      else if (particles.x[index] <= midpointx
	       && particles.y[index] > midpointy &&
	       particles.z[index] > midpointz)
	v[6].push_back(index);
      else if (particles.x[index] > midpointx
	       && particles.y[index] > midpointy &&
	       particles.z[index] > midpointz)
	v[7].push_back(index);
    }
  
  size_t seq = start;
  for (size_t j = 0; j < 8; j++)
    {
      size_t size = v[j].size();
      
      if (size >= 1)
	{
	  for (size_t k = 0; k < size; k++)
	    {
	      if (k == 0)
		child[j].members[0] = seq;
	      if (k == size - 1)
		child[j].members[1] = seq;
	      
	      index = v[j][k];
	      // This uses an address table
	      size_t pos = addr_table[index - start];
	      size_t out = particles.index[seq];
	      Swap(pos, seq, lambda, particles);
	      addr_table[index - start] = seq;
	      addr_table[out - start] = pos;
	      
	      seq++;
	    }
	  
	  node_count++;
	  tree[panel_index].children.push_back(node_count - 1);
	  tree.push_back(child[j]);
	  v[j].clear();
	}
    }
  
  delete[] addr_table;
}

//***********************************************************************//
void build_tree_3D_Recursive(size_t panel_index,
			     double *lambda[3],
			     struct xyz &particles,
			     int level) {
  if (level > max_level)
    max_level = level;
  
  size_t n = tree[panel_index].members[1]
    - tree[panel_index].members[0] + 1;
  
  if (n >= (size_t)N0) {
    split_tree_node(panel_index, lambda, particles);
    
    for (size_t i = 0; i < tree[panel_index].children.size(); i++) {
      size_t panel_index_new = tree[panel_index].children[i];
      build_tree_3D_Recursive(panel_index_new,
			      lambda,
			      particles,
			      level + 1);
    }
  }
  else
    leaf.push_back(panel_index);
}

//***********************************************************************//
void Panel_Moments(size_t panel_index,
		   double *lambda[3],
		   struct xyz &particles,
		   double m[][Pflat]) {
  double t1[P + 1];
  double t2[P + 1];
  double t3[P + 1];
  
  for (int i = 0; i < P + 1; i++) {
    t1[i] = tree[panel_index].t1[i];
    t2[i] = tree[panel_index].t2[i];
    t3[i] = tree[panel_index].t3[i];
  }
  
  double w1i[P + 1];
  double w2j[P + 1];
  double w3k[P + 1];
  double dj[P + 1];
  dj[0] = 0.5;
  dj[P] = 0.5;
  for (int j = 1; j<P; j++)
    dj[j] = 1;
  
  for (int j = 0; j < P + 1; j++)
    w3k[j] = w2j[j] = w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
  
  int tp0 = tree[panel_index].members[0];
  int tp1 = tree[panel_index].members[1];
  int cluster_size;
  cluster_size = tp1 - tp0 + 1;
  
  double **a1i;
  double **a2j;
  double **a3k;
  double *D;
  
  a1i = (double**)calloc(P + 1, sizeof(double*));
  a2j = (double**)calloc(P + 1, sizeof(double*));
  a3k = (double**)calloc(P + 1, sizeof(double*));
  D = (double*)calloc(cluster_size + 1, sizeof(double));
  for (int i = 0; i < (P + 1); i++) {
    a1i[i] = (double*)calloc(cluster_size + 1, sizeof(double));
    a2j[i] = (double*)calloc(cluster_size + 1, sizeof(double));
    a3k[i] = (double*)calloc(cluster_size + 1, sizeof(double));
  }
  

  for (size_t tp_j = tp0; tp_j <= tp1 ; tp_j++) {
    
    double x = particles.x[tp_j];
    double y = particles.y[tp_j];
    double z = particles.z[tp_j];

    int flag1 = -1;
    int flag2 = -1;
    int flag3 = -1;
    double SumA1 = 0.0;
    double SumA2 = 0.0;
    double SumA3 = 0.0;
    
    for (int j = 0; j < P + 1; j++) {
      if ( abs(x-t1[j])<1.e-07) { //removable singularity
	flag1 = j;
      }
      else {
	double temp1 = w1i[j] / (x - t1[j]);
	a1i[j][tp_j - tp0] = temp1;
	SumA1 += temp1;
      }
	
      if ( abs(y-t2[j])<1.e-07) { //removable singularity
	flag2 = j;
      }
      else {
	double temp2 = w2j[j] / (y - t2[j]);
	a2j[j][tp_j - tp0] = temp2;
	SumA2 += temp2;
      }
      
      if ( abs(z-t3[j])<1.e-07) { //removable singlarity
	flag3 = j;
      }
      else {
	double temp3 = w3k[j] / (z - t3[j]);
	a3k[j][tp_j - tp0] = temp3;
	SumA3 += temp3;
      }
      
    }

    // if a flag was set, adjust sum and a to handle removable singularity
    if (flag1 > -1) {
      SumA1 = 1.0;
      for (int j = 0; j < P + 1; j++) a1i[j][tp_j - tp0] = 0.0;
      a1i[flag1][tp_j - tp0] = 1.0;
    }
    if (flag2 > -1) {
      SumA2 = 1.0;
      for (int j = 0; j < P + 1; j++) a2j[j][tp_j - tp0] = 0.0;
      a2j[flag2][tp_j - tp0] = 1.0;
    }
    if (flag3 > -1) {
      SumA3 = 1.0;
      for (int j = 0; j < P + 1; j++) a3k[j][tp_j - tp0] = 0.0;
      a3k[flag3][tp_j - tp0] = 1.0;
    }
    
    D[tp_j - tp0] = 1.0 / (SumA1 * SumA2 * SumA3);
  }
  
  int kk = -1;
  for (int i = 0; i < P + 1; i++) {
    for (int j = 0; j < P + 1; j++) {
      for (int k = 0; k < P + 1; k++) {	
	kk = kk + 1;
	
	double sum[6] = {0};
	for (size_t tp_j = tp0; tp_j <= tp1; tp_j++) {
	  
	  double nx = particles.nx[tp_j];
	  double ny = particles.ny[tp_j];
	  double nz = particles.nz[tp_j];
	  
	  double area = particles.area[tp_j];
    
	  double s = a1i[i][tp_j - tp0] * a2j[j][tp_j - tp0]
	    * a3k[k][tp_j - tp0] * D[tp_j - tp0];
	  
	  s *= area;
	  
	  sum[0] += s * lambda[0][tp_j];
	  sum[1] += s * lambda[1][tp_j];
	  sum[2] += s * lambda[2][tp_j];
	  
	  sum[3] += s * nx;
	  sum[4] += s * ny;
	  sum[5] += s * nz;
	}
	for (int jj=0; jj<6; jj++) m[jj][kk] = sum[jj];
      }
    }
  }
  
  for (int i = 0; i < P + 1; i++) {
    free(a1i[i]);
    free(a2j[i]);
    free(a3k[i]);
  }
  
  free(a1i);
  free(a2j);
  free(a3k);
  free(D); 
}

//**********************************************************************//

vector<double> Call_Treecode_noDel(double x, double y, double z,
				   double g0_dot_n0,
				   int panel_index) {
  
  vector<double> del_x(3,0), temp_mom(3,0), u(3,0);
  
  double dx[P + 1];
  double dy[P + 1];
  double dz[P + 1];

  for (int i = 0; i < P + 1; i++) {
    dx[i] = x - tree[panel_index].t1[i];
    dy[i] = y - tree[panel_index].t2[i];
    dz[i] = z - tree[panel_index].t3[i];
  }
  
  int kk = -1;
  for (int i = 0; i < P + 1; i++) {
    double temp_i = dx[i] * dx[i];
    for (int j = 0; j < P + 1; j++) {
      double temp_j = dy[j] * dy[j];
      for (int k = 0; k < P + 1; k++) {
	double temp_k = dz[k] * dz[k];
	kk = kk + 1;
	
	del_x[0] = dx[i];
	del_x[1] = dy[j];
	del_x[2] = dz[k];
	
	temp_mom[0] = tree[panel_index].moments[0][kk]
	  - g0_dot_n0 * tree[panel_index].moments[3][kk];
	temp_mom[1] = tree[panel_index].moments[1][kk]
	  - g0_dot_n0 * tree[panel_index].moments[4][kk];
	temp_mom[2] = tree[panel_index].moments[2][kk]
	  - g0_dot_n0 * tree[panel_index].moments[5][kk];

	double r2 = temp_i + temp_j + temp_k;
	double r = sqrt(r2);
	
	double temp = dot_product(temp_mom, del_x) / (r2 * r);
	
	u[0] += temp_mom[0] / r + temp * dx[i];
	u[1] += temp_mom[1] / r + temp * dy[j];
	u[2] += temp_mom[2] / r + temp * dz[k];
      }
    }
  }
  return u;
}

//**********************************************************************//
void Cluster_Chev_Points(size_t tree_size) {
  double h;
  h = 3.14159265358979323846/P;
  double t[P + 1] = {0.0};
  for (int i = 0; i < P + 1; i++)
    t[i] = cos(i * h);  //Chebyshev interpolation points [-1,1]
  
  double x1,x2,y1,y2,z1,z2;
  size_t tree_index;
  
  for (tree_index = 0; tree_index < tree_size ; tree_index++) {
    x1 = tree[tree_index].xinterval[0];
    x2 = tree[tree_index].xinterval[1];
    y1 = tree[tree_index].yinterval[0];
    y2 = tree[tree_index].yinterval[1];
    z1 = tree[tree_index].zinterval[0];
    z2 = tree[tree_index].zinterval[1];
    
    for (int i = 0; i < P + 1; i++) {// map to the cluster
      tree[tree_index].t1[i] =  x1 + (t[i] + 1)/2 * (x2 - x1);
      tree[tree_index].t2[i] =  y1 + (t[i] + 1)/2 * (y2 - y1);
      tree[tree_index].t3[i] =  z1 + (t[i] + 1)/2 * (z2 - z1);
    }
  }
}

//***********************************************************************//
vector<double> Call_Ds(int limit_1, int limit_2,
		       double p_x, double p_y, double p_z,
		       double g0_dot_n0,
		       struct xyz &particles,
		       double *lambda[3],
		       double DEL) {    

  vector<double> dx(3,0), nrml(3,0), g(3,0), u(3,0);
  
  for (size_t j = limit_1; j <= limit_2; j++) {
    
    dx[0] = p_x - particles.x[j];
    dx[1] = p_y - particles.y[j];
    dx[2] = p_z - particles.z[j];
    
    nrml[0] = particles.nx[j];
    nrml[1] = particles.ny[j];
    nrml[2] = particles.nz[j];

    double area = particles.area[j];
    
    g[0] = lambda[0][j] - g0_dot_n0 * nrml[0]; //use subtraction in Stokeslet
    g[1] = lambda[1][j] - g0_dot_n0 * nrml[1];
    g[2] = lambda[2][j] - g0_dot_n0 * nrml[2];

    double r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

    double H1, H2;
    H_gauss(r, DEL, H1, H2); 
    H2 *= dot_product(g, dx);
    
    u[0] += (g[0] * H1  +  dx[0] * H2) * area;
    u[1] += (g[1] * H1  +  dx[1] * H2) * area;
    u[2] += (g[2] * H1  +  dx[2] * H2) * area;      
  }
  return u;
}

//***********************************************************************//
void Call_Ds_3del(int limit_1, int limit_2,
		  double p_x, double p_y, double p_z,
		  double g0_dot_n0,
		  struct xyz &particles,
		  double *lambda[3],
		  double DEL1, double DEL2, double DEL3,
		  vector<double>& u1,
		  vector<double>& u2,
		  vector<double>& u3) {    
  
  vector<double> dx(3,0), nrml(3,0), g(3,0), u(3,0);
  
  for (size_t j = limit_1; j <= limit_2; j++) {
    
    dx[0] = p_x - particles.x[j];
    dx[1] = p_y - particles.y[j];
    dx[2] = p_z - particles.z[j];
    
    nrml[0] = particles.nx[j];
    nrml[1] = particles.ny[j];
    nrml[2] = particles.nz[j];

    double area = particles.area[j];
    
    g[0] = (lambda[0][j] - g0_dot_n0 * nrml[0]) * area; //use subtraction in Stokeslet
    g[1] = (lambda[1][j] - g0_dot_n0 * nrml[1]) * area;
    g[2] = (lambda[2][j] - g0_dot_n0 * nrml[2]) * area;

    double r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    
    double g_dot_dx = dot_product(g, dx);

    double H1, H2;
    H_gauss(r, DEL1, H1, H2); 
    H2 *= g_dot_dx;

    u1[0] += g[0] * H1  +  dx[0] * H2;
    u1[1] += g[1] * H1  +  dx[1] * H2;
    u1[2] += g[2] * H1  +  dx[2] * H2;      

    H_gauss(r, DEL2, H1, H2); 
    H2 *= g_dot_dx;
        
    u2[0] += g[0] * H1  +  dx[0] * H2;
    u2[1] += g[1] * H1  +  dx[1] * H2;
    u2[2] += g[2] * H1  +  dx[2] * H2;      

    H_gauss(r, DEL3, H1, H2); 
    H2 *= g_dot_dx;
        
    u3[0] += g[0] * H1  +  dx[0] * H2;
    u3[1] += g[1] * H1  +  dx[1] * H2;
    u3[2] += g[2] * H1  +  dx[2] * H2;      
  }
}

//**********************************************************************//

void eval_Stokes_Sum_cluster_3del(double *lambda[3],
				  struct xyz &particles,
				  double p_x, double p_y, double p_z,
				  double f0_dot_n0,
				  size_t panel_index,
				  double DEL1, double DEL2, double DEL3,
				  vector<double>& Sum1,
				  vector<double>& Sum2,
				  vector<double>& Sum3) {	
  
  
  size_t limit_1 = tree[panel_index].members[0];
  size_t limit_2 = tree[panel_index].members[1];  
  
  double xc = tree[panel_index].xc;
  double yc = tree[panel_index].yc;
  double zc = tree[panel_index].zc;
    
  double tpx = p_x - xc;
  double tpy = p_y - yc;
  double tpz = p_z - zc;
  
  double R_sq = tpx * tpx + tpy * tpy + tpz * tpz;

  
  if (tree[panel_index].MAC < R_sq) { //particle and cluster well-separated

    if (use_tree_approximate == 1) {
      vector<double> tree_result = Call_Treecode_noDel(p_x, p_y, p_z,
						       f0_dot_n0,
						       panel_index);
      for (int k=0; k<3; k++) { //use the sum with del1 for all three delta's
	Sum1[k] += tree_result[k];
	Sum2[k] += tree_result[k];
	Sum3[k] += tree_result[k];
      }
    }
    else {
      vector<double> DS_result = Call_Ds(limit_1, limit_2,
					 p_x, p_y, p_z,
					 f0_dot_n0, 
					 particles,
					 lambda,
					 DEL1);
      for (int k=0; k<3; k++) { //use the sum with del1 for all three delta's
	Sum1[k] += DS_result[k];
	Sum2[k] += DS_result[k];
	Sum3[k] += DS_result[k];
      }
    }
  }
  else {
    if ((limit_2 - limit_1) < N0) { //if cluster is a leaf, use direct sum
      
      // this is the "nearby" case, and the del values have to be done separately

      vector<double> DS_result_del1(3,0), DS_result_del2(3,0), DS_result_del3(3,0);
      Call_Ds_3del(limit_1, limit_2,
		   p_x, p_y, p_z,
		   f0_dot_n0,
		   particles,
		   lambda,
		   DEL1, DEL2, DEL3,
		   DS_result_del1,
		   DS_result_del2,
		   DS_result_del3);

      for (int k=0; k<3; k++)  {
	Sum1[k] += DS_result_del1[k];
	Sum2[k] += DS_result_del2[k];
	Sum3[k] += DS_result_del3[k];
      }
    }
    
    else {//if cluster is not a leaf, look at children recursively
     
      for (int k=0; k<3; k++) {
	Sum1[k] = 0.0;
	Sum2[k] = 0.0;
	Sum3[k] = 0.0;
      }
      
      size_t length = tree[panel_index].children.size();
      for (size_t i = 0; i < length; i++) {
	
	size_t index = tree[panel_index].children[i];
	vector<double> temp_res1(3,0), temp_res2(3,0), temp_res3(3,0);
	eval_Stokes_Sum_cluster_3del(lambda,
				     particles,
				     p_x, p_y, p_z,
				     f0_dot_n0,
				     index,
				     DEL1, DEL2, DEL3,
				     temp_res1,
				     temp_res2,
				     temp_res3);
	for (int k=0; k<3; k++) {
	  Sum1[k] += temp_res1[k];
	  Sum2[k] += temp_res2[k];
	  Sum3[k] += temp_res3[k];
	}
      }
    }
  }
}

//****************************************************************************//
void eval_Stokes_SL_offSurf_3del_tree(double h,  const vector<double>& rho,
				      double DEL1, double DEL2, double DEL3,
				      int N_quad, const vector<Surf_point>& Surfc,
				      int N_target, const vector<Target_point>& Target,
				      vector<double>& SL) {

  cout << "Treecode with theta = " << sqrt(sq_theta) << "  N0 = " << N0
       << "  p = " << P << endl;
  
  double theta = sqrt(sq_theta);
  double EightPI = 1.0/(8.0*PI);

  int N_cube = N_quad;
  
  tree.reserve(5000);
  leaf.reserve(5000);
  struct xyz particles(N_cube);
  double *lambda[3];
  lambda[0] = new double[N_cube];
  lambda[1] = new double[N_cube];
  lambda[2] = new double[N_cube];

  //create a physical copy of the input data
  //so the treecode can modify the data
  for (int i=0; i<N_cube; i++) {

    particles.x[i] = Surfc[i].x[0];
    particles.y[i] = Surfc[i].x[1];
    particles.z[i] = Surfc[i].x[2];
    particles.index[i] = -1;
    particles.old_index[i] = i;

    particles.nx[i] = Surfc[i].Nrml[0];
    particles.ny[i] = Surfc[i].Nrml[1];
    particles.nz[i] = Surfc[i].Nrml[2];

    particles.area[i] = Surfc[i].Area;

    lambda[0][i] = Surfc[i].f[0];
    lambda[1][i] = Surfc[i].f[1];
    lambda[2][i] = Surfc[i].f[2];
  }
    
  //***************** Set up tree *******************************
  clock_t Start_total, Start_btree;
  clock_t End_total, End_btree;
  
  Start_total = clock(); // Get currenct CPU time
  Start_btree = clock();
  
  build_tree_init(N_cube);
  build_tree_3D_Recursive(0, lambda, particles, 0);
 
  End_btree = clock();
  
  // //***************** Compute moments for each panel ************
  size_t size = tree.size();
  Cluster_Chev_Points(size);


  clock_t moments_tm = clock();
  for (size_t i = 1; i < size; i++) {// skip root
    Panel_Moments(i, lambda, particles, tree[i].moments);
  }
  moments_tm = clock() - moments_tm;
  cout << "     CPU for moments " << ((float)moments_tm)/CLOCKS_PER_SEC
       << " seconds" << endl;

  
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions

#pragma omp parallel 
{
  cout << "Thread number " << omp_get_thread_num() << endl;

#pragma omp for
  for (int i=0; i<N_target; i++) {    

    double p_x = Target[i].x;
    double p_y = Target[i].y;
    double p_z = Target[i].z;

    double f0_dot_n0 = dot_product(Target[i].nrst_f, Target[i].nrst_Nrml);

    vector<double> u1(3,0), u2(3,0), u3(3,0);
    eval_Stokes_Sum_cluster_3del(lambda,
				 particles,
				 p_x, p_y, p_z,
				 f0_dot_n0,
				 0,
				 DEL1, DEL2, DEL3,
				 u1, u2, u3);

    vector<double> u_comp(3,0);
    solve_3by3(Target[i].b, rho,
	       DEL1, DEL2, DEL3,
	       u1, u2, u3,
	       u_comp);
    
    for (int k=0; k<3; k++) {
      SL[3*i+k] = u_comp[k] * EightPI;
    }
  }
 }
  End_total = clock(); 
  
  cout << "     build tree time = " << ((float)(End_btree - Start_btree))/CLOCKS_PER_SEC
       << " seconds" << endl;  
  cout << "     total treecode time = " << ((float)(End_total - Start_total))/CLOCKS_PER_SEC
       << " seconds" << endl;
  
  delete [] lambda[0];
  delete [] lambda[1];
  delete [] lambda[2];
  tree.clear();
  leaf.clear();
}


