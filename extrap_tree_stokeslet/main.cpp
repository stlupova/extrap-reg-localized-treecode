#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <omp.h>

#include "stokeslets.h"

using namespace std;

void test_Stokes_SL_5ord(double h,
			 int N_quad, vector<Surf_point>& Surfc,
			 int N_target, vector<Target_point>& Target);

void set_rho_del(double h,
		 vector<double>& rho,
		 vector<double>& DEL);

void cout_norms(int n,
		const vector<double>& S_ex,
		const vector<double>& S_comp);

//*********************************************************************

int main(int argc, char** argv) {
  int N = N_gridlines;
  double h = 1.0 / N; // grid spacing

  double tm = omp_get_wtime();
  
  vector<Surf_point> Surfc; // array of surface (quadrature) points                    
  Generate_Surface(N, h, &Surfc);
  int N_quad = Surfc.size();
  cout << "Number of quadrature points = " << N_quad << endl;
  
  vector<Target_point> Target; // array of target points
  Generate_Targets(h, N_quad, Surfc, &Target);
  int N_target = Target.size();
  cout << "Number of target points = " << N_target << endl;
  
  test_Stokes_SL_5ord(h, N_quad, Surfc, N_target, Target);
  
  tm = omp_get_wtime() - tm;
  cout << "CPU time = " << tm << " seconds" << endl;
  
  return 0;
}

//*********************************************************************
//*********************************************************************

void test_Stokes_SL_5ord(double h,
			 int N_quad, vector<Surf_point>& Surfc,
			 int N_target, vector<Target_point>& Target) {

  cout << "Computing Stokeslet integral by extrapolated reg'n with 3 del's..." << endl;
  
  int N_quad3 = 3 * N_quad;
  int N_target3 = 3 * N_target;
  
  for (int i=0; i<N_quad; i++) {
    Surfc[i].f = stokeslet_density(Surfc[i].x,Surfc[i].Nrml);
  }

  vector<double> SL_ex(N_target3,0);
  Exact_solution(N_target, Target, SL_ex);

  vector<double> rho(3,0), DEL(3,0);
  set_rho_del(h, rho, DEL);

  double tm = omp_get_wtime();
  
  vector<double> SL_comp(N_target3,0);
  if (use_tree == 1) { // use treecode to speed up calculations
    eval_Stokes_SL_offSurf_3del_tree(h, rho,
				     DEL[0], DEL[1], DEL[2],
				     N_quad, Surfc,
				     N_target, Target,
				     SL_comp); 
  }
  else {
    eval_Stokes_SL_offSurf_3del(h, rho,
				DEL[0], DEL[1], DEL[2],
				N_quad, Surfc,
				N_target, Target,
				SL_comp); 
  }  
  
  tm = omp_get_wtime() - tm;
  cout << "CPU time to compute the Stokeslet integral = " << tm << " seconds" << endl;

  // display norms of solution and error
  cout_norms(N_target, SL_ex, SL_comp);
}

//*********************************************************************

void set_rho_del(double h,
		 vector<double>& rho,
		 vector<double>& DEL) {
  
  if (del_h == 1) {
    rho[0] = 3.0; rho[1] = 4.0; rho[2] = 5.0;
    for (int i=0; i<3; i++) {
      DEL[i] = rho[i] * h;
    }
  }
  else {
    rho[0] = 2.0; rho[1] = 3.0; rho[2] = 4.0;
    for (int i=0; i<3; i++) {
      DEL[i] = rho[i] * pow(1.0/64.0,1.0/5.0) * pow(h,4.0/5.0);
    }
  }
}

//*********************************************************************

void cout_norms(int n,
		const vector<double>& S_ex,
		const vector<double>& S_comp) {
		
  double S_max, S_l2;
  int ind_S_max;  
  norms_3d_vector(n, S_ex, S_max, ind_S_max, S_l2);
  cout << "max exact value = " << S_max << "   l2 exact values = " << S_l2 << endl;
  
  double err_max, err_l2;
  int ind_max;  
  error_3d_vector(n, S_comp, S_ex, err_max, ind_max, err_l2);
  cout << "max error = " << err_max << "   l2 error = " << err_l2 << endl;

}
