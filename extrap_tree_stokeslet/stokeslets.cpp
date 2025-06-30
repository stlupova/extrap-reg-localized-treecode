#include <iostream>
#include <vector>
#include <iomanip>
#include <omp.h>

#include "stokeslets.h" 

using namespace std;

//*********************************************************************
// evaluates stokes single layer potential (Surfc.f is density)

void eval_Stokes_SL_offSurf_3del(double h, const vector<double>& rho,
				 double DEL1, double DEL2, double DEL3,
				 int N_quad, const vector<Surf_point>& Surfc,
				 int N_target, const vector<Target_point>& Target,
				 vector<double>& SL) {

  double EightPI = 1.0/(8.0*PI);

  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions

#pragma omp parallel 
{
  vector<double> pt(3, 0), dx(3, 0);
  double f0_dot_n0;

  cout << "Thread number " << omp_get_thread_num() << endl;

#pragma omp for
  for (int i=0; i<N_target; i++) {    

    pt[0] = Target[i].x;
    pt[1] = Target[i].y;
    pt[2] = Target[i].z;

    f0_dot_n0 = dot_product(Target[i].nrst_f, Target[i].nrst_Nrml);
    
    vector<double> u1(3,0), u2(3,0), u3(3,0);
    for (int j=0; j<N_quad; j++) {


      for (int k=0; k<3; k++)  dx[k] = pt[k] - Surfc[j].x[k];      

      double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      
      vector<double> f = Surfc[j].f;

      for (int k=0; k<3; k++)  f[k] -= f0_dot_n0 * Surfc[j].Nrml[k];
      
      double f_dot_dx = dot_product(f, dx);

      double H1_d1, H2_d1, H1_d2, H2_d2, H1_d3, H2_d3;
      H_gauss(r, DEL1, H1_d1, H2_d1); 
      H_gauss(r, DEL2, H1_d2, H2_d2); 
      H_gauss(r, DEL3, H1_d3, H2_d3); 

      for (int k=0; k<3; k++) {
	double m1 = f[k] * Surfc[j].Area;
	double m2 = dx[k] * f_dot_dx * Surfc[j].Area;
	
	u1[k] += m1 * H1_d1 + m2 * H2_d1;
	u2[k] += m1 * H1_d2 + m2 * H2_d2;
	u3[k] += m1 * H1_d3 + m2 * H2_d3;
      }
    }

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
}
