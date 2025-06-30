#include <iostream>
#include <vector>

#include "utilities.h"

using namespace std;

//*********************************************************************
// Level set function representing the surface

double phi(const vector<double>& x) {
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;
  
  return x[0]*x[0]/a2[0] + x[1]*x[1]/a2[1] + x[2]*x[2]/a2[2] - 1.0; 
}

//*********************************************************************
// D_phi/D_x(i): i-th derivative of phi

double Dphi(int i, const vector<double>& x) {
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;
  
  return 2.0*x[i]/a2[i]; 
}

//*********************************************************************
// Second derivatives of phi

void D2phi(const vector<double>& x, double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23, double& phi31,
	   double& phi32, double& phi33) {
  
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;

  phi11 = 2.0/a2[0]; phi12 = 0.0;       phi13 = 0.0;
  phi21 = 0.0;       phi22 = 2.0/a2[1]; phi23 = 0.0;
  phi31 = 0.0;       phi32 = 0.0;       phi33 = 2.0/a2[2];
}

//*********************************************************************
//*********************************************************************

// Regularization
void H_gauss(double r, double d, double& H1, double& H2) {
  if (r < 1e-14) {
    H1 = 2.0/rootPI/d;
    H2 = 4.0/3.0/rootPI/(d*d*d);
  }
  else if ( r/d > 6.0) {
    H1 = 1.0/r;
    H2 = 1.0/(r*r*r);
  }
  else {
    double s1, s2;
    s_opt2(r/d,s1,s2);
    H1 = s1/r;
    H2 = s2/(r*r*r);
  }
}

//*********************************************************************

void s_opt2(double r, double& s1, double& s2) {
  s1 = erf(r);
  s2 = s1 - 2.0*r*exp(-r*r)/rootPI; 
}
//*********************************************************************
//*********************************************************************

vector<double> stokeslet_density(const vector<double>& pt,
				 const vector<double>& nl) {
  vector<double> f(3,0);

  // TEST: Translating spheroid
  double a_sq = ellipse_a * ellipse_a;
  double b_sq = ellipse_b * ellipse_b;
  double e = sqrt( 1.0 - b_sq / a_sq );
  double e2 = e * e;
  double Le = log( (1.0 + e) / (1.0 - e) );
  
  f[0] = 4.0 * e2 * e * ellipse_a / ellipse_b / ( (1.0 + e2) * Le - 2.0 * e )
    / sqrt( a_sq - e2 * pt[0] * pt[0] );

  return f;
}

//*********************************************************************

void Exact_solution(int N_target, const vector<Target_point>& Target,
		    vector<double>& u_ex) {
  
  // TEST: Translating spheroid
  double a_sq = ellipse_a * ellipse_a;
  double b_sq = ellipse_b * ellipse_b;

  for (int i=0; i<N_target; i++) {
    
    double x = Target[i].x;
    double y = Target[i].y;
    double z = Target[i].z;
    
    double r_sq = y * y + z * z;
    double r = sqrt(r_sq);
    
    double c = sqrt( a_sq - b_sq );
    double R1 = sqrt( (x+c) * (x+c) + r_sq );
    double R2 = sqrt( (x-c) * (x-c) + r_sq );
    
    double R1c = R1 - (x+c);
    double R2c = R2 - (x-c);
    
    double B10;
    if (r_sq < 1.e-14) {
      if ( x<0 ) B10 = log( (x-c) / (x+c) );
      else B10 = log( (x+c) / (x-c) );
    }
    else B10 = log( R2c / R1c );
    
    double B30 = (x+c)/R1 - (x-c)/R2;
    
    double dB1dx = 1.0/R1 - 1.0/R2;
    
    double dB1dyz = 0.0;
    if (r_sq > 1.e-14) {
      dB1dyz = 1.0/R2/R2c - 1.0/R1/R1c;
    }
    double dB3dyz = 1.0/R2 - 1.0/R1 + x * dB1dyz;
    
    double e = c/ellipse_a;
    double e2 = e * e;
    double Le = log( (1.0+e) / (1.0-e) );
    double a1 = e2 / ( (1.0 + e2) * Le - 2.0 * e );
    double b1 = a1 * (1.0 - e2) / (2.0 * e2);
    
    u_ex[3*i  ] = 2.0 * a1 * B10 - a1 * B30 + 2.0 * b1
      *( (x-c)/R2 - (x+c)/R1 + B10 + x * dB1dx );
    double Er = a1 * (1.0/R2-1.0/R1) + 2.0 * b1 * dB3dyz;
    u_ex[3*i+1] = Er * y;
    u_ex[3*i+2] = Er * z;
  } 
}

//*********************************************************************
//*********************************************************************

// Signed distance between two points in 3D

double distance(const vector<double>& x,
		const vector<double>& y,
		const vector<double>& normal) {
  
  return (x[0]-y[0])*normal[0] + (x[1]-y[1])*normal[1] + (x[2]-y[2])*normal[2];
}

//*********************************************************************

double dot_product(const vector<double>& x,
		   const vector<double>& y) {
  
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

//*********************************************************************
// Solve linear system through LU decomposition with pivoting.                
//                                                                           
// Factorize PA = LU with pivoting:                                           
//   The lower and upper triangular matrices are still stored in the original 
//   matrix and the permutation matrix "P" is stored in the vector "int *p".  
//*********************************************************************
 
int LU_factorization(vector< vector<double> >& A,
		     vector<int>& p) {
  int n = A.size();

  for (int j=0; j<n; j++) p[j] = j;

  for (int j=0; j<n; j++) {
    int k = j;
    double m = A[p[j]][j];
    // Search for maximum in this column
    for (int i=j+1; i<n; i++) {
      if (abs(A[p[i]][j]) > abs(m)) {
        m = A[p[i]][j];
        k = i;
      }
    }

    // "Swap" maximum row with current row (using the permutation vector)
    if (k != j) {
      int temp = p[j];
      p[j] = p[k];
      p[k] = temp;
    }

    double ajj = A[p[j]][j];
    if (abs(ajj) < 1.0e-15) {
      return (-1);
    }

    for (int i=j+1; i<n; i++) {
      double lij = A[p[i]][j] / ajj;
      A[p[i]][j] = lij; // lower triangular elements
      for (int k=j+1; k<n; k++) {
        A[p[i]][k] -= lij * A[p[j]][k]; // upper triangular elements
      }
    }
  }
  return 0;
}

//*********************************************************************

vector<double> solveByLU(const vector< vector<double> >& A,
			 const vector<int>& p,
			 const vector<double>& b) {
  int n = A.size();
  vector<double> x(n,0);

  // Solve Ly=b by forward substitution
  x[0] = b[p[0]];
  for (int i=1; i<n; i++) {
    x[i] = b[p[i]];
    double rowsum = 0.0;
    for (int j=0; j<i; j++) {
      rowsum += A[p[i]][j] * x[j];
    }
    x[i] -= rowsum;
  }

  // Solve Ux=y by back substitution
  x[n-1] = x[n-1] / A[p[n-1]][n-1];
  for (int i=n-2; i>=0; i--) {
    double rowsum = 0.0;
    for (int j = n - 1; j > i; j--) {
      rowsum += A[p[i]][j] * x[j];
    }    
    x[i] = (x[i] - rowsum) / A[p[i]][i];
  }

  return x;
}
//*********************************************************************

void initialize_vector(int n, vector<double>& vec) {
  for (int i=0; i<n; i++)  vec[i] = 0.0;
}

//*********************************************************************

void initialize_Matrix(int n, int m,
		       vector< vector<double> >& Mat) {
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      Mat[i][j] = 0.0;
    }
  }   
}

//*********************************************************************
//*********************************************************************

void error_vector(int n,
		  const vector<double>& x_comp,
		  const vector<double>& x_ex,
		  double &err_max,
		  int &ind_max,
		  double &err_l2) {

  err_max = abs(x_comp[0] - x_ex[0]);
  ind_max = 0;
  err_l2 = err_max * err_max;
  
  for (int i=1; i<n; i++) {
    double err_temp = abs(x_comp[i] - x_ex[i]);
    if ( err_temp > err_max) {
      err_max = err_temp;
      ind_max = i;
    }
    err_l2 += err_temp * err_temp;
  }
  err_l2 = sqrt(err_l2 / n);

}
//*********************************************************************

void error_3d_vector(int n,
		     const vector<double>& x_comp,
		     const vector<double>& x_ex,
		     double &err_max,
		     int &ind_max,
		     double &err_l2) {

  double val = (x_comp[0] - x_ex[0]) * (x_comp[0] - x_ex[0])
             + (x_comp[1] - x_ex[1]) * (x_comp[1] - x_ex[1])
             + (x_comp[2] - x_ex[2]) * (x_comp[2] - x_ex[2]);
  
  err_max = sqrt(val);
  ind_max = 0;
  err_l2 = val;
  
  for (int i=1; i<n; i++) {
    double err_temp = (x_comp[3*i] - x_ex[3*i]) * (x_comp[3*i] - x_ex[3*i])
      + (x_comp[3*i+1] - x_ex[3*i+1]) * (x_comp[3*i+1] - x_ex[3*i+1])
      + (x_comp[3*i+2] - x_ex[3*i+2]) * (x_comp[3*i+2] - x_ex[3*i+2]);
    if ( sqrt(err_temp) > err_max) {
      err_max = sqrt(err_temp);
      ind_max = i;
    }
    err_l2 += err_temp;
  }
  err_l2 = sqrt(err_l2 / n);

}

//*********************************************************************

void norms_3d_vector(int n,
		     const vector<double>& x,
		     double &x_max,
		     int &ind_max,
		     double &x_l2) {

  double val = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];  

  x_max = sqrt(val);
  ind_max = 0;
  x_l2 = val;
  
  for (int i=1; i<n; i++) {
    double x_temp = x[3*i] * x[3*i] + x[3*i+1] * x[3*i+1] + x[3*i+2] * x[3*i+2];
    if ( sqrt(x_temp) > x_max) {
      x_max = sqrt(x_temp);
      ind_max = i;
    }
    x_l2 += x_temp;
  }
  x_l2 = sqrt(x_l2 / n);

}

//*********************************************************************

void solve_3by3(double b, const vector<double>& rho,
		double DEL1, double DEL2, double DEL3,
		const vector<double>& S_d1,
		const vector<double>& S_d2,
		const vector<double>& S_d3,
		vector<double>& S_comp) {

  vector<double> lam(3,0);    
  lam[0] = b / DEL1;
  lam[1] = b / DEL2;
  lam[2] = b / DEL3;
  
  vector<double> Line(3,0), RHS(3,0), SOL(3,0);
  vector< vector<double> > Mat(3,Line);
  
  for (int i=0; i<3; i++) {
    Mat[i][0] = 1.0;
    Mat[i][1] = rho[i] * I0(lam[i]);
    Mat[i][2] = rho[i] * rho[i] * rho[i] * I2(lam[i]);
  }
  
  vector<int> p(3,0);
  int result = LU_factorization(Mat, p);
  
  // first component
  RHS[0] = S_d1[0];
  RHS[1] = S_d2[0];
  RHS[2] = S_d3[0];
  SOL = solveByLU(Mat, p, RHS);
  
  S_comp[0] = SOL[0];
  
  // second component
  RHS[0] = S_d1[1];
  RHS[1] = S_d2[1];
  RHS[2] = S_d3[1];
  SOL = solveByLU(Mat, p, RHS);
  
  S_comp[1] = SOL[0];
  
  // third component
  RHS[0] = S_d1[2];
  RHS[1] = S_d2[2];
  RHS[2] = S_d3[2];
  SOL = solveByLU(Mat, p, RHS);
  
  S_comp[2] = SOL[0];
}

//*********************************************************************

double I0(double lam) {
  return exp(-lam*lam)/rootPI - abs(lam) * erfc(abs(lam));
}

double I2(double lam) {
  return 2.0/3.0 * ((0.5-lam*lam) * exp(-lam*lam)/rootPI
		    + abs(lam)*lam*lam * erfc(abs(lam)));
}

