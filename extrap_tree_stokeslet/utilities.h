#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

static const int N_gridlines = 128;
static const int del_h = 1; // 1:    del=rho*h,       rho = (3,4,5)
                            // else: del=rho*h^(4/5), rho=(2,3,4)*(1/64)^(1/5)
static const int use_tree = 1; //1: use tree structure to localize calculations
static const int use_tree_approximate = 1; //1: use treecode to approximate well-separated evaluations

static const double PI = 3.14159265358979323846;
static const double rootPI = sqrt(PI);
static const double theta = 70.0*PI/180.0; //angle in Beale et al (2.2)
static const double tol = 1e-14; //tolerance in the search of quadrature points (Newton's/bisection method)
static const double ellipse_a = 1.0;
static const double ellipse_b = 0.5;
static const double ellipse_c = 0.5;

using namespace std;

struct Surf_point
{
  Surf_point() : x(3,0), f(3,0), Nrml(3,0), Area(0) {}
  vector<double> x;
  vector<double> f;  // usually used for Stokeslet density
  vector<double> Nrml;
  double Area;
};

struct Target_point
{
  Target_point() : x(0), y(0), z(0), nrst_x(0), nrst_y(0), nrst_z(0), b(0), nrst_Nrml(3,0), nrst_f(3,0) {}
  double x, y, z;
  double nrst_x, nrst_y, nrst_z, b;
  vector<double> nrst_Nrml, nrst_f;
};

double phi(const vector<double>& x);
double Dphi(int i, const vector<double>& x);
void D2phi(const vector<double>& x, double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23, double& phi31,
	   double& phi32, double& phi33);

void H_gauss(double r, double d, double& H1, double& H2);
void s_opt2(double r, double& s1, double& s2);

vector<double> stokeslet_density(const vector<double>& pt,
				 const vector<double>& nl);

void Exact_solution(int N_target, const vector<Target_point>& Target,
		    vector<double>& u_ex);

double distance(const vector<double>& x,
		const vector<double>& y,
		const vector<double>& normal);

double dot_product(const vector<double>& x,
		   const vector<double>& y);

int LU_factorization(vector< vector<double> >& A,
		     vector<int>& p);

vector<double> solveByLU(const vector< vector<double> >& A,
			 const vector<int>& p,
			 const vector<double>& b);

void initialize_vector(int n, vector<double>& vec);
void initialize_Matrix(int n, int m,
		       vector< vector<double> >& Mat);

void error_vector(int n,
		  const vector<double>& x_comp,
		  const vector<double>& x_ex,
		  double &err_max,
		  int &ind_max,
		  double &err_l2);

void error_3d_vector(int n,
		     const vector<double>& x_comp,
		     const vector<double>& x_ex,
		     double &err_max,
		     int &ind_max,
		     double &err_l2);

void norms_3d_vector(int n,
		     const vector<double>& x,
		     double &x_max,
		     int &ind_max,
		     double &x_l2);

void solve_3by3(double b, const vector<double>& rho,
		double DEL1, double DEL2, double DEL3,
		const vector<double>& S_d1,
		const vector<double>& S_d2,
		const vector<double>& S_d3,
		vector<double>& S_comp);

double I0(double lam);
double I2(double lam);

