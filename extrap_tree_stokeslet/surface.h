#include <iostream>
#include <vector>
#include <cmath>

#include "utilities.h"

using namespace std;

void Generate_Surface(int N, double h, vector<Surf_point>* Surfc);
int sign_phi(const vector<double>& x);
double Find_hypersurf_pt(const vector<double>& A, double H_seg, int i);
double Newton(const vector<double>& PT, double H_seg,
	      int i, double a, double b);
double bisection(const vector<double>& PT, double H_seg,
		 int i, double a, double b);
double Part_Unity(int i, const vector<double>& Nrml);
double bump(double r);

int find_Nearest(const vector<double>& pt,
		 int N_quad, const vector<Surf_point>& Surfc,
		 const vector<double>& x0,
		 Surf_point& nrst);

void Generate_Targets(double h,
		      int N_quad, const vector<Surf_point>& Surfc,
		      vector<Target_point>* Target);
