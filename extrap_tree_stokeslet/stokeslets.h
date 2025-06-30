#include <iostream>
#include <vector>
#include <cmath>

#include "KITC.h"

using namespace std;

void eval_Stokes_SL_offSurf_3del(double h, const vector<double>& rho,
				 double DEL1, double DEL2, double DEL3,
				 int N_quad, const vector<Surf_point>& Surfc,
				 int N_target, const vector<Target_point>& Target,
				 vector<double>& SL);
