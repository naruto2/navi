#ifndef _EST_NAVIERSTOKES_HPP_
#define _EST_NAVIERSTOKES_HPP_
int navierstokes_init(const char *filename, double Re, double dt,
		      double (*u)(double x, double y),
		      double (*v)(double x, double y));
void navierstokes(matrix<double>&A, std::vector<double>&U,
		  std::vector<double>&b);
void plotuv(std::vector<double>&U);
void fprintuv(std::vector<double>&U);

int islabel(const char *label);
#endif
