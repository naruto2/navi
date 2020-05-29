#ifndef _EST_XYC2MSH_HPP_
#define _EST_XYC2MSH_HPP_

typedef struct{ int a, b, c, A, B, C; } nde;
typedef struct{ double x, y; char *label; } xyc;

extern void f2mesh(FILE*, std::vector<xyc>&, std::vector<nde>&);
extern long dimp2(std::vector<nde>&N);
extern double delta(int, std::vector<xyc>&, std::vector<nde>&);
#endif
