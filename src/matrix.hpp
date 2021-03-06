#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

template<typename Real>
class matrix : public std::vector< std::map<long,Real> > {

public:
  matrix()       : std::vector< std::map<long, Real> >(){}
  matrix(long n) : std::vector< std::map<long, Real> >(n){}
};

int solver(matrix<double>&A,std::vector<double>&x,const std::vector<double>&b);
void plotmatrix(matrix<double>&A);
void printmatrix(matrix<double>&A);
void printvector(const std::vector<double>&b);
#endif
