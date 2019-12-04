#include "../../include/CRModel.h"
#include <iostream>
#include <fstream>
#include <iomanip>

std::ostream& operator<<(std::ostream& os, const Metaparameters& m){
  os << "gamma0 = " << m.gamma0 << std::endl;
  os << "alpha0 = " << m.alpha0 << std::endl;
  os << "sigma0 = " << m.sigma0 << std::endl;
  os << "p = " << m.p << std::endl;
  os << "R0 = " << m.R0 << std::endl;
  os << "S0 = " << m.S0 << std::endl;
  os << "NR = " << m.NR << std::endl;
  os << "NS = " << m.NS << std::endl;
  return os;
}
std::ostream& operator<<(std::ostream& os, const nctype& number){
  os << std::real(number);
  if(std::imag(number)>0.){
    os << "+" << std::imag(number)<<"j";
  }else if(std::imag(number)<0.){
    os << std::imag(number)<<"j";
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const ncvector& vec){
  for(size_t i = 0; i <vec.size(); ++i){
    os << vec[i] << " ";
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const CRModel& M){
  M.display(os);
  return os;
}
std::ostream& operator<<(std::ostream& os, const Parameter_set& p){
  display_matrix_w_name(os, "gamma", p.gamma);
  //os << "mean = " << mean(p.gamma) << std::endl;
  display_matrix_w_name(os, "alpha", p.alpha);
  //os << "mean = " << mean(p.alpha) << std::endl;
  display_matrix_w_name(os, "tau", p.tau);
  //os << "mean = " << mean(p.tau) << std::endl;
  display_matrix_w_name(os, "sigma", p.sigma);
  //os << "mean = " << mean(p.sigma) << std::endl;
  display_vector_w_name(os, "l", p.l);
  display_vector_w_name(os, "m", p.m);
  display_vector_w_name(os, "d", p.d);
  os << "NR = " << p.NR << std::endl;
  os << "NS = " << p.NS << std::endl;
  return os;
}
std::ostream& operator<<(std::ostream& os, const Model_parameters& M){
  os << *(M.get_parameters()) << std::endl;
  return os;
}
std::ostream& operator<<(std::ostream& os, const ntensor& T){
  for(size_t i = 0; i < T.size(); ++i){
    os << "Element " << i << " of tensor :" << std::endl << T[i] << std::endl;
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const nvector& v){
  for (size_t i=0; i < v.size(); ++i){
    os << std::right << std::fixed  << std::setprecision(print_precision) << v[i] << " ";
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const nmatrix& M){
  for (size_t i=0; i < M.size(); ++i){
    os << M[i] << std::endl;
  }
  return os;
}

std::ostream& display_matrix_w_name(std::ostream& os, std::string mat_name, const nmatrix & mat){
  os << mat_name << " = " << mat[0] << std::endl;
  for(size_t i = 1; i < mat.size(); ++i){
    os << std::setw(mat_name.size()+4+4) <<std::setfill(' ')<< mat[i] << std::endl;
  }
  return os;
}
std::ostream& display_vector_w_name(std::ostream& os, std::string vec_name, const nvector& vect){
  os << vec_name << " = " << vect << std::endl;
  return os;
}
