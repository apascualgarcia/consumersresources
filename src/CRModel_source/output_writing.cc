#include "../../include/CRModel.h"
#include <iostream>
#include <fstream>
#include <iomanip>

std::ostream& operator<<(std::ostream& os, const Metaparameters& m){
  os << "gamma0 = " << m.gamma0 << "; ";
  os << "alpha0 = " << m.alpha0 << "; ";
  os << "sigma0 = " << m.sigma0 << "; ";
  os << "l0 = " << m.l0 << ";";
  os << "p = " << m.p << "; ";
  os << "R0 = " << m.R0 << "; ";
  os << "S0 = " << m.S0 << "; ";
  os << "NR = " << m.NR << "; ";
  os << "NS = " << m.NS  << "; ";
  os << "matrix path = " << m.foodmatrixpath <<"; ";
  os << "convergence threshold = " << m.convergence_threshold;
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
  os << "NS = " << p.NS ;
  return os;
}

std::ostream& operator<<(std::ostream& os, const Model_parameters& M){
  M.display(os);
  return os;
}
std::ostream& operator<<(std::ostream& os, const ntensor& T){
  for(size_t i = 0; i < T.size()-1; ++i){
    os << "Element " << i << " of tensor :" << std::endl << T[i] << std::endl;
  }
  os << "Element " << T.size()-1 << " of tensor :" << std::endl << T[T.size()-1];
  return os;
}
std::ostream& operator<<(std::ostream& os, const nvector& v){
  for (size_t i=0; i < v.size(); ++i){
    os << std::scientific << std::setprecision(print_precision) << v[i] << " ";
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const nmatrix& M){
  for (size_t i=0; i < M.size()-1; ++i){
    os << M[i] << std::endl;
  }
  os << M[M.size()-1];
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
std::ostream& operator<<(std::ostream& os, const stability_metrics& m){
  os << std::endl;
  os << "Resilience : " << m.resilience<< std::endl;
  os << "Number of extinctions : " << m.extinctions << std::endl;
  os << "Angle between equilibria : " << m.angle_between_equilibria<< std::endl;
  os << "Distance between equilibria : " << m.distance_between_equilibria;
  return os;
}

std::ostream& operator<<(std::ostream& os, const statistics& stats){
  os << stats.mean_ << "+/-" << stats.std_deviation_ << " (median : "<<stats.median_ << ")";
  return os;
}


std::ostream& operator<<(std::ostream& os, const stabilitymode& stab){
  switch(stab){
    case dynamical:{
      os << "dynamical";
      break;
    }
    case structural:{
      os << "structural";
      break;
    }
    default:
      break;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const interval & interval_){
  os <<  "[" << interval_.begin << ";" << interval_.end << "]";
  return os;
}

std::ostream& operator<<(std::ostream& os, const fitmode & fit){
  switch(fit){
    case sigmoidal:
      os << "sigmoidal";
      break;
    case sigmoidal_erf:
      os << "sigmoidal erf";
      break;
    case polynomial:
      os << "polynomial";
      break;
    default:
      os << "unknown fit";
      break;
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const systemstability& stab){
  switch(stab){
    case stable:{
      os << "stable";
      break;
    }
    case marginal:{
      os << "marginally stable";
      break;
    }
    case unstable:{
      os << "unstable";
      break;
    }
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const stability& stab){
  os << stab.unstable << " " << stab.marginally_stable << " " << stab.stable << " ";
  return os;
}

std::ostream& operator<<(std::ostream& os, const error& e){
  os << "Error during runtime : " << e.message << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const buildingmode & b){
  switch(b){
    case use_l:{
      os << "use_l";
      break;
    }
    case use_m:{
      os << "use_m";
      break;
    }
    default:{
      break;
    }
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const alphamode& a){
  switch(a){
    case random_structure:{
      os << "random structure";
      break;
    };

    case no_release_when_eat:{
      os << "no intraspecific syntrophy";
      break;
    };

    case optimal_matrix:{
      os << "optimized matrix";
      break;
    };

    case fully_connected:{
      os << "fully connected";
      break;
    };

    default:{
      throw error("This mode of alpha_mode has not been implemented yet in the << operator",1);
      break;
    }
  }
  return os;
}

std::ostream& display_food_matrix(std::ostream& os, const foodmatrix& f){
  for(size_t i=0; i < f.size(); ++i){
    for(size_t j=0; j < f[i].size(); ++j){
      if(f[i][j]>0){
        os << "1 ";
      }else{
        os << "0 ";
      }
    }
    if(i < f.size()-1){
      os << std::endl;
    }
  }
  return os;
}
