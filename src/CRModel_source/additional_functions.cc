#include "../../include/CRModel.h"
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <algorithm>

std::mt19937 random_engine;

foodmatrix load_food_matrix(const Metaparameters& m){
  foodmatrix f(m.NS,nvector(m.NR, 0.));
  if(m.verbose > 1){
    std::cout << "\t Loading food matrix from " << m.foodmatrixpath << std::endl;
  }
  std::ifstream in(m.foodmatrixpath);
  if (!in) {
    std::cerr << "Cannot open file for the matrix " << m.foodmatrixpath << std::endl;
    std::cerr << "Now aborting the simulation" << std::endl;
    abort();
    return f;
  }
  for (unsigned int x = 0; x < m.NS; x++) {
    for (unsigned int y = 0; y < m.NR; y++) {
      in >> f[x][y];
    }
  }
  in.close();
  return f;
}
ntype norm(const nvector& v){
  ntype a=0.;
  for (size_t i =0 ; i < v.size(); ++i){
    a+=pow(v[i], 2.0);
  }
  return pow(a, 0.5);
}
ntype mean(const nmatrix& m){
  ntype mean(0.);
  size_t R = m.size(), C = m[0].size();
  for(size_t i=0; i < R; ++i){
    for(size_t j=0; j < C; ++j){
      mean+=(m[i][j]/(C*R));
    }
  }
  return mean;
}
nmatrix random_uniform_matrix(const unsigned int& R, const unsigned int& C, const ntype& mean_){
  nmatrix mat(R, nvector(C, 0.));
  std::uniform_real_distribution<ntype> unif_distrib(0., 1.);
  ntype total=0.;
  for(size_t i=0; i < R; ++i){
    for(size_t j=0; j < C; ++j){
      mat[i][j] = unif_distrib(random_engine);
      total+=mat[i][j];
    }
  }


  for(size_t i=0; i < R; ++i){
    for(size_t j=0; j < C; ++j){
      mat[i][j] *=(R*C*mean_/total);
    }
  }
  return mat;
}
void rescale_mean(nmatrix& M, const ntype& mean){
  ntype total=0.;
  for(size_t i=0; i < M.size(); ++i){
    for(size_t j=0; j < M[i].size(); ++j){
      total+=M[i][j];
    }
  }
  if(total!=0.){
    for(size_t i =0; i < M.size(); ++i){
      for(size_t j=0; j < M[i].size(); ++j){
        M[i][j] *= (M.size()*M[i].size()*mean/total);
      }
    }
  }else{
    std::cerr << "Error, matrix is zero" << std::endl;
  }
  return;
}

bool non_neg_elements(const nmatrix& m){
  for(size_t i = 0; i < m.size();++i){
    for(size_t j=0; j< m[i].size();++j){
      if(m[i][j] < 0.){
        return false;
      }
    }
  }
  return true;
}
bool non_neg_elements(const nvector& v){
  for(size_t i =0; i < v.size();++i){
    if(v[i] < 0.){
      return false;
    }
  }
  return true;
}

void initialize_random_engine(const Metaparameters& m){
  random_engine.seed(m.seed_number);
  return;
}
void print_rand_number(){
  std::uniform_real_distribution<ntype> dist(0., 1.);
  std::cout << dist(random_engine) << std::endl;
}

ntype mean(const nvector& x){
  ntype mean = 0.;
  size_t Nsimul = x.size();
  for(size_t i = 0; i < Nsimul; ++i){
    mean += ntype(x[i]/Nsimul);
  }
  return mean;
}

ntype median(const nvector& x){
  ntype median = 0.;
  nvector sort_x = x;
  std::sort(sort_x.begin(), sort_x.end());
  unsigned int N = x.size();
  if(N%2==0){
    median = 0.5*(sort_x[N/2-1]+sort_x[N/2]);
  }else{
    median = sort_x[(N-1)/2];
  }
  return median;
}

ntype standard_dev(const nvector& x){
  ntype std = 0.;
  ntype variance = 0.;
  ntype m = mean(x);
  size_t N = x.size();
  for(size_t i=0; i < N; ++i){
    variance += ntype(1./(N-1))*(x[i]-m)*(x[i]-m);
  }
  std=sqrt(variance);
  return std;
}

ntype angle(const nvector& v1, const nvector& v2){
  ntype costheta = 0.;
  costheta = (v1*v2)/(norm(v1)*norm(v2));
  return acos(costheta);
}

ntype distance(const nvector& v1, const nvector& v2){
  return norm(v1-v2);
}

nvector operator+(const nvector& v1, const nvector& v2){
  nvector diff;
  size_t N1 = v1.size(), N2 = v2.size(), N=N1;
  if(N1!=N2){
    std::cerr << "Trying to add vectors of two different sizes " << std::endl;
    std::cerr << "Only the first common components will be taken into account" << std::endl;
    if(N2 < N1){
      N = N2;
    }
  }

  for(size_t i = 0; i < N; ++i){
    diff.push_back(v1[i]+v2[i]);
  }
  return diff;
}

nvector operator-(const nvector& v1){
  nvector opp;
  for(size_t i=0; i < v1.size(); ++i){
    opp.push_back(-v1[i]);
  }
  return opp;
}

statistics::statistics(const nvector& v){
  this->mean_ = mean(v);
  this->std_deviation_ = standard_dev(v);
  this->median_ = median(v);
}

statistics::statistics(const statistics & s){
  this->mean_ = s.mean_;
  this->std_deviation_=s.std_deviation_;
  this->median_ = s.median_;
}

statistics::statistics(){
  this->mean_ = NULL;
  this->std_deviation_ = NULL;
  this->median_ = NULL;
}

nvector operator-(const nvector& v1, const nvector& v2){
  return v1+(-v2);
}

ntype operator*(const nvector& v1, const nvector& v2){
  ntype dot_prod=0.;
  for(size_t i=0; i < v1.size(); ++i){
    dot_prod += v1[i]*v2[i];
  }
  return dot_prod;
}