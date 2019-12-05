#include "../../include/CRModel.h"
#include <iostream>
#include <fstream>
#include <random>

std::mt19937 random_engine;

foodmatrix load_food_matrix(const Metaparameters& m){
  foodmatrix f(m.NS,nvector(m.NR, 0.));
  if(m.verbose > 1){
    std::cout << "  Loading food matrix from " << m.foodmatrixpath << std::endl;
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
