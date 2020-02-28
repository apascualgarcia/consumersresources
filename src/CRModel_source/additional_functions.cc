#include "../../include/CRModel.h"
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>

std::mt19937 random_engine;

foodmatrix load_food_matrix(const Metaparameters& m){
  foodmatrix f(m.NS,nvector(m.NR, 0.));
  nmatrix input;
  if(m.verbose > 1){
    std::cout << "\t Loading food matrix from " << m.foodmatrixpath << std::endl;
  }
  std::ifstream in(m.foodmatrixpath);
  if (!in) {
    error err("Cannot open file for the matrix "+m.foodmatrixpath);
    throw err;
  }
  if(in.good()){
    std::string line;
    unsigned int index=0;
    while(getline(in, line)){
      std::istringstream iss(line);
      input.push_back(nvector());
      unsigned int element;
      while(iss>>element){
        input[index].push_back(element);
      }
      index+=1;
    }
  }
  in.close();

  for (unsigned int x = 0; x < m.NS; x++) {
    for (unsigned int y = 0; y < m.NR; y++) {
      f[x][y]=input[x][y];
    }
  }

  if(m.verbose > 1){
    std::cout << "\t Note that the food matrix taken is " << m.foodmatrixpath << " but rearranged in the optimal way" << std::endl;
  }

  return foodmatrix(order_matrix_by_column_degree(order_matrix_by_column_degree(f)));
}

nmatrix load_syntrophy_matrix(const Metaparameters& m){
  nmatrix a(m.NR,nvector(m.NS, 0.));
  nmatrix input;
  std::string syntrophy_path = optimal_alpha_matrix_path(m.foodmatrixpath);
  if(m.verbose > 1){
    std::cout << "\t Loading alpha matrix from " << syntrophy_path << std::endl;
  }
  std::ifstream in(syntrophy_path);
  if (!in) {
    error err("Cannot open file for the alpha matrix "+syntrophy_path);
    throw err;
  }
  if(in.good()){
    std::string line;
    unsigned int index=0;
    while(getline(in, line)){
      std::istringstream iss(line);
      input.push_back(nvector());
      unsigned int element;
      while(iss>>element){
        input[index].push_back(element);
      }
      index+=1;
    }
  }
  in.close();

  for (unsigned int x = 0; x < m.NR; x++) {
    for (unsigned int y = 0; y < m.NS; y++) {
      a[x][y]=input[x][y];
    }
  }

  if(m.verbose > 1){
    std::cout << "\t Note that the syntrophy matrix taken is " << syntrophy_path << " but rearranged in the optimal way" << std::endl;
  }

  return a;
}
nmatrix order_matrix_by_row_degree(const nmatrix& m){
  nmatrix ordered_mat;

  unsigned int rows=m.size(), cols=m[0].size();
  std::vector<unsigned int> row_degs(rows, 0);

  for(size_t i=0; i < rows; ++i){
    for(size_t j=0; j < cols; ++j){
      if(m[i][j]*m[i][j]>0.){
        row_degs[i]+=1;
      }
    }
  }

  /* contains the inverted sorted indices of row_degs i.e. sorted_indices[rows-1] is the row with the largest degree and so on*/
  std::vector<size_t> sorted_indices=sort_indices(row_degs);

  for(size_t i=0; i < rows; ++i){
    ordered_mat.push_back(m[sorted_indices[rows-i-1]]);
  }

  return ordered_mat;
}
nmatrix order_matrix_by_column_degree(const nmatrix& m){
  return transpose(order_matrix_by_row_degree(transpose(m)));
}
std::string optimal_alpha_matrix_path(const std::string& g_path){
  /* IT IS IMPORTANT THAT THE FILE EXTENSION IS .txt */
  std::string alpha_path=g_path;
  if(g_path.size()<=3){
    throw error("Invalid path to gamma matrix (name is too short)");
  }
  for(size_t i=0; i<3; ++i){
    alpha_path.pop_back();
  }
  return alpha_path+"_optimal_alpha.txt";

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
ntype det(const nmatrix& m){
  Eigen::MatrixXf mat(m.size(), m[0].size());
  for(size_t i=0; i<m.size(); ++i){
    for(size_t j=0; j < m[0].size();++j){
      mat(i,j)=m[i][j];
    }
  }
  return ntype(mat.determinant());
}
ntype connectance(const nmatrix& m){
  ntype connectance=0.;
  unsigned int rows=m.size(), cols=m[0].size();
  for(size_t i=0; i < rows; ++i){
    for(size_t j=0; j < rows; ++j){
      if(m[i][j]*m[i][j]>0.){
        connectance+=1.;
      }
    }
  }
  connectance/=(rows*cols);
  return connectance;
}
ntype nestedness(const nmatrix& mat){
  ntype eta=0.;
  unsigned int rows = mat.size(), cols=mat[0].size();

  /* first build binary matrix out of mat */
  nmatrix m(rows, nvector(cols, 0.));
  std::vector<unsigned int> degrees;

  /* and get the degree distribution */
  for(size_t i=0; i < rows; ++i){
    unsigned int local_degree=0.;
    for(size_t j=0; j < cols; ++j){
      if(mat[i][j]*mat[i][j]>0){
        local_degree+=1;
        m[i][j]=1;
      }
    }
    degrees.push_back(local_degree);
  }

  nmatrix overlap=m*transpose(m);
  ntype eta_num=0., eta_denom=0.;
  for(size_t i=0; i < rows;++i){
    for(size_t j=0; j < i; ++j){
      eta_num+=overlap[i][j];
      if(degrees[i]>degrees[j]){
        eta_denom+=degrees[j];
      }else{
        eta_denom+=degrees[i];
      }
    }
  }

  eta=eta_num/eta_denom;

  return eta;
}
ntype trace(const nmatrix& m){
  ntype trace=0.;
  for(size_t i=0; i < m.size() || i < m[0].size(); ++i){
    trace+=m[i][i];
  }
  return trace;
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

nmatrix operator+(const nmatrix & A , const nmatrix& B){
  unsigned int rows = A.size(), cols = A[0].size();
  if((rows!=B.size())||(cols!=B[0].size())){
    error e("Impossible to add matrices, the dimensions do not match.");
    throw e;
  }
  nmatrix Add(rows, nvector(cols, 0.));
  for(size_t i=0; i < rows; ++i){
    for(size_t j=0; j < cols; ++j){
      Add[i][j] = A[i][j]+B[i][j];
    }
  }
  return Add;
}
nmatrix operator-(const nmatrix& A){
  nmatrix Opp(A.size(), nvector(A[0].size(), 0.));
  for(size_t i=0; i < A.size(); ++i){
    for(size_t j=0; j < A[0].size();++j){
      Opp[i][j]= -A[i][j];
    }
  }
  return Opp;
}

nmatrix operator-(const nmatrix& A, const nmatrix& B){
  return A+(-B);
}

nmatrix operator*(const nmatrix& A, const nmatrix& B){
  unsigned int rows_A = A.size(), cols_A = A[0].size();
  unsigned int rows_B = B.size(), cols_B = B[0].size();

  nmatrix P(rows_A, nvector(cols_B, 0.));

  if(cols_A!=rows_B){
    error e("Error in the matrix multiplication, the dimensions do not match.");
    throw e;
  }

  for(size_t i=0; i < rows_A; ++i){
    for(size_t j=0; j < cols_B; ++j){
      for(size_t k=0; k < rows_B; ++k){
        P[i][j] += A[i][k]*B[k][j];
      }
    }
  }

  return P;
}

nmatrix transpose(const nmatrix& m){
  unsigned int rows_m = m.size(), cols_m = m[0].size();
  nmatrix transp(cols_m, nvector(rows_m, 0.));
  for(size_t i=0; i<cols_m; ++i){
    for(size_t j=0; j < rows_m; ++j){
      transp[i][j] = m[j][i];
    }
  }
  return transp;
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
  if(m.verbose > 0){
    std::cout << "Initiated random engine with seed " << m.seed_number << std::endl;
  }
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

nvector linear_interval(const ntype& begin, const ntype& end, unsigned int Npoints){
  nvector interval(Npoints, 0.);
  if(Npoints==1){
    interval[0]=begin;
  }else{
    for(size_t i=0; i < Npoints; ++i){
      interval[i] = begin+(end-begin)*i/(Npoints-1);
    }
  }
  return interval;
}

nvector log_interval(const ntype& begin, const ntype& end, unsigned int Npoints){
  nvector interval(Npoints,0.);
  nvector log_interval=linear_interval(log10(begin), log10(end), Npoints);
  for(size_t i=0; i < Npoints; ++i){
    interval[i]=pow(10., log_interval[i]);
  }
  return interval;
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


writemode::writemode(bool write_, std::ostream& os):write(write_),write_path(os){}
writemode::writemode():write(false), write_path(std::cout){}


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


std::vector<std::string> load_food_matrix_list(std::string path_to_list){
  std::vector<std::string> matrices;
  std::ifstream in(path_to_list);
  if (!in) {
    std::cerr << "Cannot open file containing the list of matrices " << path_to_list << std::endl;
  }else{
    do{
      std::string a;
      in >> a;
      matrices.push_back(a);
    }while(!in.eof());
    /* remove last string if white space */
    std::string str = matrices[matrices.size()-1];
    if(str.find_first_not_of(' ') == std::string::npos){
      matrices.pop_back();
    }
  }
  in.close();
  return matrices;
}

void error::handle(){
  //std::cerr << (*this) << std::endl;
  switch(this->category){
    case 0:{
      std::cerr << *this << "This will abort the simulation." << std::endl;
      abort();
      break;
    }
    case 1:{
      std::cerr << *this ;
      std::cerr << "The simulation will keep running but please keep that in mind when interpreting results." << std::endl;
      break;
    }
    default:{
      break;
    }
  }
  return;
}

bool is_an_error(ntype a){
  return isnan(a);
}

std::ofstream open_external_file_append(std::string path){
  std::ofstream myfile;
  myfile.open(path, std::ios::app);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << path << std::endl;
    error err("Could not open external file "+path);
    throw err;
  }
  myfile << "# File opened on " << current_time() <<  std::setprecision(print_precision);
  return myfile;
}

std::ofstream open_external_file_truncate(std::string path){
  std::ofstream myfile;
  myfile.open(path, std::ios::trunc);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << path << std::endl;
    std::cerr << "Could not open " << path << std::endl;
    error err("Could not open external file "+path);
    throw err;
  }
  myfile << "# File created on " << current_time() << std::setprecision(print_precision);
  return myfile;
}

std::string current_time(){
  auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::string time_now(std::ctime(&timenow));
  return time_now;
}

bool compare_complex(const nctype& a, const nctype& b){
  return(real(a)<real(b));
}

bool operator<(const nctype& a, const nctype& b){
  return(real(a)<real(b));
}
bool operator>(const nctype& a, const nctype& b){
  return(real(a)>real(b));
}

std::ostream& display_food_matrix(std::ostream& os, const foodmatrix& f){
  os << std::fixed << std::setprecision(0);
  for(size_t i=0; i < f.size(); ++i){
    for(size_t j=0; j < f[0].size();++j){
      os << f[i][j] << " ";
    }
    os << std::endl;
  }
  return os;
}
