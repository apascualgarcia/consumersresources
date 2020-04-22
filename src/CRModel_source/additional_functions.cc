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
    while(std::getline(in, line)){
      /* skip line if it starts with a comment */
      if(line[0]!='#'){
        std::istringstream iss(line);
        input.push_back(nvector());
        unsigned int element;
        while(iss>>element){
          input[index].push_back(element);
        }
        index+=1;
      }
    }
  }
  in.close();

  for (unsigned int x = 0; x < m.NS; x++) {
    for (unsigned int y = 0; y < m.NR; y++) {
      f[x][y]=input[x][y];
    }
  }
  return f;
}

nmatrix load_syntrophy_matrix(const Metaparameters& m){
  nmatrix a(m.NR,nvector(m.NS, 0.));
  nmatrix input;
  //old version, we now allow a syntrophy matrix in whatever path
  //std::string syntrophy_path = optimal_alpha_matrix_path(m.foodmatrixpath);
  std::string syntrophy_path=m.syntrophy_matrix_path;
  if(m.verbose > 2){
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
      if(line[0]!='#'){
        std::istringstream iss(line);
        input.push_back(nvector());
        unsigned int element;
        while(iss>>element){
          input[index].push_back(element);
        }
        index+=1;
      }
    }
  }
  in.close();

  for (unsigned int x = 0; x < m.NR; x++) {
    for (unsigned int y = 0; y < m.NS; y++) {
      a[x][y]=input[x][y];
    }
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
  if(g_path.size()<=4){
    throw error("Invalid path to gamma matrix (name "+g_path+" is too short)");
  }
  for(size_t i=0; i<4; ++i){
    alpha_path.pop_back();
  }
  return alpha_path+"_optimal_alpha.txt";

}
std::string optimal_alpha_matrix_path_from_syntrophy_folder(const Metaparameters& m){
  std::string folder_path=m.syntrophy_matrix_path;
  size_t index_start=m.foodmatrixpath.find_last_of("/");
  size_t length=m.foodmatrixpath.size()-6-index_start+1;
  std::string mat_name=m.foodmatrixpath.substr(index_start+1, length);

  return folder_path+'/'+mat_name+"_optimal_alpha.txt";
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

unsigned int number_of_links(const nmatrix & mat){
  unsigned int links=0;
  for(size_t i=0; i < mat.size(); ++i){
    for(size_t j=0; j < mat[i].size(); ++j){
      if(mat[i][j]*mat[i][j]>0.){
        links+=1;
      }
    }
  }
  return links;
}

ntype connectance(const nmatrix& m){
  unsigned int rows=m.size(), cols=m[0].size();
  ntype connectance=ntype(number_of_links(m))/ntype(rows*cols);
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

nvector operator*(const nmatrix& M, const nvector& v){
  nvector result;
  for(size_t i=0; i < M.size(); ++i){
    ntype value=0.;
    for(size_t j=0; j < M[0].size(); ++j){
      if(M[i].size()!=v.size()){
        throw error("Dimensions do not match for matrix*vector multiplication");
      }
      value+= M[i][j]*v[j];
    }
    result.push_back(value);
  }

  return result;
}

nvector operator*(const nvector& v, const nmatrix& M){
  nvector result;
  if(M.size()!=v.size()){
    throw error("Dimensions do not match for matrix*vector multiplication");
  }
  for(size_t j=0; j < M[0].size(); ++j){
    ntype value=0.;
    for(size_t i=0; i < M.size(); ++i){
      value+= M[i][j]*v[i];
    }
    result.push_back(value);
  }

  return result;
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
  myfile << "# File last opened on " << current_time() <<  std::setprecision(print_precision);
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

bool is_there_coprophagy(const nmatrix& alpha, const nmatrix& gamma){
  unsigned int NR=alpha.size();
  unsigned int NS=alpha[0].size();

  if(gamma.size()!=NS || gamma[0].size()!=NR){
    throw error("Please throw alpha and gamma of appropriate dimensions when checking if there is coprophagy");
  }

  for(size_t mu=0; mu < NR; ++mu){
    for(size_t i=0; i < NS ;++i){
      if(alpha[mu][i]*gamma[i][mu]>0.){
        return false;
      }
    }
  }

  return false;
}

bool has_an_empty_row(const nmatrix& gamma){
  for(size_t i=0; i < gamma.size(); ++i){
    bool row_is_empty=true;
    for(size_t j=0; j < gamma[i].size() && row_is_empty ; ++j){
      if(gamma[i][j]*gamma[i][j]>0){
        row_is_empty=false;
      }
    }
    if(row_is_empty){
      return true;
    }
  }
  return false;
}

bool has_an_empty_column(const nmatrix& gamma){
  return has_an_empty_row(transpose(gamma));
}

std::vector<unsigned int> row_degrees(const nmatrix& M){
  std::vector<unsigned int> degs;
  for(size_t i=0; i < M.size(); ++i){
    unsigned int local_deg=0;
    for(size_t j=0; j < M[i].size();++j){
      if(M[i][j]*M[i][j]>0.){
        local_deg+=1;
      }
    }
    degs.push_back(local_deg);
  }
  return degs;
}

std::vector<unsigned int> columns_degrees(const nmatrix& M){
  return row_degrees(transpose(M));
}

ntype assortativity(const nmatrix & mat){
  std::vector<std::array<size_t ,2>> edges;
  unsigned int Nrows=mat.size(), Ncolumns=mat[0].size();
  ntype ass =0.;
  for(size_t i=0; i < Nrows; ++i){
    for(size_t j=0; j < Ncolumns; ++j){
      if(mat[i][j]*mat[i][j]>ZERO){
        edges.push_back(std::array<size_t,2>({i,j}));
        std::cout << "(" <<i << ","<< j <<")" << std::endl;
      }
    }
  }

  for(size_t i=0; i < edges.size();++i){

  }

  return ass;
}

Eigen::Matrix<ntype, Eigen::Dynamic, Eigen::Dynamic> convert_nmatrix_to_eigen_matrix(const nmatrix& mat){
  Eigen::Matrix<ntype, Eigen::Dynamic, Eigen::Dynamic> eigen_mat;
  unsigned int rows = mat.size(), columns=mat[0].size();
  eigen_mat.resize(rows, columns);

  for(size_t i=0; i < rows; ++i){
    for(size_t j=0; j < columns; ++j){
      eigen_mat(i,j)=mat[i][j];
    }
  }
  return eigen_mat;
}

unsigned int rank(const nmatrix& mat){
  Eigen::Matrix<ntype, Eigen::Dynamic, Eigen::Dynamic> eigen_mat = convert_nmatrix_to_eigen_matrix(mat);
  Eigen::FullPivLU<Eigen::Matrix<ntype, Eigen::Dynamic, Eigen::Dynamic>> lu(eigen_mat);
  return lu.rank();
}

ntype mean_non_zero_elements(const nmatrix & mat){
  ntype mean=0.;
  ntype index=0;
  for(size_t i=0; i < mat.size(); ++i){
    for(size_t j=0; j < mat[i].size(); ++j){
      if(mat[i][j]*mat[i][j]>ZERO){
        mean+=mat[i][j];
        index+=1;
      }
    }
  }
  if(index>=1.){
    return mean/index;
  }
  return mean;
}

nmatrix operator*(const ntype& lambda, const nmatrix& mat){
  nmatrix new_mat=nmatrix(mat.size(), nvector(mat[0].size(), 0.));
  for(size_t i=0; i < new_mat.size();++i){
    for(size_t j=0; j < new_mat[i].size(); ++j){
      new_mat[i][j]=lambda*mat[i][j];
    }
  }
  return new_mat;
}

nmatrix operator*(const nmatrix& mat, const ntype& lambda){
  return lambda*mat;
}

nmatrix operator/(const nmatrix& mat, const ntype& lambda){
  return (1./lambda)*mat;
}

nmatrix binary_matrix_no_intraspecific_syntrophy(const nmatrix& g){
  nmatrix a=random_binary_matrix_with_connectance(g[0].size(), g.size(), connectance(g));
  std::uniform_int_distribution<unsigned int> mu_distrib(0, g[0].size()-1);
  std::uniform_int_distribution<unsigned int> i_distrib(0, g.size()-1);
  while(is_there_coprophagy(a,g)){
    for(size_t mu=0; mu < a.size(); ++mu){
      for(size_t i=0; i < a[0].size(); ++i){
        /* if there is coprophagy, we move the problematic links */
        if(a[mu][i]*g[i][mu]>0){
          unsigned int new_mu = mu_distrib(random_engine);
          unsigned int new_i = i_distrib(random_engine);
          while(g[new_i][new_mu]>0){
            new_mu = mu_distrib(random_engine);
            new_i = i_distrib(random_engine);
          }
          a[mu][i]=0;
          a[new_mu][new_i]=1;
        }
      }
    }
  }
  return a;
}


nmatrix random_binary_matrix_with_connectance(const unsigned int& rows, const unsigned int& columns, const ntype& conn){
  nmatrix mat(rows, nvector(columns, 0.));
  std::uniform_real_distribution<ntype> unif_distrib(0., 1.);
  for(size_t i=0; i < rows; ++i){
    for(size_t j=0; j < columns;++j){
      if(unif_distrib(random_engine)< conn){
        mat[i][j]=1;
      }
    }
  }

  return mat;
}


nmatrix build_LRI_matrix(const nmatrix& g, const Metaparameters& m, const ntype& target_conn){

  Metaparameters metaparams=m;
  MonteCarloSolver mcsolv;
  mcsolv.T=5.;
  mcsolv.max_steps=1000000;
  mcsolv.max_fails=1000;
  mcsolv.annealing_freq=1000;
  mcsolv.annealing_const=1.-1e-2;
  mcsolv.display_stride=10000;
  mcsolv.cost_function=quadratic_form_low_intra_resource_interaction;
  mcsolv.additional_params=&metaparams;

  metaparams.alpha0=metaparams.NR*metaparams.sigma0*metaparams.R0*metaparams.gamma0;
  bool allow_coprophagy=true;
  foodmatrix alpha = optimal_syntrophy_from_consumption(g, allow_coprophagy, mcsolv, target_conn);
  return alpha;
}
