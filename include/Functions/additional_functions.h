#ifndef ADDITIONAL_FUNCTIONS_H
#define ADDITIONAL_FUNCTIONS_H

#include "../Classes/Metaparameters.h"
#include "../Classes/Custom_types.h"
#include "../Classes/Extinction.h"
#include<algorithm>
#include<random>
#include<string>
#include<iostream>

/* loads the food matrix and relabels columns and resources such that gamma is "most triangular" */
foodmatrix load_food_matrix(const Metaparameters&);
nmatrix load_syntrophy_matrix(const Metaparameters&);
/* relabels the row of  m such that degree of row i > degree of row j if i < j */
nmatrix order_matrix_by_row_degree(const nmatrix&);
/* relabels the columns of m such that degree of column i > degree of column j if i < j */
nmatrix order_matrix_by_column_degree(const nmatrix& m);
/* returns path of optimal alpha matrix for a given path to the gamma matrix */
/*    use this one when path_to_syntrophy_matrix is a path to file */
std::string optimal_alpha_matrix_path(const std::string&);
/*    and this one when it is a path to a folder */
std::string optimal_alpha_matrix_path_from_syntrophy_folder(const Metaparameters& m);
ntype norm(const nvector&);
ntype mean(const nmatrix &);
ntype det(const nmatrix&);

/* computes number of links in a matrix */
unsigned int number_of_links(const nmatrix&);

std::vector<unsigned int> row_degrees(const nmatrix& M);
std::vector<unsigned int> columns_degrees(const nmatrix& M);

/* here we take connectance as the number of links divided by the number of total possible links */
ntype connectance(const nmatrix&);
/* here we take the nestedness as defined in Bastolla's 2013 paper */
ntype nestedness(const nmatrix &);
ntype trace(const nmatrix&);
nmatrix random_uniform_matrix(const unsigned int&, const unsigned int&, const ntype&);
void rescale_mean(nmatrix&, const ntype&);

bool is_there_coprophagy(const nmatrix& alpha, const nmatrix& gamma);
/* returns true if gamma has one row filled with zeros only */
bool has_an_empty_row(const nmatrix& gamma);
/* returns true if gamma has one column filled with zeros only */
bool has_an_empty_column(const nmatrix& gamma);

nmatrix operator+(const nmatrix&, const nmatrix&);
nmatrix operator-(const nmatrix&, const nmatrix &);
nmatrix operator*(const nmatrix&, const nmatrix&);
nmatrix operator-(const nmatrix&);

nmatrix transpose(const nmatrix & m);

bool non_neg_elements(const nmatrix&);
bool non_neg_elements(const nvector&);

void initialize_random_engine(const Metaparameters&);
void print_rand_number();

ntype mean(const nvector&);
ntype standard_dev(const nvector&);
ntype median(const nvector&);

nvector linear_interval(const ntype& begin, const ntype& end, unsigned int Npoints);
nvector log_interval(const ntype& begin, const ntype& end, unsigned int Npoints);

ntype angle(const nvector&, const nvector&);
ntype distance(const nvector&, const nvector&);

nvector operator+(const nvector&, const nvector&);
nvector operator-(const nvector&);
nvector operator-(const nvector&, const nvector&);
ntype operator*(const nvector&, const nvector&);
nvector operator*(const nmatrix&, const nvector&);

/* /!\ THIS CORRESPONDS TO v^T M for v a vector */
nvector operator*(const nvector& v, const nmatrix& M);


/*  returns a vector of string containing the path of the matrices (loads
    path_to_list file) */
std::vector<std::string> load_food_matrix_list(std::string path_to_list);

std::string current_time();

/* checks is number is an error */
bool is_an_error(ntype);

/* open external file */
std::ofstream open_external_file_append(std::string);
std::ofstream open_external_file_truncate(std::string);


bool compare_complex(const nctype&, const nctype&);
bool operator<(const nctype&, const nctype&);
bool operator>(const nctype&, const nctype&);


/* taken from https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes on Feb 26 2020 */
template <typename T>
std::vector<size_t> sort_indices(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

template<typename T>
T minimum(const std::vector<T> &vec){
  return *std::min_element(vec.begin(), vec.end());
};

template<typename T>
T maximum(const std::vector<T> & vec){
  return *std::max_element(vec.begin(), vec.end());
};


#endif
