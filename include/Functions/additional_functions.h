#ifndef ADDITIONAL_FUNCTIONS_H
#define ADDITIONAL_FUNCTIONS_H

#include "../Classes/Metaparameters.h"
#include "../Classes/Custom_types.h"
#include "../Classes/Extinction.h"
#include<random>
#include<string>
#include<iostream>

foodmatrix load_food_matrix(const Metaparameters&);
ntype norm(const nvector&);
ntype mean(const nmatrix &);
nmatrix random_uniform_matrix(const unsigned int&, const unsigned int&, const ntype&);
void rescale_mean(nmatrix&, const ntype&);

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

#endif
