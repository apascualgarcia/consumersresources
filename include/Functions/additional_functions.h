#ifndef ADDITIONAL_FUNCTIONS_H
#define ADDITIONAL_FUNCTIONS_H

#include "../Classes/Metaparameters.h"
#include "../Classes/Custom_types.h"
#include "../Classes/Extinction.h"
#include<random>

foodmatrix load_food_matrix(const Metaparameters&);
ntype norm(const nvector&);
ntype mean(const nmatrix &);
nmatrix random_uniform_matrix(const unsigned int&, const unsigned int&, const ntype&);
void rescale_mean(nmatrix&, const ntype&);

bool non_neg_elements(const nmatrix&);
bool non_neg_elements(const nvector&);

void initialize_random_engine(const Metaparameters&);
void print_rand_number();

ntype mean(const nvector&);
ntype standard_dev(const nvector&);

nvector operator+(const nvector&, const nvector&);
nvector operator-(const nvector&);
nvector operator-(const nvector&, const nvector&);

#endif
