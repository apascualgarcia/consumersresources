#ifndef OPTIMIZE_ALPHA_MATRIX_H
#define OPTIMIZE_ALPHA_MATRIX_H
#include "../Classes/Custom_types.h"
nmatrix optimal_syntrophy_from_consumption(const nmatrix& gamma, bool coprophagy, const MonteCarloSolver& mcs);
ntype quadratic_form(const nmatrix& alpha, const nmatrix& gamma, const nvector& v);
ntype probability_density(const nmatrix& alpha, const nmatrix& gamma, const nvector& u, const ntype& T);
nmatrix proposed_new_alpha(const nmatrix & alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps);
void choose_next_alpha(nmatrix& alpha, const nmatrix& gamma, const nvector& u, bool coprophagy, const ntype& T, unsigned int steps, unsigned int& fails);
void apply_MC_algorithm(nmatrix& alpha, const nmatrix& gamma, const nvector& u, bool coprophagy, const MonteCarloSolver& mcs);
nmatrix create_alpha(const ntype& connectance_in, const nmatrix& gamma);
nmatrix proposed_new_alpha_Alberto(const nmatrix& alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps);
void modify_row(nmatrix& alpha, const nmatrix& gamma, bool coprophagy);
void modify_column(nmatrix& alpha, const nmatrix& gamma, bool coprophagy);

#endif
