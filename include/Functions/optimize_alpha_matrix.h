#ifndef OPTIMIZE_ALPHA_MATRIX_H
#define OPTIMIZE_ALPHA_MATRIX_H
#include "../Classes/Custom_types.h"
nmatrix optimal_syntrophy_from_consumption(const nmatrix& gamma, bool coprophagy, MonteCarloSolver& mcs);

ntype probability_density(const nmatrix& alpha, const nmatrix& gamma, const MonteCarloSolver& mcs);
nmatrix proposed_new_alpha(const nmatrix & alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps);
bool choose_next_alpha(nmatrix& alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps, unsigned int& fails, const MonteCarloSolver& mcs);
void apply_MC_algorithm(nmatrix& alpha, const nmatrix& gamma, bool coprophagy, MonteCarloSolver& mcs);
nmatrix create_alpha(const ntype& connectance_in, const nmatrix& gamma);
nmatrix proposed_new_alpha_Alberto(const nmatrix& alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps);
void modify_row(nmatrix& alpha, const nmatrix& gamma, bool coprophagy);
void modify_column(nmatrix& alpha, const nmatrix& gamma, bool coprophagy);


ntype quadratic_form(const nmatrix& alpha, const nmatrix& gamma, void* params);
ntype quadratic_form_Alberto(const nmatrix& alpha, const nmatrix& gamma, void* params);
ntype quadratic_form_low_intra_resource_interaction(const nmatrix&, const nmatrix&, void*);
#endif
