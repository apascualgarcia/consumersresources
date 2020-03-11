#ifndef OPTIMIZE_ALPHA_MATRIX_H
#define OPTIMIZE_ALPHA_MATRIX_H
#include "../Classes/Custom_types.h"

// CAREFUL, HERE COPROPHAGY variable IS TO BE TAKEN IN THE SENSE (coprophagy is allowed or not)
// if coprophagy, we can observe it, otherwise we do not

nmatrix optimal_syntrophy_from_consumption(const nmatrix& gamma, bool coprophagy, MonteCarloSolver& mcs);
ntype probability_density(const nmatrix& alpha, const nmatrix& gamma, const MonteCarloSolver& mcs);
nmatrix proposed_new_alpha(const nmatrix & alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps);
bool choose_next_alpha(nmatrix& alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps, unsigned int& fails, const MonteCarloSolver& mcs);
void apply_MC_algorithm(nmatrix& alpha, const nmatrix& gamma, bool coprophagy, MonteCarloSolver& mcs);
nmatrix create_alpha(const ntype& connectance_in, const nmatrix& gamma, bool coprophagy_allowed);
nmatrix proposed_new_alpha_Alberto(const nmatrix& alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps);
void modify_row(nmatrix& alpha, const nmatrix& gamma, bool coprophagy);
void modify_column(nmatrix& alpha, const nmatrix& gamma, bool coprophagy);
/* flips a random element from zero to one or the other way around. The variable coprophagy decides if coprophagy is allowed or not */
nmatrix flip_one_element(const nmatrix& alpha, const nmatrix& gamma, bool coprophagy);


ntype quadratic_form(const nmatrix& alpha, const nmatrix& gamma, void* params);
ntype quadratic_form_Alberto(const nmatrix& alpha, const nmatrix& gamma, void* params);
ntype quadratic_form_low_intra_resource_interaction(const nmatrix&, const nmatrix&, void*);
#endif
