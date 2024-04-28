#ifndef OPTIMIZE_MATRIX_H
#define OPTIMIZE_MATRIX_H
#include "../Classes/Custom_types.h"

class MonteCarloSolver{
protected:
  MCS_Running_Parameters run_params;
  ntype(*cost_function) (const Model_parameters&, void*);

private:
  virtual bool max_steps_reached() const;
  virtual bool energy_converged(const Model_parameters&, void*);
  virtual void adjust_MCS_params();

  virtual void initialize_parameters(Model_parameters&, void*);
  virtual void choose_next_parameters(Model_parameters&, void*);
  virtual Model_parameters propose_new_parameters(const Model_parameters&, void*);

  ntype probability_density(const Model_parameters&, void*) const;

public:

  // default constructors
  MonteCarloSolver();

  virtual void optimization_procedure(Model_parameters&, void*);
  virtual void display(std::ostream&);
};

class MCSv2 : public MonteCarloSolver{
private:
  void initialize_parameters(Model_parameters&, void*);
  Model_parameters propose_new_parameters(const Model_parameters&, void*);

};

/* list of cost functions */
ntype quadratic_form(const Model_parameters &, void*);

























Model_parameters choice_next_parameters_legacy(const Model_parameters& old_params, void* params);
ntype probability_density(const Model_parameters& eco_net, const MonteCarloSolver& mcs);


void apply_MC_algorithm(Model_parameters& model_params, MonteCarloSolver& mcs, void* extra_params);
bool choose_next_parameters(Model_parameters& model_params, unsigned int fails, const MonteCarloSolver& mcs,void* params);


nmatrix create_alpha(const ntype& connectance_in, const nmatrix& gamma, bool coprophagy_allowed);
nmatrix create_gamma(unsigned int NR, unsigned int NS, const ntype& conn_targ);

nmatrix proposed_new_matrix_Alberto(const nmatrix& alpha, unsigned int steps);
nmatrix proposed_new_alpha_Leo(const nmatrix& alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps);
void modify_row(nmatrix& mat);
void modify_column(nmatrix& mat);
/* flips a random element from zero to one or the other way around. The variable coprophagy decides if coprophagy is allowed or not */
nmatrix flip_one_element(const nmatrix& alpha, const nmatrix& gamma, bool coprophagy);


ntype quadratic_form(const nmatrix& alpha, const nmatrix& gamma, void* params);
ntype quadratic_form_corrected_AlbertoMay2021(const nmatrix& alpha, const nmatrix& gamma, void* params);



#endif
