#include "../../include/CRModel.h"

void MonteCarloSolver::optimization_procedure(Model_parameters& model_params, void* extra_params){
  bool stop=false;

  initialize_parameters(model_params, extra_params);

  while(!stop){
    stop = this->max_steps_reached() || this->energy_converged(model_params, extra_params);
    this->adjust_MCS_params();
    this->choose_next_parameters(model_params, extra_params);
    this->run_params.current_step+=1;

    if(this->run_params.current_step % this->run_params.display_stride ==0){
      std::cout << *this << " Energy = " << this->cost_function(model_params, extra_params) << std::endl;
    }

    if(stop){
      std::cout << "---> convergence achieved, stopping the algorithm. " ;
      std::cout << *this << " Energy = " << this->cost_function(model_params, extra_params) << std::endl;
    }

  }

  return;
}
void MonteCarloSolver::initialize_parameters(Model_parameters& model_params, void* extra_params){
  return;
}
bool MonteCarloSolver::max_steps_reached() const{
  return (this->run_params.current_step>=this->run_params.max_steps);
}
void MonteCarloSolver::adjust_MCS_params(){
  
  /* when the move has not been accepted too many times, increase temp */
  if(this->run_params.current_fails>=this->run_params.max_fails){
    this->run_params.current_temperature /= this->run_params.annealing_const;
  }

  /* at a given frequency, the temperature is reduced */
  if(this->run_params.current_step % this->run_params.annealing_freq==0){
    this->run_params.current_temperature *= this->run_params.annealing_const;
  }
  return;
}
bool MonteCarloSolver::energy_converged(const Model_parameters& model_params, void* extra_params) {

  double current_energy = this->cost_function(model_params, extra_params);

  /* we check if the new matrix (computed at the previous loop or the initial one) == the old matrix. If it is, count it as a "fail" */
  if(this->run_params.changed_matrix){
    this->run_params.current_fails=0;
    this->run_params.last_changed_elements.push_back(current_energy);
    }else{
      this->run_params.current_fails+=1;
    }

  /* then compute if with the previous move the energy is converging */
  if(this->run_params.last_changed_elements.size() > this->run_params.convergence_Naverage){
    this->run_params.last_changed_elements.erase(this->run_params.last_changed_elements.begin());
  }

  double mean_energy=mean(this->run_params.last_changed_elements);

  if(this->run_params.changed_matrix && this->run_params.last_changed_elements.size()==this->run_params.convergence_Naverage){
    if(abs(current_energy-mean_energy) <= this->run_params.eps*abs(mean_energy)){
      this->run_params.current_convergence+=1;
    }else{
      this->run_params.current_convergence=0;
    }
  }

  return this->run_params.current_convergence>=this->run_params.required_convergence;

}
void MonteCarloSolver::display(std::ostream& os){
  this->run_params.display(os);
  return;
}
void MonteCarloSolver::choose_next_parameters(Model_parameters& current_params, void* extra_params){
  
  Model_parameters new_params = this->propose_new_parameters(current_params, extra_params);

  ntype proba_ratio = this->probability_density(new_params, extra_params)/this->probability_density(current_params, extra_params);
  
  std::uniform_real_distribution<ntype> real_distrib(0., 1.);
  
  if(real_distrib(random_engine)<proba_ratio){
      current_params=new_params;
      this->run_params.changed_matrix = true;
  }else{
      this->run_params.changed_matrix = false;
  }
  return;
}
Model_parameters MonteCarloSolver::propose_new_parameters(const Model_parameters& current_params, void*){
  Model_parameters new_params = current_params;
  flip_one_binary_matrix_element(new_params.get_parameters()->alpha);
  return new_params;
}
ntype MonteCarloSolver::probability_density(const Model_parameters& model_params, void* extra_params) const{
  return exp(-this->cost_function(model_params, extra_params)/this->run_params.current_temperature);
}
ntype quadratic_form(const Model_parameters& model_params, void* extra_params){
  /* the goal is to minimize the maximal sum of LHS in the intra resource regime */
  Metaparameters* m= (Metaparameters*)(extra_params);

  nmatrix alpha = (model_params.get_parameter_set().alpha);
  nmatrix gamma = (model_params.get_parameter_set().gamma);

  unsigned int NR=alpha.size();
  nmatrix AG=alpha*gamma;
  nmatrix GG=transpose(gamma)*gamma;
  nvector Z = nvector(NR, 0.);

  /* we want the absolute trace to be as close to zero as possible*/
  /* and we want the rest to be as close to zero as possible*/
  for(size_t mu=0; mu < NR;++mu){
    Z[mu]+=(AG[mu][mu]-m->R0*GG[mu][mu]);
    for(size_t nu=0; nu < NR;++nu){
      if(nu!=mu){
        Z[mu]+=abs(AG[mu][nu]-m->R0*GG[mu][nu]);
      }
    }
  }


  ntype energy = 0;

  /* version with Heaviside function */
  /* ntype X = -1.*int(NR);
  for(size_t mu = 0; mu < NR; ++mu){
    X+=Heaviside(-Z[mu]);
  } */
  
  /* without Heaviside function */
  ntype coeff = 1;
  for(size_t mu=0; mu < NR; ++mu){
    energy+=coeff*Z[mu];
  }
  return energy;
}
MonteCarloSolver::MonteCarloSolver(){
  this->cost_function = quadratic_form;
}
void MCSv2::initialize_parameters(Model_parameters& mod_params, void* extra_params){
  
  Metaparameters* m = (Metaparameters*)(extra_params);
  ntype max_alpha0 = m->gamma0*m->R0*m->NR;
  if(m->sigma0 < 0.5){
    max_alpha0 *= m->sigma0;
  }else{
    max_alpha0 *= (1-m->sigma0);
  }

  std::uniform_real_distribution<ntype> real_distrib(0., max_alpha0);

  /* create a random binary matrix*/
  nmatrix new_alpha = create_random_binary_matrix(m->NR,m->NS);
  for(size_t mu = 0; mu < new_alpha.size(); ++mu){
    for(size_t i=0; i < new_alpha[mu].size();++i){
      new_alpha[mu][i] *= real_distrib(random_engine); 
    }
  }

  mod_params.get_parameters()->alpha = new_alpha;

  return;
}
Model_parameters MCSv2::propose_new_parameters(const Model_parameters& current_params, void* extra_params){
  Model_parameters new_params = current_params;
  
  /* ===================

  pick a random element : 
  if it is not zero, draw an alpha0 value, otherwise make it zero
  
  =================== */

  Metaparameters* m = (Metaparameters*)(extra_params);
  ntype max_alpha0 = m->gamma0*m->R0*m->NR;
  if(m->sigma0 < 0.5){
    max_alpha0 *= m->sigma0;
  }else{
    max_alpha0 *= (1-m->sigma0);
  }

  std::uniform_real_distribution<ntype> real_distrib(0., max_alpha0);
  std::uniform_int_distribution<unsigned int> distrib_NR(0, m->NR-1);
  std::uniform_int_distribution<unsigned int> distrib_NS(0, m->NS-1);

  /* pick random element */
  unsigned int mu = distrib_NR(random_engine), i=distrib_NS(random_engine);
  if(new_params.get_parameters()->alpha[mu][i]!=0){
    new_params.get_parameters()->alpha[mu][i]=0;
  }else{
    new_params.get_parameters()->alpha[mu][i] = real_distrib(random_engine);
  }

  return new_params;
}