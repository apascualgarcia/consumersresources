#include "../../include/CRModel.h"


nmatrix optimal_syntrophy_from_consumption(const nmatrix& gamma, bool coprophagy, MonteCarloSolver& mcs){
  return optimal_syntrophy_from_consumption(gamma, coprophagy, mcs, connectance(gamma));
}
nmatrix optimal_syntrophy_from_consumption(const nmatrix& gamma, bool coprophagy, MonteCarloSolver& mcs, const ntype& target_conn){
  nmatrix unit_gamma(gamma.size(), nvector(gamma[0].size(), 0.));
  for(size_t i=0; i < gamma.size(); ++i){
    for(size_t mu=0; mu < gamma[0].size(); ++mu){
      if(gamma[i][mu]>0){
        unit_gamma[i][mu]=1.;
      }
    }
  }
  // the first choice of alpha by definition will have no intra specific syntrophy
  nmatrix alpha = flip_whole_binary_matrix(transpose(unit_gamma));
  EcologicalNetwork eco_net;
  eco_net.A = alpha;
  eco_net.G = unit_gamma;
  apply_MC_algorithm(eco_net, mcs);
  return alpha;
}
ntype probability_density(const nmatrix& alpha, const nmatrix& gamma, const MonteCarloSolver& mcs){
  return exp(-mcs.cost_function(alpha, gamma, mcs.additional_params)/mcs.T);
}
ntype probability_density(const EcologicalNetwork& eco_net, const MonteCarloSolver& mcs){
  return probability_density(eco_net.A, eco_net.G, mcs);
}
nmatrix proposed_new_alpha(const nmatrix & alpha, const nmatrix& gamma, bool coprophagy_allowed, unsigned int steps, const MonteCarloSolver& mcs){
  nmatrix new_alpha;
  /* As long as the leaving condition is not fulfilled, we do not leave the loop */
  bool leave_loop=false;
  do{
    switch(mcs.mcmode){
      case unconstrained: {
        /* July 2nd 2021: Leo proposes this new version: the same as Alberto but with a 0.5 probability of making a 0->1 or 1->0 in the matrix */
        new_alpha=proposed_new_alpha_Leo(alpha, gamma, coprophagy_allowed, steps);
        break;
      }
      case constant_connectance : {
        /* AS OF JULY 2ND proposed_new_alpha_Alberto VERSION WORKS, it only changes the nestedness without changing the connectance */
        new_alpha=proposed_new_matrix_Alberto(alpha, steps);
        break;
      }
      default : {
        throw error("Unknown MC mode in the proposed_new_alpha function");
        break;
      }
    }
    /* if coprophagy is allowed then we leave the loop in any case, if it is not we leave it only if
      indeed no coprophagy is observed */
    leave_loop = coprophagy_allowed or not(is_there_coprophagy(new_alpha, gamma));
  }while(not(leave_loop));
  return new_alpha;

}
EcologicalNetwork proposed_new_eco_net(const EcologicalNetwork& old_net, unsigned int steps, const MonteCarloSolver& mcs){
  EcologicalNetwork new_net;
  bool leave_loop=false;
  do{
    bool add_condition=true;
    switch(mcs.mcmode){
      case unconstrained:{
        new_net = old_net;
        flip_one_binary_matrix_element(new_net.A);
        break;
      }

      case constant_connectance:{
        new_net = old_net;
        swap_two_matrix_elements(new_net.A);
        break;
      }

      case both_modified:{
        new_net = old_net;
        swap_two_matrix_elements(new_net.G);
        flip_one_binary_matrix_element(new_net.A);
        add_condition = is_matrix_full_rank(new_net.G);
        break;
      }

      default:{
        throw error("Unknown MC mode in the proposed_new_eco_net function");
        break;
      }
    }
    leave_loop = (mcs.iss_allowed or not(is_there_coprophagy(new_net))) and add_condition;
  }while(not(leave_loop));
  return new_net;
}
/* creates an optimal consumption matrix with connectance ctarg  */
nmatrix optimal_consumption_matrix(unsigned int NR, unsigned int NS, const ntype& ctarg, MonteCarloSolver& mcs){
  nmatrix gamma = create_gamma(NR, NS, ctarg);
  nmatrix dummy(NS, nvector(NR, 0.));
  mcs.iss_allowed=true;
  EcologicalNetwork eco_net;
  eco_net.A = gamma;
  eco_net.G = dummy;
  apply_MC_algorithm(eco_net, mcs);
  return gamma;
}
nmatrix create_gamma(unsigned int NR, unsigned int NS, const ntype& ctarg){
  nmatrix gamma(NS, nvector(NR,0.));
  std::uniform_real_distribution<ntype> unif_distrib(0., 1.);
  /*  we fill gamma such that its connectance is ctarg while making sure that
      every species eats something and every resource is eaten by one species */
  unsigned int number_of_links = ctarg*NR*NS;
  unsigned int diag_elements = std::min(NR, NS);

  if(diag_elements > number_of_links){
    throw error("The target connectance is too low.");
  }

  if(diag_elements==number_of_links){
    if(has_an_empty_row(gamma)|| has_an_empty_column(gamma)){
      throw error("Please increase connectance, filling up the diagonal took all available links.");
    }
  }

  /* we fill up the diagonal elements first, in the case NR=NS this ensures that no column or row is empty */
  for(size_t k=0; k < diag_elements; ++k){
    gamma[k][k]=1;
  }

  unsigned int to_fill = number_of_links-diag_elements;
  ntype proba = ntype(to_fill)/ntype(NR*NS);

  /* finally we fill the remaining links */
  do{
    for(size_t i=0; i < NS; ++i){
      for(size_t mu=0; mu < NR ;++mu){
        if(i!=mu){
          gamma[i][mu]=0.;
          if(unif_distrib(random_engine) < proba){
            gamma[i][mu]=1.;
          }
        }
      }
    }

  }while(has_an_empty_row(gamma)|| has_an_empty_column(gamma));
  std::cout << "Created gamma with connectance " << connectance(gamma) << std::endl;
  return gamma;
}
nmatrix flip_one_element(const nmatrix& alpha, const nmatrix& gamma, bool allowed_coprophagy){
  unsigned int row_index, col_index;
  nmatrix new_alpha;

  std::uniform_int_distribution<unsigned int> row_dist(0, alpha.size()-1);
  std::uniform_int_distribution<unsigned int> col_dist(0, alpha[0].size()-1);
  do{
    /* set new_alpha to alpha so that in any case, new_alpha and alpha differ by one element */
    new_alpha=alpha;
    /* pick random element and flip it */
    row_index=row_dist(random_engine);
    col_index=col_dist(random_engine);

    new_alpha[row_index][col_index]=1-alpha[row_index][col_index];

  }while(not(allowed_coprophagy) && is_there_coprophagy(new_alpha, gamma));

  return new_alpha;
}
bool choose_next_matrix(nmatrix& alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps, unsigned int& fails, const MonteCarloSolver& mcs){
  nmatrix new_alpha=proposed_new_alpha(alpha, gamma, coprophagy, steps, mcs);
  ntype proba_ratio=probability_density(new_alpha, gamma, mcs)/probability_density(alpha, gamma, mcs);
  if(proba_ratio>1){
    alpha=new_alpha;
    fails=0;
    return true;
  }else{
    fails+=1;
    std::uniform_real_distribution<ntype> real_distrib(0., 1.);
    if(real_distrib(random_engine)<proba_ratio){
      alpha=new_alpha;
      return true;
    }
  }
  return false;
}
bool choose_next_ecological_network(EcologicalNetwork& eco_net, unsigned int steps, unsigned int fails, MonteCarloSolver& mcs){
  EcologicalNetwork new_net = proposed_new_eco_net(eco_net, steps, mcs);
  ntype proba_ratio = probability_density(new_net, mcs)/probability_density(eco_net, mcs);
  if(proba_ratio>1){
    eco_net=new_net;
    fails=0;
    return true;
  }else{
    fails+=1;
    std::uniform_real_distribution<ntype> real_distrib(0., 1.);
    if(real_distrib(random_engine)<proba_ratio){
      eco_net=new_net;
      return true;
    }
  }
  return false;
}
void apply_MC_algorithm(EcologicalNetwork& eco_net, MonteCarloSolver& mcs){
  unsigned int convergence=0, fails=0, steps=0;
  bool max_steps_reached=false, changed=true, stop=false, energy_converged=false, increase_temp=false;
  nvector last_changed_elements;
  // these vector keep track of the last required_convergence elements (same order as writing elements)
  nmatrix previous_results;
  double mean_energy=0.;
  Metaparameters* m = (Metaparameters*) mcs.additional_params;

  std::ofstream energy_file;
  if(mcs.write_mode=="all"){
    energy_file=open_external_file_truncate(mcs.energy_file);
  }

  if(mcs.write_mode=="converged_only"){
    energy_file=open_external_file_append(mcs.energy_file);
  }

  if(mcs.write_mode!="none"){
    energy_file << "# The following metaparameters were used for this matrix optimization : " << *m << std::endl;
    energy_file << "# Energy optimization to find the best syntrophy matrix for a given consumption matrix " << std::endl;
    energy_file << "# Are given, in that order: energy, A-nestedness, A-connectance, G-nestedness, G-connectance and temperature of the syntrophy matrix" << std::endl;
  }


  if(mcs.write_mode=="all"){
    energy_file << "# Each new line is a further step of the MCS optimization algorithm" << std::endl;
  }
  if(mcs.write_mode=="converged_only"){
    energy_file << "# This line is the converged value at the end of the MCS optimization algorithm" << std::endl;
  }
  const unsigned int Naverage=50, required_convergence=2000;
  const double eps=1e-3;

  while(!stop){
    double current_energy = mcs.cost_function(eco_net.A, eco_net.G, mcs.additional_params);
    double nestA, connA, nestG, connG, eff_comp;
    Metaparameters* m = (Metaparameters*)(mcs.additional_params);


    nestA = nestedness(eco_net.A);
    connA = connectance(eco_net.A);
    nestG = nestedness(eco_net.G);
    connG = connectance(eco_net.G);
    eff_comp = eco_net.effective_competition(*m);

    if(mcs.write_mode=="all"){
      energy_file << current_energy << " " << nestA << " " << connA << " " << nestG;
      energy_file << " "<< connG <<" "<< mcs.T << std::endl;
    }


    /* we check if the new matrix (computed at the previous loop or the initial one) == the old matrix. If not, count it as a "fail" */
    if(changed){
      fails=0;
      last_changed_elements.push_back(current_energy);
      if(mcs.write_mode=="converged_only"){
        previous_results.push_back(nvector{current_energy, nestA, connA, nestG, connG, mcs.T});
      }
    }else{
      fails+=1;
    }

    /* then compute if with the previous move the energy is converging */
    if(last_changed_elements.size() > Naverage){
      last_changed_elements.erase(last_changed_elements.begin());
    }

    if(mcs.write_mode=="converged_only"){
      if(previous_results.size() > required_convergence){
        previous_results.erase(previous_results.begin());
      }
    }

    mean_energy=mean(last_changed_elements);

    if(changed && last_changed_elements.size()==Naverage){
      if(abs(current_energy-mean_energy) <= eps*abs(mean_energy)){
        convergence+=1;
      }else{
        convergence=0;
      }
    }

    energy_converged = (convergence>=required_convergence);
    max_steps_reached = (steps>=mcs.max_steps);
    stop = max_steps_reached || energy_converged;
    increase_temp = fails >= mcs.max_fails;

    if(mcs.write_mode=="converged_only" && stop){
      nmatrix results = transpose(previous_results);
      for(size_t i = 0; i < results.size(); ++i){
        energy_file << mean(results[i]) << " ";
      }
      for(size_t i=0; i < results.size(); ++i){
        energy_file << standard_dev(results[i]);
        if(i==results.size()-1){
          energy_file << std::endl;
        }else{
          energy_file << " ";
        }
      }
    }

    if(steps%mcs.display_stride==0 || stop ){
      std::cout << "\t Step " << steps;
      std::cout << "\t T=" << mcs.T;
      std::cout <<"\t cost function=" << current_energy;
      std::cout << "\t nestA=" << nestA;
      std::cout << "\t connA=" << connA;
      std::cout << "\t nestG=" << nestG;
      std::cout << "\t connG=" << connG;
      std::cout << "\t eff_comp=" << eff_comp;
      if(stop){
        std::cout << std::endl << "-------" << std::endl << "STOPPING THE ALGORITHM ";
      }
      /* not displaying this message anymore since consider also negative cost functions
      if(reached_zero){
        std::cout << " -> reached zero on the cost function." << std::endl;
      }
      */

      if(max_steps_reached){
        std::cout << " -> max number of steps reached." << std::endl;
        if(mcs.write_mode=="converged_only"){
          energy_file << "# According to our algorithm, there was no convergence so please take the results with a grain of salt" << std::endl;
        }
      }

      if(energy_converged){
        std::cout << " -> cost function converging." << std::endl;
      }


     // std::cout << ", matrix = " << std::endl;
      //display_food_matrix(std::cout, alpha);
      std::cout << std::endl;
    }
    /* before getting the next matrix, we readjust the temperature if needed */
    /* when the move has not been accepted too many times, increase temp */
    if(increase_temp){
      mcs.T = mcs.T/mcs.annealing_const;
      fails=0;
    }

    /* at a given frequency, the temperature is reduced */
    if(steps%mcs.annealing_freq==0){
      mcs.T*=mcs.annealing_const;
    }

    /* Finally we get the next matrix and increase the amount of steps by 1 */
    changed=choose_next_ecological_network(eco_net, steps, fails, mcs);
    steps+=1;
  }
  energy_file.close();
  return;
}
nmatrix create_alpha(const ntype& connectance_in, const nmatrix& gamma, bool allowed_coprophagy){
  nmatrix alpha(gamma[0].size(),nvector(gamma.size(), 0));
  std::uniform_real_distribution<ntype> real_distrib(0., 1.);
  do{
    for(size_t mu=0; mu < alpha.size(); ++mu){
      for(size_t i=0; i < alpha[mu].size();++i){
        if(real_distrib(random_engine)<connectance_in){
          alpha[mu][i]=1;
        }
      }
    }
  }while(not(allowed_coprophagy) && is_there_coprophagy(alpha, gamma));
  std::cout << "Created initial alpha guess with nestedness " << nestedness(alpha) << " and connectance " << connectance(alpha) << std::endl;
  return alpha;
}
nmatrix proposed_new_matrix_Alberto(const nmatrix& alpha, unsigned int steps){
  nmatrix new_alpha=alpha;
  if(steps%2==0){
    modify_row(new_alpha);
  }else{
    modify_column(new_alpha);
  }
  return new_alpha;
}
nmatrix proposed_new_alpha_Leo(const nmatrix& alpha, const nmatrix& gamma, bool coprophagy, unsigned int steps){
  nmatrix new_alpha = alpha;
  flip_one_binary_matrix_element(new_alpha);
  return new_alpha;
}
void modify_row(nmatrix& alpha){
  /*  we first choose a random row (i.e. resource) that is neither completely empty nor
      completely filled */
  const unsigned int NR = alpha.size();
  const unsigned int NS = alpha[0].size();
  std::vector<unsigned int> sum_row(NR, 0.);


  for(size_t nu=0; nu < NR; ++nu){
    for(size_t i=0; i < NS; ++i){
      sum_row[nu]+=(unsigned int)(alpha[nu][i]);
    }
  }

  std::uniform_int_distribution<size_t> random_resources(0, NR-1);

  unsigned int mu=random_resources(random_engine);

  while(sum_row[mu]==0 || sum_row[mu]==NS){
    mu=random_resources(random_engine);
  }

  /*  then we find the indices of the columns (i.e. consumers) in that row
      (i.e. resource) which are zeros and which are ones */
  std::vector<size_t> zero_els;
  std::vector<size_t> one_els;

  for(size_t i=0; i < NS; ++i){
    if(alpha[mu][i]<1e-15){
      zero_els.push_back(i);
    }else if(alpha[mu][i]==1){
      one_els.push_back(i);
    }
  }

  /*  we choose one index which has one and another which has zero, may have a problem here :
      we did not explicitly check that zero_els.size() was at least 1 !! could have the setting where
      zero_els is empty and one_els is full or vice-versa */
  if(zero_els.size()==0){
    std::cout << "NR = " << NR << ", NS=" << NS << std::endl;
    std::cout << "Picked row " << mu << std::endl;
    std::cout << "The zero elements are " << zero_els << std::endl;
    throw error("zero_els has size 0 (row) ", 1);
  }

  if(one_els.size()==0){
    std::cout << "NR = " << NR << ", NS=" << NS << std::endl;
    std::cout << "Picked row " << mu << std::endl;
    std::cout << "The one elements are " << one_els << std::endl;
    throw error("one_els has size 0 (row) ", 1);
  }
  std::uniform_int_distribution<size_t> zero_els_indices(0, zero_els.size()-1);
  std::uniform_int_distribution<size_t> one_els_indices(0, one_els.size()-1);

  size_t zero_index=zero_els[zero_els_indices(random_engine)];
  size_t one_index=one_els[one_els_indices(random_engine)];

  if(one_index>=NS){
    display_food_matrix(std::cout, alpha);
    std::cout << "NR = " << NR << ", NS=" << NS << std::endl;
    std::cout << "one_index=" << one_index << std::endl;
    throw error("one_index has an index larger than it should have (row)",1);
  }

  if(zero_index>=NS){
    display_food_matrix(std::cout, alpha);
    std::cout << "NR = " << NR << ", NS=" << NS << std::endl;
    std::cout << "zero_index=" << zero_index << std::endl;
    throw error("zero_index has an index larger than it should have (row)",1);
  }

  /*  we check if there is another non empty value in the column (i.e.
      check if that species releases to something else) */
  bool conditionRel=false;
  for(size_t nu=0; nu < NR && !conditionRel; ++nu){
    if(mu!=nu){
      if(alpha[nu][one_index]==1){
        conditionRel=true;
      }
    }
  }

  // Remove the coprophagy part, this is taken care of in a loop above if needed
  // /* we check if with the chosen zero index there is coprophagy */
  // bool conditionCopr=false;
  // if(gamma[zero_index][mu]!=1){
  //   conditionCopr=true;
  // }else if(coprophagy){
  //   conditionCopr=true;
  // }

  /* if both conditions are fulfilled, we can swap the two elements */
  if(conditionRel){
    alpha[mu][zero_index]=1;
    alpha[mu][one_index]=0;
  }
  return;
}
void modify_column(nmatrix& alpha){
  /*  we first choose a random column (i.e consumer) that is not empty and
      not completely filled */
  const unsigned int NR = alpha.size();
  const unsigned int NS = alpha[0].size();
  std::vector<unsigned int> sum_column(NS, 0.);
  for(size_t mu=0; mu < NR; ++mu){
    for(size_t i=0; i < NS; ++i){
      sum_column[i]+=(unsigned int)(alpha[mu][i]);
    }
  }

  std::uniform_int_distribution<size_t> random_consumers(0, NS-1);
  unsigned int k=random_consumers(random_engine);
  while(sum_column[k]==0 || sum_column[k]==NR){
    k=random_consumers(random_engine);
  }

  /*  then we find the indices of the rows in that column which are zeros
      or which are one */
  std::vector<size_t> zero_els;
  std::vector<size_t> one_els;
  for(size_t mu=0; mu < NR; ++mu){
    if(alpha[mu][k]<1e-15){
      zero_els.push_back(mu);
    }else if(alpha[mu][k]==1){
      one_els.push_back(mu);
    }
  }

  if(zero_els.size()==0){
    std::ofstream myfile=open_external_file_truncate("data_output/error_matrix.out");
    std::cout << std::endl;
    display_food_matrix(myfile, alpha);
    myfile.close();
    display_food_matrix(std::cout, alpha);
    std::cout << std::endl;
    std::cout << "NR = " << NR << ", NS=" << NS << std::endl;
    std::cout << "Picked column " << k << std::endl;
    std::cout << "The zero elements are " << zero_els << std::endl;
    throw error("zero_els has size 0 (column)", 1);
  }

  if(one_els.size()==0){
    std::cout << std::endl;
    display_food_matrix(std::cout, alpha);
    std::cout << std::endl;
    std::cout << "NR = " << NR << ", NS=" << NS << std::endl;
    std::cout << "Picked column " << k << std::endl;
    std::cout << "The one elements are " << one_els << std::endl;
    throw error("one_els has size 0 (column)", 1);
  }

  /*  we then choose randomly one of the zero elements and one of the one elements */
  std::uniform_int_distribution<size_t> zero_els_indices(0, zero_els.size()-1);
  std::uniform_int_distribution<size_t> one_els_indices(0, one_els.size()-1);

  size_t zero_index=zero_els[zero_els_indices(random_engine)];
  size_t one_index=one_els[one_els_indices(random_engine)];

  if(one_index>=NR){
    std::cout << std::endl;
    display_food_matrix(std::cout, alpha);
    std::cout << "one_index=" << one_index << std::endl;
    throw error("one_index has an index larger than it should have (column)",1);
  }

  if(zero_index>=NR){
    std::cout << std::endl;
    display_food_matrix(std::cout, alpha);
    std::cout << "zero_index=" << zero_index << std::endl;
    throw error("zero_index has an index larger than it should have (column)",1);
  }

  /*  we check if there is another non empty value in the row with the one, i.e.
      if the resource is being released by another species */
  bool conditionRel=false;
  for(size_t i=0; i < NS && !conditionRel; ++i){
    if(i!=k){
      if(alpha[one_index][i]==1){
        conditionRel=true;
      }
    }
  }

  // Remove the coprophagy part, this is taken care of in a loop above if needed
  /* we check if with the chosen zero index there is coprophagy */
  // bool conditionCopr=false;
  // if(gamma[k][zero_index]!=1){
  //   conditionCopr=true;
  // }else if(coprophagy){
  //   conditionCopr=true;
  // }

  /* if both conditions are fulfilled we can swap the two elements */
  if(conditionRel){
    alpha[zero_index][k]=1;
    alpha[one_index][k]=0;
  }

  return;
}
ntype quadratic_form_low_intra_resource_interaction(const nmatrix& alpha, const nmatrix& gamma, void* params){
  Metaparameters* m = (Metaparameters*)(params);
  return m->quadratic_form_low_intra_resource_interaction(alpha, gamma);
}
ntype quadratic_form_nestedness(const nmatrix& A, const nmatrix& G, void*params){
  ntype* target = (ntype*)(params);
  return abs(nestedness(G)-(*target));
}
ntype quadratic_form_nestedness_rank(const nmatrix& A, const nmatrix& G, void*params){
  int max_rank = G.size();
  if(G[0].size() < max_rank){
    max_rank =G[0].size();
  }
  return quadratic_form_nestedness(A, G, params)+max_rank-rank(G);
}
ntype quadratic_form_LRI_with_critical_radius(const nmatrix& alpha, const nmatrix& gamma, void* params){
  Metaparameters* m= (Metaparameters*)(params);
  return m->accurate_quadratic_form_LRI(alpha, gamma);
}
ntype quadratic_form_LRI_newly_corrected(const nmatrix& alpha, const nmatrix& gamma, void* params){
  Metaparameters* m= (Metaparameters*)(params);
  return m->newly_corrected_quadratic_form_LRI(alpha, gamma);
}
ntype quadratic_form_effective_competition(const nmatrix & A, const nmatrix & G, void* params){
  EcologicalNetwork eco_net(A, G);
  Metaparameters* m = (Metaparameters*)(params);
  return eco_net.effective_competition(*m);
}

ntype quadratic_form_Alberto(const nmatrix& alpha, const nmatrix& gamma, void* params){
  const nvector& u = *(nvector*)(params);
  return u*(alpha*gamma)*u;
}

// this function is the same as "quadratic_form" except corrected according to Alberto's ideas we discussed in May 2021
ntype quadratic_form_corrected_AlbertoMay2021(const nmatrix& alpha, const nmatrix& gamma, void* params){
  unsigned int NR=alpha.size();
  nmatrix alpha_gamma=alpha*gamma;
  nmatrix gamma_square=transpose(gamma)*gamma;

  ntype to_minimize=0.;
  /* we want the absolute trace to be as close to zero as possible*/
  /* and we want the rest to be as close to zero as possible*/
  ntype off_diag=0., trace=0.;
  for(size_t mu=0; mu < NR;++mu){
    trace+=(alpha_gamma[mu][mu]-gamma_square[mu][mu]);
    for(size_t nu=0; nu < NR;++nu){
      if(nu!=mu){
        off_diag+=abs(alpha_gamma[mu][nu]-gamma_square[mu][nu]);
      }
    }
  }
  to_minimize=trace+off_diag;
  return to_minimize;
}
ntype quadratic_form(const nmatrix& A, const nmatrix& G, void* params){
  /* the goal is to minimize the maximal sum of LHS in the intra resource regime */
  Metaparameters* m= (Metaparameters*)(params);
  unsigned int NR=A.size();
  nmatrix AG=A*G;
  nmatrix GG=transpose(G)*G;
  nvector Z = nvector(NR, 0.);

  /* we want the absolute trace to be as close to zero as possible*/
  /* and we want the rest to be as close to zero as possible*/
  for(size_t mu=0; mu < NR;++mu){
    Z[mu]+=(m->alpha0*AG[mu][mu]-m->gamma0*m->R0*GG[mu][mu]);
    for(size_t nu=0; nu < NR;++nu){
      if(nu!=mu){
        Z[mu]+=abs(m->alpha0*AG[mu][nu]-m->gamma0*m->R0*GG[mu][nu]);
      }
    }
  }


  /* version with Heaviside function */
  ntype energy = 0;
  ntype X = -1.*int(NR);
  for(size_t mu = 0; mu < NR; ++mu){
    X+=Heaviside(-Z[mu]);
  }
  /* without Heaviside function */
  ntype coeff = 1;
  for(size_t mu=0; mu < NR; ++mu){
    energy+=coeff*Z[mu];
  }
  return energy;
}
