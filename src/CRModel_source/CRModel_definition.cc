#include "CRModel.h"
#include "ConfigFile.tpp"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>
#include <string>
#include <iomanip>
#include <ctime>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>


//std::random_device rdevice;
std::mt19937 random_engine;

/*** ALL THE CONSTRUCTORS for the classes ****/
Metaparameters::Metaparameters(int argc, char *argv[]){
  std::string inputPath("configuration.in"); // Fichier d'input par defaut
  if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice3
                // config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans
                                    // une "map" de strings.

  for (int i(2); i < argc; ++i) // Input complementaires ("./Onde
                                // config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  this->gamma0 = configFile.get<ntype>("gamma_0");
  this->alpha0 = configFile.get<ntype>("alpha_0");
  this->sigma0 = configFile.get<ntype>("sigma_0");
  this->p = configFile.get<ntype>("p");
  this->R0 = configFile.get<ntype>("R0");
  this->S0 = configFile.get<ntype>("S0");
  this->NR = configFile.get<unsigned int>("NR");
  this->NS = configFile.get<unsigned int>("NS");
  this->gamma_mode = string_to_gamma_mode(configFile.get<std::string>("gamma_mode"));
  this->tau_mode = string_to_tau_mode(configFile.get<std::string>("tau_mode"));
  this->alpha_mode = string_to_alpha_mode(configFile.get<std::string>("alpha_mode"));
  this->foodmatrixpath = configFile.get<std::string>("path_to_food_matrix");
  this->verbose=configFile.get<unsigned int>("verbose-level");
  this->energy_constraint=configFile.get<bool>("energy_constraint");
  this->nb_attempts = configFile.get<unsigned int>("number_of_attempts");
  this->save_path = configFile.get<std::string>("path_to_save_file");
  this->epsilon = configFile.get<ntype>("epsilon");
  this->l0 = configFile.get<ntype>("l_0");
  this->budget_constraint = configFile.get<bool>("budget_constraint");
  this->seed_number = configFile.get<unsigned int>("seed_number");
  this->tf = configFile.get<ntype>("tf");
  this->perturb_eq = configFile.get<ntype>("perturbation_equilibrium");
  this->perturb_parameters = configFile.get<ntype>("perturbation_parameters");
  if(this->verbose > 0){
    std::cout << "Loading metaparameters from " << inputPath << std::endl;
  }
  return;
}
Dynamical_variables::Dynamical_variables(nvector* resources_, nvector* consumers_){
  resources=resources_;
  consumers=consumers_;
  return;
}
Model_parameters::Model_parameters(){
  params = new Parameter_set();
  return;
}
CRModel::CRModel(Model_parameters* mod_params){
  model_param=mod_params;
  eq_vals=NULL;
  metaparameters=NULL;
  return;
}
CRModel::CRModel(const foodmatrix& F, Metaparameters& meta){
  unsigned int attempts(0);
  this->model_param = new Model_parameters();
  Parameter_set* p = model_param->get_parameters();

  do{
    attempts+=1;

    nvector Req = build_resources(meta);
    nvector Seq = build_consumers(meta);

    this->eq_vals = new ntensor();
    nmatrix equilibria;
    equilibria.push_back(Req);
    equilibria.push_back(Seq);
    eq_vals->push_back(equilibria);

    p->NR = meta.NR;
    p->NS = meta.NS;

    // first sigma, Req, Seq are drawn randomly
    p->sigma=build_sigma(meta);

    // then we build gamma according to the food matrix
    p->gamma=build_gamma(F,meta);

    // then we build alpha according to the other parameters
    p->alpha=build_alpha(p, meta, Req, attempts);
    p->tau = build_tau(p, meta, attempts);

    // d is then set
    nvector d;
    for (size_t i=0; i < meta.NS; ++i){
      ntype result = 0.;
      for (size_t mu =0 ; mu < meta.NR; ++mu){
        result+=(p->sigma)[i][mu]*(p->gamma)[i][mu]*Req[mu]-(p->tau)[mu][i];
      }
      d.push_back(result);
    }
    p->d = d;

    // still have to set l and m
    nvector l, m;
    l = build_l(meta);
    for(size_t nu=0; nu < p->NR; ++nu){
      ntype C = 0.;
      for(size_t j = 0; j < p->NS; ++j){
        C+=(p->alpha[nu][j]*Seq[j]-p->gamma[j][nu]*Req[nu]*Seq[j]);
      }
      m.push_back(ntype(l[nu]+C)/Req[nu]);
    }

    // std::exponential_distribution<ntype> exp_distrib(1.);
    // for(size_t nu = 0; nu < p->NR; ++nu){
    //   ntype  C = 0.;
    //   for(size_t j = 0; j < p->NS; ++j){
    //     C+= (p->alpha[nu][j]-p->gamma[j][nu]*Req[nu])*Seq[j];
    //   }
    //   if(C > 0.){
    //     l.push_back(exp_distrib(random_device));
    //   }else{
    //     l.push_back(exp_distrib(random_device) - C);
    //   }
    //   m.push_back(ntype(l[nu]+C)/Req[nu]);
    // }

    p->l = l;
    p->m = m;
  }while(not(this->constraints_fulfilled(meta)));

  this->metaparameters = &meta;
  if(meta.verbose > 1){
    std::cout << "Feasible system build in "<<attempts<<" iteration(s). ";
    if(attempts > meta.nb_attempts){
      std::cout << "The metaparameters had to be changed.";
    }else{
      std::cout << "It was possible to use the initial metaparameters.";
    }
    std::cout << std::endl;
  }


  return;
}


/*** ALL THE DESTRUCTORS for the classes ****/
Model_parameters::~Model_parameters(){
  delete this->params;
  return;
}

/*** ALL THE FUNCTIONS RELATED TO THE ALGORITHMIC BUILDING OF THE PARAMETER SET ***/
nmatrix build_sigma(const Metaparameters& m){
  nmatrix sigma(m.NS, nvector(m.NR,0.));
  std::uniform_real_distribution<ntype> sigma_distribution((1.-m.epsilon)*m.sigma0, (1.+m.epsilon)*m.sigma0);
  for(size_t i=0; i < sigma.size(); ++i){
    for(size_t mu=0; mu < sigma[i].size(); ++mu){
      sigma[i][mu] = sigma_distribution(random_engine);
    }
  }
  //rescale_mean(sigma, m.sigma0);
  return sigma;
}
nvector build_resources(const Metaparameters& m){
  nvector v(m.NR, 0.);
  /* OLD CODE for old distribution
  std::uniform_real_distribution<ntype> distribution(0., 1.);
  ntype avg = 0.;
  for(size_t i=0; i < v.size(); ++i){
    v[i] = distribution(random_device);
    avg += v[i]/m.NR;
  }
  for(size_t i=0; i < v.size(); ++i){
    v[i] *= (m.R0/avg);
  }*/

  std::uniform_real_distribution<ntype> R_distrib((1.-m.epsilon)*m.R0, (1.+m.epsilon)*m.R0);
  for(size_t i=0; i < v.size(); ++i){
    v[i] = R_distrib(random_engine);
  }
  return v;
}
nvector build_consumers(const Metaparameters& m){
  nvector v(m.NS, 0.);
  /* OLD CODE for old distribution
  std::uniform_real_distribution<ntype> distribution(0., 1.);
  ntype avg = 0.;
  for(size_t i=0; i < v.size(); ++i){
    v[i] = distribution(random_device);
    avg += v[i]/m.NS;
  }
  for(size_t i=0; i < v.size(); ++i){
    v[i] *= (m.S0/avg);
  }*/
  std::uniform_real_distribution<ntype> S_distrib((1.-m.epsilon)*m.S0, (1.+m.epsilon)*m.S0);
  for(size_t i=0; i < v.size(); ++i){
    v[i] = S_distrib(random_engine);
  }
  return v;
}
nmatrix build_gamma(const foodmatrix& F, const Metaparameters& m){
  nmatrix gamma(F.size(), nvector(F[0].size(), 0.));

  if(F.size()!=m.NS or F[0].size()!=m.NR){
    std::cerr << "Meta parameters not compatible with food matrix, please change this " << std::endl;
  }

  /* NEW CODE (implements gamma0 as interaction strength and not as budget)*/
  if(m.gamma_mode==random_val){
    std::uniform_real_distribution<ntype> gamma_distrib((1.-m.epsilon)*m.gamma0, (1.+m.epsilon)*m.gamma0);
    for(size_t i=0; i < gamma.size(); ++i){
      bool species_existing = false;
      for(size_t mu=0; mu < gamma[i].size(); ++mu){
        if(F[i][mu] > 0.){
          species_existing = true;
          gamma[i][mu] = F[i][mu] * gamma_distrib(random_engine);
        }
      }
      if(not(species_existing)){
        std::cerr << "Problem in the food matrix, species " << i << " does not eat anything." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  /* OLD CODE (implements gamma0 as budget and not as interaction strength)
  // build the gamma matrix
  if(m.gamma_mode==random_val){
    std::normal_distribution<ntype> distribution(1., 0.1);
    for(size_t i=0; i < gamma.size(); ++i){
      for(size_t j=0; j < gamma[i].size(); ++j){
        if(pow(F[i][j], 2.)  > 0.){
          do{
            gamma[i][j] = distribution(random_device);
          } while(gamma[i][j] < 0.);
        }
      }
    }
  }else{
    std::cerr << "different behaviours for different gamma mode not implemented yet, gamma will be returned as 0." << std::endl;
    return gamma;
  }*/

  // rescale the gamma matrix to satisfy budget constraint
  if(m.budget_constraint){
    for(size_t i=0; i < gamma.size(); ++i){
      ntype eating_cap=0.;
      bool species_existing = false;
      for(size_t j = 0; j < gamma[i].size(); ++j){
        if(gamma[i][j] > 0.){
          species_existing = true;
          eating_cap += gamma[i][j];
        }
      }
      if(species_existing){
        for(size_t j=0; j < gamma[i].size(); ++j){
          gamma[i][j] /= eating_cap;
          gamma[i][j] *= (m.gamma0*m.NR);
        }
      }else{
        std::cerr << "Problem in the food matrix, species " << i << " does not eat anything." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  return gamma;
}
nmatrix build_alpha(const Parameter_set* p, Metaparameters& m, const nvector& Req, unsigned int attempts){
  std::uniform_real_distribution<ntype> empty_or_not_distrib(0., 1.);
  std::uniform_real_distribution<ntype> alpha_distrib((1-m.epsilon)*m.alpha0, (1+m.epsilon)*m.alpha0);
  nmatrix alpha = nmatrix(p->NR, nvector(p->NS, 0.));
  if(attempts <= m.nb_attempts){
    switch(m.alpha_mode){
      case random_structure:{
        for(size_t mu=0; mu < p->NR; ++mu){
          for(size_t i=0; i < p->NS; ++i){
            alpha[mu][i] = alpha_distrib(random_engine);
          }
        }
        break;
      }

      case no_release_when_eat:{
        for(size_t i=0; i < p->NS; ++i){
          for(size_t mu=0; mu < p->NR; ++mu){
            if(p->gamma[i][mu]>0. and ntype(empty_or_not_distrib(random_engine)) < m.p){
              alpha[mu][i] = alpha_distrib(random_engine);
            }
          }
        }
        //rescale_mean(alpha, m.alpha0);
      }
      default:{
        std::cerr << " This case of alpha mode has not been implemented yet " << std::endl;
      }
    }
  }else{
    if(m.tau_mode==tau0){
      if(m.energy_constraint){
        for(size_t i = 0; i < p->NS; ++i){
          nvector v;
          for(size_t nu = 0; nu < p->NR; ++nu){
            v.push_back(p->sigma[i][nu]*Req[nu]);
          }
          ntype max_value = *std::min_element(v.begin(), v.end());
          std::uniform_real_distribution<ntype> distrib(0., max_value) ;
          for(size_t nu = 0; nu < p->NR; ++nu){
            alpha[nu][i] = distrib(random_engine);
          }
        }

      }
    }else if(m.tau_mode==taualpha){
      if(m.energy_constraint){
        for(size_t i = 0; i < p->NS; ++i){
          nvector v, w;
          for(size_t nu = 0; nu < p->NR; ++nu){
            v.push_back(p->sigma[i][nu]*Req[nu]);
            w.push_back((1-p->sigma[i][nu])*Req[nu]);
          }
          ntype max_value = *std::min_element(v.begin(), v.end());
          if (*std::min_element(w.begin(), w.end())<max_value){
            max_value=*std::min_element(w.begin(), w.end());
          }
          std::uniform_real_distribution<ntype> distrib(0., max_value);
          for(size_t nu = 0; nu < p->NR; ++nu){
            alpha[nu][i] = distrib(random_engine);
          }
        }
    }
  }
    m.alpha0 = mean(alpha);
  }
  /* OLD CODE for when the Metaparameters might not be consistent
  nvector alpha_upperbound;

  // we compute the upperbound of each column sum of alpha
  switch(m.tau_mode){
    case tau0:{
      for(size_t i = 0; i < p->NS; ++i){
        ntype mass_upperbound = 0.;
        for(size_t nu = 0; nu < p->NR; ++nu){
          mass_upperbound+=(1.-p->sigma[i][nu])*p->gamma[i][nu]*Req[nu];
        }
        alpha_upperbound.push_back(mass_upperbound);
      }
      break;
    }
    case taualpha:{
      for(size_t i = 0; i < p->NS; ++i){
        ntype mass_upperbound(0.), feasability_upperbound(0.);
        for(size_t nu = 0; nu < p->NR; ++nu){
          mass_upperbound+=(1.-p->sigma[i][nu])*p->gamma[i][nu]*Req[nu];
          feasability_upperbound+=p->sigma[i][nu]*p->gamma[i][nu]*Req[nu];
        }
        alpha_upperbound.push_back(std::min(mass_upperbound, feasability_upperbound));
      }
      break;
    }
    default:{
      std::cerr << "build_alpha for this value of tau_mode not implemented yet";
      break;
    };
  };
  */
  return alpha;
}
nmatrix build_tau(Parameter_set* p, Metaparameters& m, unsigned int attempts){
  nmatrix tau;
  switch(m.tau_mode){
    case tau0:{
      tau = nmatrix(p->NR, nvector(p->NS, 0.));
      break;
    }
    case taualpha:{
      tau = p->alpha;
      break;
    }
    default: {
      std::cerr << "tau_mode value not implemented yet";
      break;
    }
  };
  return tau;
}
nvector build_l(const Metaparameters & m){
  std::uniform_real_distribution<ntype> l_distrib((1.-m.epsilon)*m.l0, (1.+m.epsilon)*m.l0);
  nvector v;
  for(size_t nu=0; nu < m.NR; ++nu){
    v.push_back(l_distrib(random_engine));
  }
  return v;
}

/*** ALL THE ADDITIONAL FUNCTIONS related to CRModel ****/
nvector CRModel::equations_of_evolution(const Dynamical_variables& dyn_var) const{
  nvector v;
  const Parameter_set* p = this->model_param->get_parameters();
  const nvector R = *dyn_var.get_resources();
  const nvector S = *dyn_var.get_consumers();

  for (size_t nu=0; nu < p->NR; ++nu){
    ntype result(0.);
    result+=p->l[nu];
    result-=p->m[nu]*R[nu];
    for (size_t j=0; j < p->NS; ++j){
      result-=p->gamma[j][nu]*R[nu]*S[j];
      result+=p->alpha[nu][j]*S[j];
    }
    v.push_back(result);
  }

  for (size_t i=0; i<p->NS; ++i){
    ntype result=0.;
    for (size_t mu=0; mu < p->NR; ++mu){
      result+=p->sigma[i][mu]*p->gamma[i][mu]*S[i]*R[mu];
      result-=p->tau[mu][i]*S[i];
    }
    result-=p->d[i]*S[i];
    v.push_back(result);
  }
  return v;
}
nmatrix CRModel::jacobian(const Dynamical_variables& dyn_var) const{
  const Parameter_set* p = this->model_param->get_parameters();
  const nvector R = *dyn_var.get_resources();
  const nvector S = *dyn_var.get_consumers();

  nmatrix J(p->NR+p->NS, nvector(p->NR+p->NS, 0.));
  // fills the dR/dR part (should always be diagonal since a resource never directly depends on another)
  for(size_t nu=0; nu < p->NR; ++nu){
    J[nu][nu] = -p->m[nu];
    for(size_t j=0; j < p->NS; ++j){
      J[nu][nu] -= p->gamma[j][nu]*S[j];
    }
  }

  // fills the dR/dS part (upper right part of jacobian)
  for(size_t nu=0; nu < p->NR;++nu){
    for(size_t i=0; i < p->NS; ++i){
      J[nu][i+p->NR] = p->alpha[nu][i]-p->gamma[i][nu]*R[nu];
    }
  }

  // fills the dS/dR part (lower left)
  for(size_t i=0; i < p->NS; ++i){
    for(size_t nu=0; nu < p->NR; ++nu){
      J[p->NR+i][nu] = p->sigma[i][nu]*p->gamma[i][nu]*S[i];
    }
  }

  // finally, fills the dR/dR part (lower right). It is diagonal for the same reason
  for (size_t i=0; i < p->NS; ++i){
    for(size_t mu=0; mu < p->NR; ++mu){
      J[i+p->NR][i+p->NR] += p->sigma[i][mu]*p->gamma[i][mu]*R[mu]-p->tau[mu][i];
    }
    J[i+p->NR][i+p->NR] -= p->d[i];
  }

  return J;
}
nmatrix CRModel::jacobian_at_equilibrium() const {
  nvector Req = (*eq_vals)[0][0];
  nvector Seq = (*eq_vals)[0][1];

  Dynamical_variables dyn_var(&Req, &Seq);
  nmatrix jac_eq = jacobian(dyn_var);
  return jac_eq;
}
ncvector CRModel::eigenvalues_at_equilibrium() const{
  ncvector v;
  nmatrix jac_eq = this->jacobian_at_equilibrium();
  ntype min_element(std::abs(jac_eq[0][0]));

  for(size_t i = 0; i < jac_eq.size(); ++i){
    for(size_t j=0; j < jac_eq[i].size(); ++j){
      if(jac_eq[i][j]*jac_eq[i][j] > 0. and std::abs(jac_eq[i][j]) < min_element){
        min_element = std::abs(jac_eq[i][j]);
      }
    }
  }

  // for testing purpose
  //min_element = 1.;
  //std::cout << " put min_element as 1" << std::endl;

  const unsigned int NR = this->model_param->get_parameters()->NR;
  const unsigned int NS = this->model_param->get_parameters()->NS;

  Eigen::Matrix<ntype, Eigen::Dynamic, Eigen::Dynamic> jacob;
  jacob.resize(NR+NS, NR+NS);
  // we rescale the jacobian such that even the smallest value is of order 1
  for(size_t i=0; i < NR+NS; ++i){
    for(size_t j=0; j < NR+NS; ++j){
      jacob(i,j) = jac_eq[i][j]/min_element;
    }
  }
  Eigen::Matrix<nctype, Eigen::Dynamic, 1> eivals = jacob.eigenvalues();
  for(size_t i=0; i < eivals.rows(); ++i){
    v.push_back(eivals(i)*min_element);
    //v.push_back(eivals(i));
  }
  return v;
}
void CRModel::save(std::ostream& os) const{
  std::cout << " Still have to implement this, bye" << std::endl;
  return;
}
std::ostream& CRModel::display(std::ostream& os) const{
  os << "The model is defined by the following parameters :" << std::endl;
  os << *model_param << std::endl;
  os << "The following equilibrium values have been found :" << std::endl;
  os << *eq_vals << std::endl;
  os << "The model is characterised by the following metaparameters" << std::endl;
  os << *metaparameters << std::endl;
  return os;
}
bool CRModel::constraints_fulfilled(const Metaparameters& m) const{
  if(not(this->positive_parameters())){
    return false;
  }
  if(m.energy_constraint and not(this->energy_constraint())){
    return false;
  }
  return true;
}
bool CRModel::energy_constraint() const{
  Parameter_set* p = this->model_param->get_parameters();
  for(size_t j = 0 ; j < this->eq_vals->size();++j){
    nvector Req = (*eq_vals)[j][0];
    for(size_t i = 0; i < p->NS; ++i){
      ntype somme = 0.;
      for(size_t nu=0; nu < p->NR; ++nu){
        somme+=(1-p->sigma[i][nu])*p->gamma[i][nu]*Req[nu];
        somme-=p->alpha[nu][i];
      }
      if(somme < 0.){
        return false;
      }
    }
  }
  return true;
}
bool CRModel::positive_parameters() const{
  Parameter_set* p = this->model_param->get_parameters();
  if(not(non_neg_elements(p->d))){
    std::cerr << "d contains negative elements" << std::endl;
    return false;
  }
  if(not(non_neg_elements(p->l))){
    std::cerr << "l contains negative elements" << std::endl;
    return false;
  }
  if(not(non_neg_elements(p->m))){
    std::cerr << "m contains negative elements : ";
    std::cerr << p->m << std::endl;
    return false;
  }
  return true;
}
bool CRModel::dynamically_stable() const{
  ncvector v = this->eigenvalues_at_equilibrium();
  for(size_t i = 0; i< v.size(); ++i){
    if(norm(v[i])>EIGENSOLVER_PRECISION){
      return false;
    }
  }
  return true;
}
void CRModel::save_simulation() const{
  std::ofstream myfile;
  Metaparameters* m = this->metaparameters;
  myfile.open(m->save_path,std::ios_base::app);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open "<<m->save_path <<" for saving the simulation "<< std::endl;
  }else{
    save_success=true;
    time_t now=time(0);
    /*
    tm* ltm= localtime(&now);
    myfile << 1900+ltm->tm_year << " " << 1+ltm->tm_mon << " " << ltm->tm_mday<<" ";
    myfile << 1+ltm->tm_hour << ":"<<1+ltm->tm_min << ":"<<1+ltm->tm_sec<<" ";
    */
    myfile << "# "  << m->p << " " << m->epsilon << " " << m->foodmatrixpath << " " << m->verbose << " ";
    myfile << m->energy_constraint << " " << m->budget_constraint << " " << m->nb_attempts<<" " << m->seed_number <<" ";
    myfile << m->save_path << " ";
    myfile << ctime(&now) ;
    myfile << m->gamma0 <<" " <<m->alpha0 << " " << m->sigma0 << " ";
    myfile << m->p << " " << m->R0 << " " << m->S0 << " " <<m->l0<<" "<< m->NR << " ";
    myfile << m->NS;
    ncvector eigenvalues = this->eigenvalues_at_equilibrium();
    for(size_t i=0; i<eigenvalues.size(); ++i){
      myfile << " " << eigenvalues[i];
    }
    myfile << std::endl;
  }
  if(m->verbose > 0){
    if(save_success){
      std::cout << "Successfully saved model to " << m->save_path << std::endl;
    }
  }
  myfile.close();
  return;
}
void CRModel::save_jacobian_at_equilibrium(std::string savepath) const{
  Metaparameters* m = this->metaparameters;
  std::ofstream myfile, myfile2;
  myfile.open(savepath,std::ios_base::app);
  myfile2.open(savepath+".meta", std::ios_base::app);
  bool save_success(false);
  if(not(myfile.is_open()) or not(myfile2.is_open())){
    std::cerr << "Could not open "<<savepath <<" for saving the simulation "<< std::endl;
  }else{
    save_success=true;
    time_t now=time(0);
    /*
    tm* ltm= localtime(&now);
    myfile << 1900+ltm->tm_year << " " << 1+ltm->tm_mon << " " << ltm->tm_mday<<" ";
    myfile << 1+ltm->tm_hour << ":"<<1+ltm->tm_min << ":"<<1+ltm->tm_sec<<" ";
    */
    myfile2 << "# "  << m->p << " " << m->epsilon << " " << m->foodmatrixpath << " " << m->verbose << " ";
    myfile2 << m->energy_constraint << " " << m->budget_constraint << " " << m->nb_attempts<<" " << m->seed_number <<" ";
    myfile2 << savepath << " ";
    myfile2 << ctime(&now) ;
    myfile2 << m->gamma0 <<" " <<m->alpha0 << " " << m->sigma0 << " ";
    myfile2 << m->p << " " << m->R0 << " " << m->S0 << " " <<m->l0<<" "<< m->NR << " ";
    myfile2 << m->NS << std::endl;

    myfile << this->jacobian_at_equilibrium() << "#"<< std::endl;

  }
  if(m->verbose > 0){
    if(save_success){
      std::cout << "Successfully saved model to " << savepath << std::endl;
    }
  }
  myfile.close();
  return;
}
void CRModel::write_death_rates(std::string savepath) const{
  std::ofstream myfile, myfile2, myfile3;
  Parameter_set* p = this->model_param->get_parameters();
  myfile.open(savepath+".resources", std::ios_base::app);
  if(not(myfile.is_open())){
    std::cerr << "Could not open "<<savepath <<" for saving the resources death rates "<< std::endl;
  }else{
    for(size_t nu = 0; nu < p->NR; ++nu){
      myfile << p->m[nu] << " ";
    }
    myfile<<std::endl;
  }

  myfile2.open(savepath+".consumers", std::ios_base::app);
  if(not(myfile2.is_open())){
    std::cerr << "Could not open "<<savepath <<" for saving the consumers death rates "<< std::endl;
  }else{
    for(size_t i = 0; i < p->NS; ++i){
      myfile2 << p->d[i] << " ";
    }
    myfile2 << std::endl;
  }

  Metaparameters* m = this->metaparameters;
  myfile3.open(savepath+".meta", std::ios_base::app);
  if(not(myfile3.is_open())){
    std::cerr << "Could not open "<<savepath <<" for saving the metaparameters of death rates "<< std::endl;
  }else{
    myfile3 << m->gamma0 <<" " <<m->alpha0 << " " << m->sigma0 << " ";
    myfile3 << m->p << " " << m->R0 << " " << m->S0 << " " <<m->l0<<" "<< m->NR << " ";
    myfile3 << m->NS << std::endl;
  }

  myfile.close();
  myfile2.close();
  myfile3.close();

  return;
}
void CRModel::write_time_evolution(const Dynamical_variables & dyn, ntype tf) const{
  Metaparameters* m = this->metaparameters;
  std::ofstream myfile;
  myfile.open(m->save_path, std::ios::app);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << m->save_path << " to write the temporal evolution" << std::endl;
  }else{
    save_success = true;
  }

  nmatrix evolution = this->time_evolution(dyn, tf);
  for(size_t i =0; i < evolution.size() ; ++i){
    myfile << evolution[i] << std::endl;
  }

  if(m->verbose > 0){
    if(save_success){
      std::cout << "Successfully saved temporal evolution of system to " << m->save_path << std::endl;
    }
  }
  myfile.close();
  return;
}
void CRModel::write_time_evolution_from_equilibrium() const{
  nvector* eq_resources = &(*eq_vals)[0][0];
  nvector* eq_consumers = &(*eq_vals)[0][1];

  Dynamical_variables dyn(eq_resources, eq_consumers);
  write_time_evolution(dyn, this->metaparameters->tf);
  return ;
}
Dynamical_variables CRModel::perturb_equilibrium() const{
  nvector unp_resources = (*eq_vals)[0][0];
  nvector unp_consumers = (*eq_vals)[0][1];

  nvector* p_resources = new nvector;
  nvector* p_consumers = new nvector;

  ntype eps = this->metaparameters->perturb_eq;

  // we perturb the equilibrium values
  for(size_t nu=0; nu < unp_resources.size(); ++nu){
    std::uniform_real_distribution<ntype> resource_dist((1.-eps)*unp_resources[nu], (1.+eps)*unp_resources[nu]);
    p_resources->push_back(resource_dist(random_engine));
  }

  for(size_t i = 0; i < unp_consumers.size(); ++i){
    std::uniform_real_distribution<ntype> consumer_dist((1.-eps)*unp_consumers[i], (1.+eps)*unp_consumers[i]);
    p_consumers->push_back(consumer_dist(random_engine));
  }

  return Dynamical_variables(p_resources, p_consumers);

}

void CRModel::perturb_parameters(const ntype & Delta) const{
  Parameter_set* p = this->model_param->get_parameters();
  std::uniform_real_distribution<ntype> uniform_distrib(-1., 1.);

  for(size_t mu=0; mu < p->l.size(); ++mu){
    p->l[mu] = p->l[mu]*(1+Delta*uniform_distrib(random_engine));
  }

  return;
}

void CRModel::perturb_parameters() const{
  this->perturb_parameters(this->metaparameters->perturb_parameters);
  return;
}

void CRModel::save_new_equilibrium(const Extinction& ext) const{
  Metaparameters* m = this->metaparameters;
  std::ofstream myfile;
  myfile.open(m->save_path, std::ios::app);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << m->save_path << " to write the new equilibrium of the system" << std::endl;
  }else{
    save_success = true;
  }
  myfile << "# "  << m->p << " " << m->epsilon << " " << m->foodmatrixpath << " " << m->verbose << " ";
  myfile << m->energy_constraint << " " << m->budget_constraint << " " << m->nb_attempts<<" " << m->seed_number <<" ";
  myfile << m->save_path << std::endl;
  myfile << m->NR << " " << m->NS << " " << m->perturb_parameters << " " << ext.t_eq << " " << ext.extinct  << std::endl;

  if(m->verbose > 0){
    if(save_success){
      std::cout << "Successfully saved new equilibrium of system in "<< m->save_path << std::endl;
    }
  }

  myfile.close();

  return;
}


/*** ALL THE ADDITIONAL USEFUL FUNCTIONS ****/
foodmatrix load_food_matrix(const Metaparameters& m){
  foodmatrix f(m.NS,nvector(m.NR, 0.));
  if(m.verbose > 0){
    std::cout << "Loading food matrix from " << m.foodmatrixpath << std::endl;
  }
  std::ifstream in(m.foodmatrixpath);
  if (!in) {
    std::cerr << "Cannot open file.\n" << "Returning empty food matrix" << std::endl;
    return f;
  }
  for (unsigned int x = 0; x < m.NS; x++) {
    for (unsigned int y = 0; y < m.NR; y++) {
      in >> f[x][y];
    }
  }
  in.close();
  return f;
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
gammamode string_to_gamma_mode(std::string mode){
  if(mode=="random_val"){
    return gammamode(random_val);
  }else if(mode=="nested"){
    return gammamode(nested);
  }else if(mode=="antinested"){
    return gammamode(antinested);
  }else{
    std::cerr << "Error, that value of gammamode has not been implemented yet or does not exist"<<std::endl;
    std::cerr << "Continuing with gammamode random_val" << std::endl;
    return gammamode(random_val);
  }
}
taumode string_to_tau_mode(std::string mode){
  if(mode=="tau0"){
    return taumode(tau0);
  }else if(mode=="taualpha"){
    return taumode(taualpha);
  }else{
    std::cerr << "Error, that value of taumode has not been implemented yet or does not exist"<<std::endl;
    std::cerr << "Continuing with taumode tau0" << std::endl;
    return taumode(tau0);
  }
}
alphamode string_to_alpha_mode(std::string mode){
  if(mode=="random_structure"){
    return alphamode(random_structure);
  }else if(mode=="no_release_when_eat"){
    return alphamode(no_release_when_eat);
  }else{
    std::cerr << "Error, that value of alphamode has not been implemented yet or does not exist"<<std::endl;
    std::cerr << "Continuing with alphamode random_structure" << std::endl;
    return alphamode(random_structure);
  }
}
std::ostream& display_matrix_w_name(std::ostream& os, std::string mat_name, const nmatrix & mat){
  os << mat_name << " = " << mat[0] << std::endl;
  for(size_t i = 1; i < mat.size(); ++i){
    os << std::setw(mat_name.size()+4+4) <<std::setfill(' ')<< mat[i] << std::endl;
  }
  return os;
}
std::ostream& display_vector_w_name(std::ostream& os, std::string vec_name, const nvector& vect){
  os << vec_name << " = " << vect << std::endl;
  return os;
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
  return;
}

void print_rand_number(){
  std::uniform_real_distribution<ntype> dist(0., 1.);
  std::cout << dist(random_engine) << std::endl;
}

/**** ALL THE GET FUNCTIONS for the classes *****/
Parameter_set* Model_parameters::get_parameters() const{
  return params;
}
nvector* Dynamical_variables::get_resources() const {
  return resources;
}
nvector* Dynamical_variables::get_consumers() const {
  return consumers;
}

/***** ALL THE SET FUNCTIONS for the classes *******/
void Model_parameters::set_sigma(const nmatrix& s) const {
  params->sigma = s;
  return;
}
void Model_parameters::set_alpha(const nmatrix& a) const {
  params->alpha = a;
  return;
}
void Model_parameters::set_gamma(const nmatrix& g) const {
  params->gamma = g;
  return;
}
void Model_parameters::set_tau(const nmatrix& t) const {
  params->tau = t;
  return;
}
void Model_parameters::set_l(const nvector& el) const {
  params->l = el;
  return;
}
void Model_parameters::set_d(const nvector& de) const {
  params->d = de;
  return;
}
void Model_parameters::set_m(const nvector& em) const {
  params->m = em;
  return;
}
void Model_parameters::set_NR(const unsigned int & enar) const {
  params->NR = enar;
  return;
}
void Model_parameters::set_NS(const unsigned int & enes) const {
  params->NS = enes;
  return;
}

/**** ALL THE operator<< overloads *****/
std::ostream& operator<<(std::ostream& os, const Metaparameters& m){
  os << "gamma0 = " << m.gamma0 << std::endl;
  os << "alpha0 = " << m.alpha0 << std::endl;
  os << "sigma0 = " << m.sigma0 << std::endl;
  os << "p = " << m.p << std::endl;
  os << "R0 = " << m.R0 << std::endl;
  os << "S0 = " << m.S0 << std::endl;
  os << "NR = " << m.NR << std::endl;
  os << "NS = " << m.NS << std::endl;
  return os;
}
std::ostream& operator<<(std::ostream& os, const nctype& number){
  os << std::real(number);
  if(std::imag(number)>0.){
    os << "+" << std::imag(number)<<"j";
  }else if(std::imag(number)<0.){
    os << std::imag(number)<<"j";
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const ncvector& vec){
  for(size_t i = 0; i <vec.size(); ++i){
    os << vec[i] << " ";
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const CRModel& M){
  M.display(os);
  return os;
}
std::ostream& operator<<(std::ostream& os, const Parameter_set& p){
  display_matrix_w_name(os, "gamma", p.gamma);
  //os << "mean = " << mean(p.gamma) << std::endl;
  display_matrix_w_name(os, "alpha", p.alpha);
  //os << "mean = " << mean(p.alpha) << std::endl;
  display_matrix_w_name(os, "tau", p.tau);
  //os << "mean = " << mean(p.tau) << std::endl;
  display_matrix_w_name(os, "sigma", p.sigma);
  //os << "mean = " << mean(p.sigma) << std::endl;
  display_vector_w_name(os, "l", p.l);
  display_vector_w_name(os, "m", p.m);
  display_vector_w_name(os, "d", p.d);
  os << "NR = " << p.NR << std::endl;
  os << "NS = " << p.NS << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const Model_parameters& M){
  os << *(M.get_parameters()) << std::endl;
  return os;
}
std::ostream& operator<<(std::ostream& os, const ntensor& T){
  for(size_t i = 0; i < T.size(); ++i){
    os << "Element " << i << " of tensor :" << std::endl << T[i] << std::endl;
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const nvector& v){
  for (size_t i=0; i < v.size(); ++i){
    os << std::right << std::fixed  << std::setprecision(print_precision) << v[i] << " ";
  }
  return os;
}
std::ostream& operator<<(std::ostream& os, const nmatrix& M){
  for (size_t i=0; i < M.size(); ++i){
    os << M[i] << std::endl;
  }
  return os;
}


Extinction_statistics compute_average_extinction(Metaparameters* metaparams, const ntype& Delta, unsigned int Nsimul){
  ntype convergence_threshold = 1e-6;
  Extinction_statistics av_extinct;
  av_extinct.t_eq.mean = 0.;
  av_extinct.t_eq.std_deviation= 0.;
  av_extinct.extinct.mean = 0.;
  av_extinct.extinct.std_deviation = 0.;
  av_extinct.new_Req.means = nvector(metaparams->NR, 0.);
  av_extinct.new_Seq.means = nvector(metaparams->NS, 0.);

  foodmatrix food_matrix = load_food_matrix(*metaparams);
  if(metaparams->verbose > 0){
    std::cout << "Computing average extinction for the given set of metaparameters." << std::endl;
  }

  double teq[Nsimul];
  double extinctions[Nsimul];

  for(size_t i = 0; i < Nsimul; ++i){
    CRModel model(food_matrix,*metaparams);
    model.perturb_parameters(Delta);
    Extinction new_equilib = model.evolve_until_equilibrium(convergence_threshold);

    teq[i] = new_equilib.t_eq;
    extinctions[i] = new_equilib.extinct;

    for(size_t j = 0 ; j < metaparams->NS; ++j){
      av_extinct.new_Seq.means[j] += (new_equilib.new_Seq[j]/Nsimul);
    }
    for(size_t mu = 0; mu < metaparams->NR; ++mu){
      av_extinct.new_Req.means[mu] += (new_equilib.new_Req[mu]/Nsimul);
    }
  }

  av_extinct.t_eq.mean = gsl_stats_mean(teq, 1, Nsimul);
  av_extinct.t_eq.std_deviation = gsl_stats_sd_m(teq, 1, Nsimul, av_extinct.t_eq.mean);

  av_extinct.extinct.mean = gsl_stats_mean(extinctions, 1, Nsimul);
  av_extinct.extinct.std_deviation = gsl_stats_sd_m(extinctions, 1, Nsimul, av_extinct.extinct.mean);

  if(metaparams->verbose > 0){
    std::cout << " Average extinction for Delta = " << Delta << " is " << av_extinct.extinct.mean;
    std::cout << " +/- " << av_extinct.extinct.std_deviation << std::endl;
  }
  return av_extinct;
}

void write_av_number_extinctions_delta_interval(Metaparameters* m, const nvector& deltas, unsigned int Nsimul)
{
  unsigned int Npoints = deltas.size();
  nvector av_extinctions = nvector(Npoints, 0.);
  nvector std_extinctions = nvector(Npoints, 0.);

  std::ofstream myfile;
  myfile.open(m->save_path,std::ios_base::app);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open "<<m->save_path <<" for saving the simulation "<< std::endl;
  }else{
    save_success=true;
    for(size_t i = 0 ; i < deltas.size(); ++i){

      Extinction_statistics ext_stat = compute_average_extinction(m, deltas[i], Nsimul);

      myfile << "# "  << m->p << " " << m->epsilon << " " << m->foodmatrixpath << " " << m->verbose << " ";
      myfile << m->energy_constraint << " " << m->budget_constraint << " " << m->nb_attempts<<" " << m->seed_number <<" ";
      myfile << m->save_path << " ";
      time_t now=time(0);
      myfile << ctime(&now) ;
      myfile << deltas[i] << " " << ext_stat.extinct.mean << " " << ext_stat.extinct.std_deviation << std::endl;
    }
  }
  myfile.close();

  return;
}
