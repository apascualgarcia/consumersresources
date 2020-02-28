#include "CRModel.h"
#include <algorithm>
nmatrix build_sigma(const Metaparameters& m){
  nmatrix sigma(m.NS, nvector(m.NR,0.));
  if(m.sigma0 > 0){
    std::uniform_real_distribution<ntype> sigma_distribution((1.-m.epsilon)*m.sigma0, (1.+m.epsilon)*m.sigma0);
    for(size_t i=0; i < sigma.size(); ++i){
      for(size_t mu=0; mu < sigma[i].size(); ++mu){
        sigma[i][mu] = sigma_distribution(random_engine);
      }
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

  if(m.R0 > 0){
    std::uniform_real_distribution<ntype> R_distrib((1.-m.epsilon)*m.R0, (1.+m.epsilon)*m.R0);
    for(size_t i=0; i < v.size(); ++i){
      v[i] = R_distrib(random_engine);
    }
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

  if(m.S0>0){
    std::uniform_real_distribution<ntype> S_distrib((1.-m.epsilon)*m.S0, (1.+m.epsilon)*m.S0);
    for(size_t i=0; i < v.size(); ++i){
      v[i] = S_distrib(random_engine);
    }
  }

  return v;
}
nmatrix build_gamma(const foodmatrix& F, const Metaparameters& m){
  nmatrix gamma(F.size(), nvector(F[0].size(), 0.));

  if(F.size()!=m.NS or F[0].size()!=m.NR){
    std::cerr << "Meta parameters not compatible with food matrix, please change this " << std::endl;
  }

  if(m.gamma0>0){
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
          error e("Problem in the food matrix "+ m.foodmatrixpath + ", species " + std::to_string(i) +" does not eat anything.");
          std::cerr << "The food matrix is " << F << std::endl;
          throw e;
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
          error e("Problem in the food matrix "+ m.foodmatrixpath + ", species " + std::to_string(i) +" does not eat anything.");
          std::cerr << "The food matrix is " << F << std::endl;
          throw e;
        }
      }
    }
  }
  return gamma;
}
nmatrix build_alpha(const Parameter_set* p, Metaparameters& m, const nvector& Req, unsigned int attempts){
  std::uniform_real_distribution<ntype> empty_or_not_distrib(0., 1.);
  std::uniform_real_distribution<ntype> alpha_distrib((1-m.epsilon)*m.alpha0, (1+m.epsilon)*m.alpha0);
  nmatrix alpha = nmatrix(p->NR, nvector(p->NS, 0.));
  if(m.alpha0 > 0){
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
            if(!(p->gamma[i][mu]>0.) || (ntype(empty_or_not_distrib(random_engine)) < m.p)){
              alpha[mu][i] = alpha_distrib(random_engine);
            }
          }
        }
        break;
        //rescale_mean(alpha, m.alpha0);
      }

      case one_release:{
        for(size_t i=0; i < p->NS; ++i){
          std::vector<unsigned int> uneaten_resources;
          for(size_t mu=0; mu < p->NR; ++mu){
            if(!(p->gamma[i][mu]>0)){
              uneaten_resources.push_back(mu);
            }
          }
          if(uneaten_resources.size()>0){
            std::uniform_int_distribution<size_t> unif_int_dist(0, uneaten_resources.size()-1);
            alpha[uneaten_resources[unif_int_dist(random_engine)]][i]=alpha_distrib(random_engine);
          }
        }
        break;
      }

      case optimal_matrix:{
        alpha=load_syntrophy_matrix(m);
        for(size_t mu=0; mu < p->NR; ++mu){
          for(size_t i=0; i < p->NS; ++i){
            if(alpha[mu][i]>0.){
              alpha[mu][i] = alpha_distrib(random_engine);
            }
          }
        }


      }

      default:{
        error e("This alpha mode has not been implemented yet.");
        throw e;
        break;
      }
    }

  }

  /* OLD CODE from Jan 27
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
  */
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
      throw error("Tau value not implemented yet.");
      break;
    }
  };
  return tau;
}
nvector build_l(const Metaparameters & m){
  nvector v(m.NR, 0.);
  if(m.l0>0){
    std::uniform_real_distribution<ntype> l_distrib((1.-m.epsilon)*m.l0, (1.+m.epsilon)*m.l0);
    for(size_t nu=0; nu < m.NR; ++nu){
      v[nu] = l_distrib(random_engine);
    }
  }
  return v;
}

nvector build_m(const Metaparameters& meta){
  nvector m(meta.NR, 0);
  if(meta.m0 > 0){
    std::uniform_real_distribution<ntype> m_distrib((1-meta.epsilon)*meta.m0, (1+meta.epsilon)*meta.m0);
    for(size_t nu=0; nu < meta.NR; ++nu){
      m[nu] = m_distrib(random_engine);
    }
  }
  return m;
}

nmatrix build_sigma_Butler(const Metaparameters& m){
  nmatrix epsilon(m.NS, nvector(m.NR,0));
  for(size_t i=0; i < m.NS; ++i){
    for(size_t nu=0; nu < m.NR; ++nu){
      epsilon[i][nu]=m.sigma0;
    }
  }
  /* In the model of Butler, the efficiency is actually constant, i.e. sigma_{inu}\gamma{i\nu}=epsilon \gamma_{i\nu} */

  return epsilon;
}
