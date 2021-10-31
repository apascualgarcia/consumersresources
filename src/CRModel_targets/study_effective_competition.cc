#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list = load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder = metaparams.syntrophy_matrix_path;

    ntype min_gamma0=0.01, max_gamma0=1.;
    ntype min_S0=0.01, max_S0=1.;
    ntype min_a0 = 0., max_a0 = 0.02;
    std::uniform_real_distribution<double> random_gamma(min_gamma0, max_gamma0);
    std::uniform_real_distribution<double> random_S(min_S0, max_S0);
    std::uniform_real_distribution<double> random_alpha(min_a0, max_a0);


    unsigned int Nsimuls=1e7;

    std::ofstream myfile = open_external_file_append(metaparams.save_path);
    myfile << "# We are writing, in that order, G-matrix location, proportion of feasible, stable, unstable, marginal, dominant eigenvalue of J, <C>";
    myfile << ", dominant eigenvalue of C, trace of C, dominant eigenvalue of B, effective comp defined as Eq.9 in APG 2017" << std::endl;

    for(auto mat: matrix_list){
      ntype prob_feasible=0., prob_stable = 0., prob_unstable= 0., prob_marginal = 0., av_dom_eig=0.;
      ntype av_av_comp=0., av_dom_c_eigval=0., av_c_trace=0., av_dom_b_eigval=0., av_inter_intra_ratio=0.;
      unsigned int N_dyn = 0;
      metaparams.foodmatrixpath=mat;
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
      myfile << metaparams.foodmatrixpath << " ";
      /* we are less conservative than before and will simply compute the percentage
           of feasible, stable, unstable and marginally stable systems in the [gamma0]x[S0] are */
      for(size_t i=0; i < Nsimuls; ++i){
        metaparams.gamma0=random_gamma(random_engine);
        metaparams.S0=random_S(random_engine);
        metaparams.alpha0=random_alpha(random_engine);

        CRModel model(metaparams, false);
        if(model.is_feasible()){
          prob_feasible+=1.;
          nmatrix C = model.get_effective_competition_matrix();
          av_av_comp+=mean(C);
          av_dom_c_eigval+=real(largest_eigenvalue(C));
          av_c_trace += trace(C);
          av_dom_b_eigval+=real(largest_eigenvalue(model.get_normalized_effective_competition_matrix()));
          av_inter_intra_ratio+=model.get_ratio_inter_intraspecific_competition();
          systemstability sys_stab = model.assess_dynamical_stability();
          switch(sys_stab){
            case stable:{
              prob_stable+=1.;
              N_dyn+=1;
              av_dom_eig+=real(model.largest_eigenvalue_at_equilibrium());
              break;
            }
            case marginal:{
              prob_marginal+=1.;
              break;
            }
            case unstable:{
              prob_unstable+=1.;
              break;
            }
            default:{
              throw error("Invalid system stability mode!");
              break;
            }
          }
        }
      }
      metaparams.syntrophy_matrix_path = syntrophy_folder;
      /* normalize to get probabilities */
      prob_feasible/=(1.*Nsimuls);
      prob_stable/=(1.*Nsimuls);
      prob_unstable/=(1.*Nsimuls);
      prob_marginal/=(1.*Nsimuls);
      av_dom_eig/=(1.*N_dyn);
      av_av_comp/=(1.*Nsimuls);
      av_dom_c_eigval/=(1.*Nsimuls);
      av_c_trace/=(1.*Nsimuls);
      av_dom_b_eigval/=(1.*Nsimuls);
      av_inter_intra_ratio/=(1.*Nsimuls);
      myfile << prob_feasible <<" "<< prob_stable << " " << prob_unstable << " " << prob_marginal;
      myfile << " " << av_dom_eig <<" " << av_av_comp <<  " " << av_dom_c_eigval << " " << av_c_trace;
      myfile  << " " << av_dom_b_eigval << " " << av_inter_intra_ratio<<  std::endl;

    }

    myfile.close();

  }catch(error e){
    e.handle();
  }
  return 0;
}
