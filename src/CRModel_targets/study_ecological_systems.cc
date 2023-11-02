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
    std::uniform_real_distribution<double> random_gamma(min_gamma0, max_gamma0);
    std::uniform_real_distribution<double> random_S(min_S0, max_S0);

    //unsigned int Nsimuls=1e2; // for testing run
    //unsigned int Nsimuls = 1e5; // for low precision production run
    unsigned int Nsimuls = 1e6; // for high precision production run
    unsigned int spacing = Nsimuls/10;

    if(Nsimuls < 1e5){
      std::cout << "WARNING : number of simulations smaller than 100'000. Are you sure you want to run at a lower value?" << std::endl;
    }

    std::ofstream myfile = open_external_file_append(metaparams.save_path);
    myfile << "# We are writing, in that order, G-matrix location, A-matrix location, G connectance, G nestedness, A connectance, A nestedness, A-mode, alpha0 value, total number of simulations, proportion of feasible, stable, unstable, marginal and dominant eigenvalue" << std::endl;

    ntype G_connectance = 0., G_nestedness=0.;

    for(auto mat: matrix_list){
      ntype prob_feasible=0., prob_stable = 0., prob_unstable= 0., prob_marginal = 0., av_dom_eig=0., av_A_connectance=0., av_A_nestedness=0.,percentage_run = 0.;
      unsigned int N_dyn = 0;
      metaparams.foodmatrixpath=mat;
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
            /* we are less conservative than before and will simply compute the percentage
           of feasible, stable, unstable and marginally stable systems in the [gamma0]x[S0] are */
      for(size_t i=0; i < Nsimuls; ++i){
        metaparams.gamma0=random_gamma(random_engine);
        metaparams.S0=random_S(random_engine);

        percentage_run += 100./Nsimuls;

        if(metaparams.verbose>1){
          std::cout << "Run " << i << " out of "<< Nsimuls << ": ";
        }

        if(metaparams.verbose>0 && i%spacing==0){
          std::cout << int((100.0*i)/Nsimuls) << "% of the simulations have been run." << std::endl;
        }

        CRModel model(metaparams, false);
        if(i==0){
          G_connectance = connectance(model.get_G());
          G_nestedness = nestedness(model.get_G());
        }
        av_A_nestedness+=nestedness(model.get_A())/Nsimuls;
        av_A_connectance+=connectance(model.get_A())/Nsimuls;
        if(model.is_feasible()){
          prob_feasible+=1.;
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
      prob_feasible/=Nsimuls;
      prob_stable/=Nsimuls;
      prob_unstable/=Nsimuls;
      prob_marginal/=Nsimuls;
      av_dom_eig/=N_dyn;

      /* Write everything to output */
      myfile << metaparams.foodmatrixpath << " ";
      myfile << metaparams.syntrophy_matrix_path << " ";
      myfile << G_connectance << " ";
      myfile << G_nestedness << " ";
      myfile << av_A_connectance << " ";
      myfile << av_A_nestedness << " ";
      myfile << metaparams.alpha_mode << " ";
      myfile << metaparams.alpha0 << " ";
      myfile << Nsimuls << " ";
      myfile << prob_feasible <<" "<< prob_stable << " " << prob_unstable << " " << prob_marginal;
      myfile << " " << av_dom_eig << std::endl;

      if(metaparams.verbose>0){
        std::cout << "Done for matrix G = " << metaparams.foodmatrixpath << std::endl;
      }

    }

    myfile.close();

  }catch(error e){
    e.handle();
  }
  return 0;
}
