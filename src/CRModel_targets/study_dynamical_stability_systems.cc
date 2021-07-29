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

    unsigned int Nsimuls=200;
    unsigned int Npoints=100;

    nvector gamma0_interval=linear_interval(min_gamma0, max_gamma0, Npoints);
    nvector S0_interval=linear_interval(min_S0, max_S0, Npoints);

    std::ofstream myfile = open_external_file_append(metaparams.save_path);

    unsigned int total_systems = Nsimuls*gamma0_interval.size()*S0_interval.size();

    for(auto mat: matrix_list){
      ntype prob_feasible=0., prob_stable = 0., prob_unstable= 0., prob_marginal = 0.;
      metaparams.foodmatrixpath=mat;
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
      myfile << metaparams.foodmatrixpath << " ";
      for(auto g : gamma0_interval){
        for(auto S : S0_interval){
          metaparams.gamma0=g;
          metaparams.S0=S;

          /* we are less conservative than before and will simply compute the percentage
             of feasible, stable, unstable and marginally stable systems in the [gamma0]x[S0] are */

          for(size_t i=0; i < Nsimuls; ++i){
            CRModel model(metaparams);
            if(model.is_feasible()){
              prob_feasible+=1.;
              systemstability sys_stab = model.assess_dynamical_stability();
              switch(sys_stab){
                case stable:{
                  prob_stable+=1.;
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
        }
      }
      /* normalize to get probabilities */
      prob_feasible/=total_systems;
      prob_stable/=total_systems;
      prob_unstable/=total_systems;
      prob_marginal/=total_systems;
      std::cout << "feasible = " << prob_feasible << " stable= " << prob_stable << " unstable= " << prob_unstable << " marginal=" << prob_marginal << std::endl;

    }

    myfile.close();

  }catch(error e){
    e.handle();
  }
  return 0;
}
