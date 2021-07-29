/*******************************************************************************

This file allows to create consumption matrices with the target connectance
and nestedness.

Typical usage (from main folder):

build/generate_consumption_matrices PATH_TO_CONFIG_FILE

*******************************************************************************/

#include "../../include/CRModel.h"

/*********** CUSTOMIZABLE PART : energy which is minimized ***********/

ntype energy_function(const nmatrix& A, const nmatrix& G, void* params){
  //return quadratic_form_corrected_AlbertoMay2021(A, G, params);
  return quadratic_form_nestedness(A, G, params);
}

/*********** END OF THE CUSTOMIZABLE PART  ***********/


int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices_list=load_food_matrix_list(metaparams.foodmatrixpath);
    MonteCarloSolver mcsolv;
    ntype T0=0.1;
    mcsolv.max_steps=2000000;
    mcsolv.max_fails=300;
    mcsolv.annealing_freq=1000;
    mcsolv.annealing_const=1.-1e-2;
    mcsolv.display_stride=10000;
    mcsolv.cost_function=energy_function;
    mcsolv.mcmode=metaparams.mcmode;
    mcsolv.write_mode="converged_only";


    // coprophagy is by convention allowed
    mcsolv.iss_allowed = metaparams.intra_specific_syntrophy;

    std::cout << "Running Monte Carlo Solver with the following parameters : " << std::endl;
    std::cout << metaparams << std::endl;
    std::string add_string = metaparams.save_path;

    std::cout << "Intraspecific syntrophy is";
    if(not(mcsolv.iss_allowed)){
      std::cout << " not";
    }
    std::cout << " allowed during this run." << std::endl;

    for(size_t i=0; i < matrices_list.size();++i){

      ntype target_nest = 0.9, target_conn=0.2;
      mcsolv.additional_params=&target_nest;

      metaparams.foodmatrixpath=matrices_list[i];
      std::cout << "MC Mode = " << mcmode_to_string(mcsolv.mcmode) << std::endl;
      mcsolv.energy_file =metaparams.save_path+"_energy";

      std::ofstream smatrix_file=open_external_file_truncate(metaparams.save_path);
      smatrix_file << "# The following metaparameters were used for this matrix optimization : " << metaparams << std::endl;
      std::cout << "Starting the Monte Carlo algorithm to find the optimal syntrophy matrix for ";
      std::cout << metaparams.foodmatrixpath << std::endl;
      mcsolv.T=T0;

      EcologicalNetwork eco_net(metaparams.NR, metaparams.NS, target_conn);
      eco_net.optimize(mcsolv);
      display_food_matrix(smatrix_file, eco_net.A);
      if(mcsolv.mcmode==both_modified){
        std::ofstream gmatrix_file=open_external_file_truncate(metaparams.foodmatrixpath+"_optimized_"+add_string);
        gmatrix_file << "# The following metaparameters were used for this matrix optimization : " << metaparams << std::endl;
        display_food_matrix(gmatrix_file, eco_net.G);
        std::cout << "An optimal consumption matrix was found and saved in " << metaparams.foodmatrixpath+"_optimized_"+add_string << std::endl;
        gmatrix_file.close();
      }
      std::cout << "An optimal syntrophy matrix was found and saved in " << metaparams.save_path << std::endl;

      smatrix_file.close();

    }

  }catch(error e){
    e.handle();
  }

  return 0;
}
