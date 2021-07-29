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
    mcsolv.write_mode="none";


    // coprophagy is by convention allowed
    mcsolv.iss_allowed = metaparams.intra_specific_syntrophy;

    std::cout << "Running Monte Carlo Solver with the following parameters : " << std::endl;
    std::cout << metaparams << std::endl;
    std::string add_string = metaparams.save_path;

    std::cout << "Intraspecific syntrophy is";
    if(not(mcsolv.iss_allowed)){
      std::cout << " not";
    }

    nvector target_connectance = linear_interval(0.08, 0.43, 8);
    nvector target_nestedness = linear_interval(0.1, 0.6, 11);

    std::cout << " allowed during this run." << std::endl;
    for(auto target_conn : target_connectance){
      for(size_t i=0; i < target_nestedness.size() && target_nestedness[i] > target_conn; ++i){
        ntype target_nest = target_nestedness[i];
        mcsolv.additional_params=&target_nest;
        std::cout << "MC Mode = " << mcmode_to_string(mcsolv.mcmode) << std::endl;
        mcsolv.T=T0;

        EcologicalNetwork eco_net(metaparams.NR, metaparams.NS, target_conn);
        eco_net.optimize(mcsolv);


        //save path = save folder
        std::string savepath = metaparams.save_path+"/RandTrix_Nr"+metaparams.NR+"_Nc"+metaparams.NS+"_Nest"+nestedness(eco_net.G)+"_Conn"+connectance(eco_net.G)+".txt";
        std::ofstream savestream = open_external_file_truncate(savepath);
        display_food_matrix(savestream, eco_net.G);
        std::cout << "A target G-matrix was found and saved in " << savepath << std::endl;

        savestream.close();
      }
    }
  }catch(error e){
    e.handle();
  }

  return 0;
}
