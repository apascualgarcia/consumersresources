/*******************************************************************************

This file allows to create matrices which respect the needed energy conditions.
More specifically, it transforms each matrix of a set (the list of matrices
needs to be given as the path_to_food_matrix value of the input configuration file)
into a form which minimizes a given cost function energy_function (which can be)
changed on line 20). New energy functions can be added on the optimize_matrix file
from the CRModel_source folder.

Typical usage (from main folder):

build/optimize_matrices PATH_TO_CONFIG_FILE path_to_food_matrix=PATH_OF_MATRIX_LIST

*******************************************************************************/

#include "../../include/CRModel.h"

/*********** CUSTOMIZABLE PART : energy which is minimized ***********/

ntype energy_function(const nmatrix& A, const nmatrix& G, void* params){
  return quadratic_form_corrected_AlbertoMay2021(A, G, params);
}

/*********** END OF THE CUSTOMIZABLE PART  ***********/


int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices_list=load_food_matrix_list(metaparams.foodmatrixpath);
    MonteCarloSolver mcsolv;
    ntype T0=10.;
    mcsolv.max_steps=1000000;
    mcsolv.max_fails=1000;
    mcsolv.annealing_freq=1000;
    mcsolv.annealing_const=1.-1e-2;
    mcsolv.display_stride=10000;
    mcsolv.cost_function=energy_function;
;
    mcsolv.additional_params=&metaparams;

    /* set alpha0 to its maximal possible value */
    //metaparams.alpha0=metaparams.NR*metaparams.sigma0*metaparams.R0*metaparams.gamma0;

    // coprophagy is by convention allowed
    bool allow_coprophagy=true;

    for(size_t i=0; i < matrices_list.size();++i){
      metaparams.foodmatrixpath=matrices_list[i];
      metaparams.save_path=optimal_alpha_matrix_path(metaparams.foodmatrixpath);
      foodmatrix gamma=load_food_matrix(metaparams);
      std::ofstream myfile=open_external_file_truncate(metaparams.save_path+"_other");
      std::cout << "Starting the Monte Carlo algorithm to find the optimal syntrophy matrix for ";
      std::cout << metaparams.foodmatrixpath << std::endl;
      mcsolv.T=T0;
      foodmatrix optimal_alpha=optimal_syntrophy_from_consumption(gamma, allow_coprophagy, mcsolv);
      display_food_matrix(myfile, optimal_alpha);
      std::cout << "An optimal syntrophy matrix was found and saved in " << metaparams.save_path << std::endl;

      myfile.close();

    }

  }catch(error e){
    e.handle();
  }

  return 0;
}
