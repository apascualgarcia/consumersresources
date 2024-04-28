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


/*********** END OF THE CUSTOMIZABLE PART  ***********/


int main(int argc, char* argv[]){
  try{

      Metaparameters metaparams(argc, argv);
      std::vector<std::string> matrix_list = load_food_matrix_list(metaparams.foodmatrixpath);
      std::string add_string = metaparams.save_path;

      for (auto mat : matrix_list){
        
        metaparams.foodmatrixpath = mat;
        metaparams.save_path=optimal_alpha_matrix_path(metaparams.foodmatrixpath)+"_"+add_string;
        std::ofstream smatrix_file=open_external_file_truncate(metaparams.save_path);
        smatrix_file << "# The following metaparameters were used for this matrix optimization : " << metaparams << std::endl;

        MCSv2 mcs;
        Model_parameters model_params(metaparams);
        model_params.optimize(mcs, &metaparams);

        nmatrix optimized_alpha = model_params.get_parameter_set().alpha;

        smatrix_file << optimized_alpha;
        smatrix_file.close();

        std::cout << "An optimal syntrophy matrix was found and saved in " << metaparams.save_path << std::endl;

      }
  
  }catch(error e){
    e.handle();
  }

  return 0;
}
