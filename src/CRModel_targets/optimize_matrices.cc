#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices_list=load_food_matrix_list(metaparams.foodmatrixpath);

    ntype T=0.4;
    bool coprophagy=false;
    unsigned int max_steps=100000;
    unsigned int max_fails=1000;

    for(size_t i=0; i < matrices_list.size();++i){
      metaparams.foodmatrixpath=matrices_list[i];
      metaparams.save_path=optimal_alpha_matrix_path(metaparams.foodmatrixpath);

      foodmatrix gamma=load_food_matrix(metaparams);
      std::ofstream myfile=open_external_file_truncate(metaparams.save_path);

      std::cout << "Starting the Monte Carlo algorithm to find the optimal syntrophy matrix ... " << std::endl;
      foodmatrix optimal_alpha=optimal_syntrophy_from_consumption(gamma, T, coprophagy, max_steps, max_fails);
      display_food_matrix(myfile, optimal_alpha);
      std::cout << "An optimal syntrophy matrix was found and saved in " << metaparams.save_path << std::endl;

      myfile.close();

    }

  }catch(error e){
    e.handle();
  }

  return 0;
}
