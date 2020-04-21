#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list=load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder = metaparams.syntrophy_matrix_path;
    for(auto mat : matrix_list){
      metaparams.foodmatrixpath=mat;
      std::cout << "Gamma matrix " << mat  << " has with nestedness ";
      CRModel model(metaparams);
      foodmatrix gamma = model.get_parameter_set()->gamma;
      std::cout << nestedness(gamma) << " and connectance ";
      std::cout << connectance(gamma) << std::endl;

      // std::ofstream myfile=open_external_file_truncate(metaparams.save_path);
      // display_food_matrix(myfile, alpha);
      // myfile.close();
      //std::cout << connectance(gamma) << ",";
    }

  }catch(error e){
    e.handle();
  }

  return 0;
}
