#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list=load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder = metaparams.syntrophy_matrix_path;
    std::cout << "[";
    for(auto mat : matrix_list){
      metaparams.foodmatrixpath=mat;
      std::cout << "Gamma matrix " << mat  << " gives rise to alpha matrix with nestedness ";
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
      CRModel model(metaparams);
      foodmatrix alpha = model.get_parameter_set()->alpha;
      foodmatrix gamma = model.get_parameter_set()->gamma;
      std::cout << nestedness(alpha) << std::endl;

      metaparams.save_path="data_output/optimal_LRI_for_NR"+std::to_string(metaparams.NR)+"_NS"+std::to_string(metaparams.NS)+"_Nest"+std::to_string(nestedness(gamma))+"_Conn"+std::to_string(connectance(gamma))+".out";
      // std::ofstream myfile=open_external_file_truncate(metaparams.save_path);
      // display_food_matrix(myfile, alpha);
      // myfile.close();
      //std::cout << connectance(gamma) << ",";

      metaparams.syntrophy_matrix_path=syntrophy_folder;
    }
    std::cout << "]" << std::endl;

  }catch(error e){
    e.handle();
  }

  return 0;
}
