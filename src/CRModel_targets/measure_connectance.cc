#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices=load_food_matrix_list(metaparams.foodmatrixpath);
    std::ofstream myfile=open_external_file_truncate(metaparams.save_path);
    myfile << "# First column is the location of the food consumption (gamma) matrix" << std::endl;
    myfile << "# Second column is the connectance of the corresponding syntrophy (alpha) matrix with the rule that species release to every resource they do not consume" << std::endl;

    for(size_t i=0; i < matrices.size(); ++i){
      metaparams.foodmatrixpath=matrices[i];
      myfile << metaparams.foodmatrixpath << " ";
      CRModel model(metaparams);
      myfile << connectance(model.get_parameter_set()->alpha) << std::endl;
    }
    myfile.close();

  }catch(error e){
    e.handle();
  }

  return 0;
}
