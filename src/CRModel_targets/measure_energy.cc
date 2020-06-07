#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices=load_food_matrix_list(metaparams.foodmatrixpath);
    std::ofstream myfile=open_external_file_truncate(metaparams.save_path);
    std::string syntrophy_folder = metaparams.syntrophy_matrix_path;


    myfile << "# Connectance of G, Nestedness of G, ";
    myfile << " Connectance of A, Nestedness of A, Energy E(G,A) " << std::endl;

    for(size_t i=0; i < matrices.size(); ++i){
      metaparams.foodmatrixpath=matrices[i];
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);

      CRModel model(metaparams);
      nmatrix gamma = model.get_parameter_set()->gamma;
      nmatrix alpha = model.get_parameter_set()->alpha;

      myfile << connectance(gamma) << " " << nestedness(gamma) << " ";
      myfile << connectance(alpha) << " " << nestedness(alpha) << " ";
      myfile << quadratic_form_low_intra_resource_interaction(alpha, gamma, &metaparams);

      metaparams.syntrophy_matrix_path=syntrophy_folder;
      myfile << std::endl;

    }
    myfile.close();

  }catch(error e){
    e.handle();
  }

  return 0;
}
