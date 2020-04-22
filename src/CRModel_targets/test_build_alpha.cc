#include "CRModel.h"

int main(int argc, char *argv[]){
  try{
    Metaparameters meta(argc, argv);
    std::vector<std::string> consumption_matrix_list = load_food_matrix_list(meta.foodmatrixpath);

    for(auto mat : consumption_matrix_list){
      meta.foodmatrixpath=mat;

      // RANDOM STRUCTURE BUT CONNECTANCE OF "no_release_when_eat"
      // meta.alpha_mode=no_release_when_eat;
      // CRModel model(meta);
      // nmatrix alpha = model.get_parameter_set()->alpha;
      // nmatrix new_alpha = random_binary_matrix_with_connectance(meta.NR, meta.NS, connectance(alpha));

      // NO INTRASPECIFIC SYNTROPHY WITH CONNECTANCE OF CONSUMPTION MATRIX
      // nmatrix gamma=load_food_matrix(meta);
      // nmatrix new_alpha = binary_matrix_no_intraspecific_syntrophy(gamma);
      // std::cout << "Connectance new alpha - connectance gamma =" << connectance(new_alpha)-connectance(gamma) << std::endl;

      // LRI MATRIX with connectance of "no_release_when_eat"
      meta.alpha_mode=no_release_when_eat;
      CRModel model(meta);
      nmatrix alpha=model.get_parameter_set()->alpha;
      nmatrix gamma=model.get_parameter_set()->gamma;
      nmatrix new_alpha = build_LRI_matrix(gamma, meta, connectance(alpha));

      std::ofstream myfile = open_external_file_truncate(optimal_alpha_matrix_path(mat));
      display_food_matrix(myfile, new_alpha);
      myfile.close();
    }

  }catch(error e){
    e.handle();
  }
  return 0;
}
