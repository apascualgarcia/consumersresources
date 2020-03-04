#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices = load_food_matrix_list(metaparams.foodmatrixpath);

    unsigned int N_matrices = matrices.size();
    for(size_t i=0; i < N_matrices; ++i){
      metaparams.foodmatrixpath = matrices[i];
      metaparams.save_path= matrices[i]+"_corr";
      std::ofstream myfile = open_external_file_truncate(metaparams.save_path);
      foodmatrix f = load_food_matrix(metaparams);
      display_food_matrix(myfile, order_matrix_by_column_degree(order_matrix_by_row_degree(f)));
      myfile.close();

    }

  }catch(error e){
    e.handle();
  }


  return 0;
}
