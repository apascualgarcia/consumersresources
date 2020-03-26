#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list=load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder=metaparams.syntrophy_matrix_path;
    /*
    nvector expected_S0;
    for(auto mat : matrix_list){
      metaparams.foodmatrixpath=mat;
      foodmatrix gamma = load_food_matrix(metaparams);
      expected_S0.push_back(1./maximum(columns_degrees(gamma)));
      std::cout << "Matrix nest=" << nestedness(gamma) << ", conn=" << connectance(gamma) << " expected coeff=" << 1./maximum(columns_degrees(gamma))<<std::endl;
    }
    std::cout << "Expected S0 = [";
    for(auto S : expected_S0){
      std::cout << S << ",";
    }
    std::cout << "]" << std::endl;
    std::cout << "Expected coefficient : " << metaparams.l0*minimum(expected_S0) << std::endl;
    */
    ntype max_ratio=0.;
    for(auto mat: matrix_list){
      metaparams.foodmatrixpath=mat;
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
      CRModel model(metaparams);
      nmatrix A = model.get_parameter_set()->alpha;
      nmatrix G = model.get_parameter_set()->gamma;
      std::vector<unsigned int> A_deg = columns_degrees(A);
      std::vector<unsigned int> G_deg = row_degrees(G);
      ntype max_val=0.;
      for(size_t i=0; i < metaparams.NS; ++i){
        ntype local_val = ntype(A_deg[i])/ntype(G_deg[i]);
        if(local_val > max_val){
          max_val = local_val;
        }
      }
      if(max_val > max_ratio){
        max_ratio=max_val;
      }
      metaparams.syntrophy_matrix_path=syntrophy_folder;
  }
  std::cout << "Expected coefficient : " << metaparams.sigma0*metaparams.R0/max_ratio << std::endl;


  }catch(error e){
    e.handle();
  }

  return 0;
}
