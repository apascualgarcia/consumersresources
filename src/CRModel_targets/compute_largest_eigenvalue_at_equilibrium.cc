#include "CRModel.h"

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices=load_food_matrix_list(metaparams.foodmatrixpath);
    std::ofstream myfile = open_external_file_truncate(metaparams.save_path);
    unsigned int Npoints(100), Nsimuls(1e6);

    /* first compute the largest alpha feasible for every matrix considered */
    ntype alpha_max=1e9;
    for(size_t i=0; i < matrices.size();++i){
      metaparams.foodmatrixpath=matrices[i];
      ntype local_alpha_max=metaparams.feasible_alpha_max();
      if(local_alpha_max < alpha_max){
        alpha_max=local_alpha_max;
      }
    }
    std::cout << "Maximum common alpha0 is " << alpha_max << std::endl;

    /* now that we have found alpha max we can compute the maximal eigenvalue observed according to alpha */
    nvector alpha_range = linear_interval(0., alpha_max, Npoints);
    for(size_t i=0; i < matrices.size(); ++i){
      metaparams.foodmatrixpath = matrices[i];
      myfile << metaparams.foodmatrixpath << " ";
      for(size_t j=0; j < Npoints; ++j){
        metaparams.alpha0 = alpha_range[j];
        CRModel  model(metaparams);
        nctype largest_eigenvalue_at_equilibrium=model.largest_eigenvalue_at_equilibrium();
        for(size_t k=0; k < Nsimuls; ++k){
          CRModel model_test(metaparams);
          nctype local_largest = model_test.largest_eigenvalue_at_equilibrium();
          if(local_largest > largest_eigenvalue_at_equilibrium){
            largest_eigenvalue_at_equilibrium = local_largest;
          }
        }
        myfile << metaparams.alpha0 << " " << largest_eigenvalue_at_equilibrium << " ";
      }
      myfile << std::endl;
    }
    myfile.close();

  }catch(error e){
    e.handle();
  }

  return 0;
}
