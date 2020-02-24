#include "CRModel.h"

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    metaparams.building_mode=use_m;
    metaparams.alpha0=0.;

    std::vector<std::string> matrices=load_food_matrix_list(metaparams.foodmatrixpath);
    std::ofstream myfile = open_external_file_truncate(metaparams.save_path);
    unsigned int Npoints(10), Nsimuls(1000);

    /* now that we have found alpha max we can compute the maximal eigenvalue observed according to alpha */
    nvector m_range = linear_interval(0., 2., Npoints);
    for(size_t i=0; i < matrices.size(); ++i){
      metaparams.foodmatrixpath = matrices[i];
      myfile << metaparams.foodmatrixpath << " ";
      for(size_t j=0; j < Npoints; ++j){
        metaparams.m0 = m_range[j];
        CRModel  model(metaparams);
        nctype largest_eigenvalue_at_equilibrium=model.largest_eigenvalue_at_equilibrium();
        for(size_t k=0; k < Nsimuls; ++k){
          CRModel model_test(metaparams);
          nctype local_largest = model_test.largest_eigenvalue_at_equilibrium();
          if(local_largest > largest_eigenvalue_at_equilibrium){
            largest_eigenvalue_at_equilibrium = local_largest;
          }
        }
        myfile << metaparams.m0 << " " << largest_eigenvalue_at_equilibrium << " ";
      }
      myfile << std::endl;
    }
    myfile.close();

  }catch(error e){
    e.handle();
  }

  return 0;
}
