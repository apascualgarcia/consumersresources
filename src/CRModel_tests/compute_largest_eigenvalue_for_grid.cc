#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list = load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder = metaparams.syntrophy_matrix_path;

    ntype min_gamma0=0.01, max_gamma0=1.;
    ntype min_S0=0.01, max_S0=1.;

    unsigned int Nsimuls=200;
    unsigned int Npoints=100;


    nvector gamma0_interval=linear_interval(min_gamma0, max_gamma0, Npoints);
    nvector S0_interval=linear_interval(min_S0, max_S0, Npoints);

    std::ofstream myfile_real = open_external_file_append(metaparams.save_path+"_real");
    std::ofstream myfile_imag = open_external_file_append(metaparams.save_path+"_imag");
    std::ofstream myfile_err = open_external_file_append(metaparams.save_path+"_err");
    for(auto mat: matrix_list){
      metaparams.foodmatrixpath=mat;
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
      myfile_real << metaparams.foodmatrixpath << " ";
      myfile_imag << metaparams.foodmatrixpath << " ";
      myfile_err << metaparams.foodmatrixpath << " ";

      for(auto g : gamma0_interval){
        for(auto S : S0_interval){
          metaparams.gamma0=g;
          metaparams.S0=S;
          std::vector<ntype> largest_eigenvalue_real;
          std::vector<ntype> largest_eigenvalue_imag;

          /* we only want to work in the fully feasible region */
          bool fully_feasible=true;
          CRModel model;
          model.create_model_parameters(metaparams);
          for(size_t i=0; i < Nsimuls && fully_feasible; ++i){
            model.attempt_to_build_model(load_food_matrix(metaparams), metaparams, 0);
            if(model.is_feasible()){
              nctype lv = model.largest_eigenvalue_at_equilibrium();
              largest_eigenvalue_real.push_back(real(lv));
              largest_eigenvalue_imag.push_back(imag(lv));

            }else{
              fully_feasible=false;
              myfile_real << g << " " << S << " " << "NaN" << " ";
              myfile_imag << g << " " << S << " " << "NaN" << " ";
              myfile_err << g << " " << S << " " << "NaN" << " ";
              }
            }
            if(fully_feasible){
              myfile_real << g << " " << S << " " << mean(largest_eigenvalue_real) << " ";
              myfile_imag << g << " " << S << " " << mean(largest_eigenvalue_imag) << " ";
              myfile_err << g << " " << S << " " << sqrt(variance(largest_eigenvalue_real)) << " ";
            }
        }
      }
      metaparams.syntrophy_matrix_path=syntrophy_folder;
      myfile_err << std::endl;
      myfile_real << std::endl;
      myfile_imag << std::endl;

    }

    myfile_err.close();
    myfile_real.close();
    myfile_imag.close();


  }catch(error e){
    e.handle();
  }
  return 0;
}
