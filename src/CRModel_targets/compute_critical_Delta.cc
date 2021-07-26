#include <iostream>
#include "CRModel.h"
#include <fstream>
using namespace std;

int main(int argc, char * argv[]){
  try{

    Metaparameters metaparams(argc, argv);
    initialize_random_engine(metaparams);
    std::ofstream myfile = open_external_file_append(metaparams.save_path);
    if(metaparams.verbose > 0){
      std::cout << "The critical delta of " << metaparams.foodmatrixpath << " will be computed (alpha mode is "<< metaparams.alpha_mode << " and type_of_structural_perturbation=" << metaparams.struct_pert_type << ")" ;
      std::cout << std::endl;
    }

    metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
    delta_solver solv_params = {fitmode(polynomial),eqmode(oneextinct), stabilitymode(structural), metaparams.perturb_mode};

    if(metaparams.verbose >0){
      std::cout << "The metaparameters used are " << metaparams << std::endl;
    }

    std::time_t start, end;
    std::time(&start);
    statistics delta = compute_critical_Delta(metaparams, solv_params);
    std::time(&end);
    std::cout << "Computed critical delta for " << metaparams.foodmatrixpath;
    double time_taken = double(end-start);
    std::cout << " in " << time_taken << " seconds " << std::endl;
    std::cout << "The relative error achieved with those parameters is " <<  delta.std_deviation_/delta.mean_*100. << "%" << std::endl;
    myfile << metaparams.foodmatrixpath << " " << delta.mean_ << " " << delta.std_deviation_ << std::endl;

    myfile.close();
  }catch(error e){
    e.handle();
  }

  return 0;
}
