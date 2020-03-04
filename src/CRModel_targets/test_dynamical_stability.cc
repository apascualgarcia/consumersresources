#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices = load_food_matrix_list(metaparams.foodmatrixpath);

    for(size_t i=0; i < matrices.size();++i){
      metaparams.foodmatrixpath=matrices[i];
      ButlerModel model(metaparams);

      nmatrix gamma=model.get_parameter_set()->gamma;
      nmatrix alpha=model.get_parameter_set()->alpha;
      nmatrix GammaBeta= model.get_Gamma_matrix()*model.get_Beta_matrix();

      std::cout << "---- FULL MODEL -----" << std::endl;
      //std::cout << model << std::endl;
      //std::cout << "det(gamma)=" << det(gamma) << std::endl;
      std::cout << "Largest eigenvalue at equilibrium : " << model.largest_eigenvalue_at_equilibrium() << std::endl;

    }
    //std::cout << "Maximum feasible alpha0 : " << metaparams.feasible_alpha_max(1e-7) << std::endl;
    // std::cout << "Gamma matrix : " << std::endl << gamma << std::endl;
    // std::cout << "Alpha matrix : " << std::endl << alpha << std::endl;
    /*
    std::cout << "Gamma=" << std::endl << model.get_Gamma_matrix() << std::endl;
    std::cout << "Beta=" << std::endl << model.get_Beta_matrix() << std::endl;
    std::cout << "Gamma*Beta = " << std::endl << GammaBeta << std::endl;
    std::cout << "det(Gamma*Beta) = " << det(GammaBeta) << std::endl;
    */


    /*
    EffectiveCRModel effmodel(model);
    std::cout << "---- EFFECTIVE MODEL -----" << std::endl;
    std::cout << "Eigenvalues at equilibrium : " << effmodel.eigenvalues_at_equilibrium() << std::endl;
    std::cout << "The model is " << effmodel.assess_dynamical_stability() << std::endl;
    */



    // std::cout << "Rough jacobian : " << std::endl << gamma*alpha-gamma*transpose(gamma) << std::endl;
    // std::cout << "gamma alpha : " << std::endl << gamma*alpha << std::endl;
    // std::cout << "-gamma gamma_t : " << std::endl << -gamma*transpose(gamma) << std::endl;




  }catch(error e){
    e.handle();
  }
  return 0;
}
