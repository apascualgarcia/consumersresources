#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);

    CRModel model(metaparams);
    EffectiveCRModel effmodel(model);
    nmatrix gamma=model.get_parameter_set()->gamma;
    nmatrix alpha=model.get_parameter_set()->alpha;
  // std::cout << "Maximum feasible alpha0 : " << metaparams.feasible_alpha_max(1e-4) << std::endl;
    // std::cout << "Gamma matrix : " << std::endl << gamma << std::endl;
    // std::cout << "Alpha matrix : " << std::endl << alpha << std::endl;
    std::cout << "Gamma*Beta = " << model.get_Beta_matrix()*model.get_Gamma_matrix() << std::endl;

    std::cout << "---- FULL MODEL -----" << std::endl;
    std::cout << "Eigenvalues at equilibrium : " << model.eigenvalues_at_equilibrium() << std::endl;
    std::cout << "The model is " << model.assess_dynamical_stability() << std::endl;

    std::cout << "---- EFFECTIVE MODEL -----" << std::endl;
    std::cout << "Eigenvalues at equilibrium : " << effmodel.eigenvalues_at_equilibrium() << std::endl;
    std::cout << "The model is " << effmodel.assess_dynamical_stability() << std::endl;


    // std::cout << "Rough jacobian : " << std::endl << gamma*alpha-gamma*transpose(gamma) << std::endl;
    // std::cout << "gamma alpha : " << std::endl << gamma*alpha << std::endl;
    // std::cout << "-gamma gamma_t : " << std::endl << -gamma*transpose(gamma) << std::endl;




  }catch(error e){
    e.handle();
  }
  return 0;
}
