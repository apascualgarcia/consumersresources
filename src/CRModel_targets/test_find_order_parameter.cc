#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);

    ntype average_d_rate=0.;
    ntype average_m_rate=0.;
    unsigned int Nsimul=100;
    ntype average_largest_eigenval=0.;
    ntype average_Rc=0.;
    for(size_t i=0; i < Nsimul; ++i){
      CRModel model(metaparams);
      nvector d=model.get_parameter_set()->d;
      nvector m=model.get_parameter_set()->m;
      nctype largest_eigenval=model.largest_eigenvalue_at_equilibrium();

      average_d_rate+=mean(d)/Nsimul;
      average_m_rate+=mean(m)/Nsimul;
      average_largest_eigenval+=real(largest_eigenval)/Nsimul;
      average_Rc+=model.critical_radius()/Nsimul;
    }
    std::cout << "-----------------" << std::endl;
    std::cout << "Connectance of gamma matrix : " << connectance(CRModel(metaparams).get_parameter_set()->gamma) << std::endl;
    std::cout << "Average consumers death rate : " << average_d_rate << std::endl;
    std::cout << "Average resource diffusion rate : " << average_m_rate << std::endl;
    std::cout << "Average largest real eigenvalue : " << average_largest_eigenval << std::endl;
    std::cout << "Average Rc : " << average_Rc << std::endl;


  }catch(error e){
    e.handle();
  }

  return 0;
}
