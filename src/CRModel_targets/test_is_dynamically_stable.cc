#include "CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    CRModel model(metaparams);

    std::cout << "Eigenvalues at equilibrium : " << model.eigenvalues_at_equilibrium() << std::endl;
    std::cout << "Eigenvalues computed analytically : ";

    ntype A = -(metaparams.l0+metaparams.NS*metaparams.alpha0*metaparams.S0);
    ntype B = metaparams.alpha0-metaparams.gamma0*metaparams.R0;
    ntype C = metaparams.sigma0*metaparams.gamma0*metaparams.S0;
    std::cout << A << " ";
    ntype Delta2 = A*A -4*B*C*metaparams.NR*metaparams.NS;
    if(Delta2 > 0){
      std::cout << 0.5*(A+sqrt(Delta2)) << " " << 0.5*(A-sqrt(Delta2)) ;
    }else{
      if(Delta2 <0){
        std::cout << 0.5*A << "-" << 0.5*sqrt(-Delta2) << "i "  << 0.5*A << "+" << 0.5*sqrt(-Delta2) << "i";
      }else{
        std::cout << 0.5*A;
      }
    }
    std::cout << std::endl;

  }catch(error e){
    e.handle();
  }
  return 0;
}
