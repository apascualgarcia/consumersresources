#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{
    /*
    Metaparameters metaparams(argc, argv);
    CRModel model(metaparams);

    if(model.is_dynamically_stable()){
      std::cout << "Model is dynamically stable" << std::endl;
    }else{
      std::cout << "Model is not dynamically stable " << std::endl;
    }

    Parameter_set* p=model.get_parameter_set();
    nmatrix gamma = p->gamma;
    nmatrix alpha = p->alpha;
    */

    nmatrix gamma={
      {1,1,1},
      {0,1,1},
      {0,0,1}
    };
    nmatrix alpha={
      {0,1,1},
      {0,0,1},
      {0,0,0}
    };

    std::cout << "gamma = " << std::endl<< gamma << std::endl;
    std::cout << "alpha = " << std::endl <<  alpha << std::endl;
    std::cout << "interaction matrix =" <<  std::endl << gamma*transpose(gamma)<< std::endl;
    std::cout << "interaction matrix deviation =" <<  std::endl << -gamma*alpha<< std::endl;

  }catch(error e){
    e.handle();
  }
  return 0;
}
