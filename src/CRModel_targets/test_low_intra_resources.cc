#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    /* the numbers of dynamically stable systems we want */
    unsigned int Nsimuls=1000, step=0;
    ntype proba_low_interaction=0.;

    std::cout << "Nestedness of food matrix " << nestedness(load_food_matrix(metaparams)) << std::endl;

    do{
      CRModel model(metaparams);
      if(model.is_dynamically_stable()){
        step+=1;
        if(model.is_in_low_intra_resource_interaction()){
          proba_low_interaction+=1./Nsimuls;
        }
      }
    }while(step<Nsimuls);

    std::cout << "The probability of being in low intra resources regime (given that you are already dynamically stable) is ";
    std::cout << proba_low_interaction << std::endl;


  }catch(error e){
    e.handle();
  }
  return 0;
}
