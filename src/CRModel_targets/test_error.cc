#include "CRModel.h"

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    CRModel model1(metaparams);
    std::cout << "Model 1 : " << model1.get_parameter_set()->gamma << std::endl;
    std::cout << "Model 1 Adress in memory : " << &model1 << std::endl;
    std::cout << "Model 1 parameters : " << model1.get_model_parameters() << std::endl;
    std::cout << "----------" << std::endl;

    CRModel model2(metaparams);
    std::cout << "Model 2 : " << model2.get_parameter_set()->gamma << std::endl;
    std::cout << "Model 2 Adress in memory : " << &model2 << std::endl;
    std::cout << "Model 2 parameters : " << model2.get_model_parameters() << std::endl;
    std::cout << std::endl;
    std::cout<< " --- We now assign model 2 to model 1 ---" << std::endl;
    model1 = model2;
    std::cout << "Model 1 : " << model1.get_parameter_set()->gamma << std::endl;
    std::cout << "Model 1 Adress in memory : " << &model1 << std::endl;
    std::cout << "Model 1 parameters : " << model1.get_model_parameters() << std::endl;
    std::cout << std::endl;

    std::cout << "--- We create model 3 made from model 2 ---" << std::endl;

    CRModel model3(model2);
    std::cout << "Model 3 : " << model3.get_parameter_set()->gamma << std::endl;
    std::cout << "Model 3 Adress in memory : " << &model3 << std::endl;
    std::cout << "Model 3 parameters : " << model3.get_model_parameters() << std::endl;
    std::cout << std::endl;


  }catch(error e){
    e.handle();
  }
  return 0;
}
