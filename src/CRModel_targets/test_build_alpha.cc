#include "CRModel.h"

int main(int argc, char *argv[]){
  try{
    Metaparameters meta(argc, argv);
    std::vector<std::string> consumption_matrix_list = load_food_matrix_list(meta.foodmatrixpath);

    for(auto mat : consumption_matrix_list){
      std::cout << mat << std::endl;
    }

  }catch(error e){
    e.handle();
  }
  return 0;
}
