#include "../../include/CRModel.h"
#include <iostream>
int main(){

  nvector v1 = {0., 1.};
  nvector v2 = {1., 0.};

  std::cout << "Angle between v1 and v2 : " << angle(v1, v2) << std::endl;

  return 0;
}
