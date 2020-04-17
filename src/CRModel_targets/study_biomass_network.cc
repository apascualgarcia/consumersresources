#include "CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters mets(argc,argv);
    CRModel model(mets);
    nmatrix bio_network = model.get_biomass_flux_network();

    display_food_matrix(std::cout, bio_network);
    std::cout << std::endl;
    std::cout << "Assortativity of network = " << assortativity(bio_network) << std::endl;

  }catch(error e){
    e.handle();
  }

  return 0;
}
