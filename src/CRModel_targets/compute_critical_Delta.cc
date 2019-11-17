#include <iostream>
#include "CRModel.h"
#include <fstream>
using namespace std;

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);

  vector<string> matrices_path;
  std::ifstream in(metaparams.foodmatrixpath);
  if (!in) {
    std::cerr << "Cannot open file.\n" << " to compute critical deltas" << std::endl;
  }else{
    do{
      std::string a;
      in >> a;
      matrices_path.push_back(a);
    }while(!in.eof());
  }
  in.close();
  std::ofstream myfile;
  myfile.open(metaparams.save_path, std::ios::app);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << metaparams.save_path << " to write the new equilibrium of the system" << std::endl;
  }else{
    save_success = true;
    for(size_t i = 0; i < matrices_path.size();++i){
        metaparams.foodmatrixpath = matrices_path[i];
        double delta = compute_critical_Delta(metaparams, 0.);
        myfile << matrices_path[i] << " " << delta << std::endl;
    }
  }
  myfile.close();
  return 0;
}
