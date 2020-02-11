#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{

    Metaparameters metaparams(argc, argv);

    unsigned int Nsimuls=1000;

    ntype S_min=0.1, S_max=0.1;
    ntype S_points=1;

    ntype g_min=0.1, g_max=0.3;
    ntype g_points=100;

    nvector S_interval=linear_interval(S_min, S_max, S_points);
    nvector g_interval=linear_interval(g_min, g_max, g_points);

    std::ofstream myfile=open_external_file_truncate(metaparams.save_path);
    myfile << "# S0, gamma0, stability full, stability eff" << std::endl;

    for(size_t i=0; i < S_points;++i){
      for(size_t j=0; j < g_points; ++j){
        metaparams.S0=S_interval[i];
        metaparams.gamma0=g_interval[j];
        stability stab_full = compute_proportion_stability(metaparams, Nsimuls, full);
        stability stab_eff = compute_proportion_stability(metaparams, Nsimuls, effective);
        std::cout << "For metaparams " << metaparams << std::endl;
        std::cout << "  Full stability : \t " << stab_full << std::endl;
        std::cout << "  Effective stability : " << stab_eff << std::endl;
        myfile << metaparams.S0 << " "<<  metaparams.gamma0 << " " << stab_full << " " << stab_eff << std::endl;
      }
    }
    myfile.close();

  }catch(error e){
    e.handle();
  }
  return 0;
}
