#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{

    Metaparameters metaparams(argc, argv);

    ntype S_min=0.1, S_max=0.1;
    ntype S_points=1;

    ntype g_min=0.0, g_max=1.;
    ntype g_points=100;

    nvector S_interval=linear_interval(S_min, S_max, S_points);
    nvector g_interval=linear_interval(g_min, g_max, g_points);

    std::ofstream myfile=open_external_file_truncate(metaparams.save_path);

    for(size_t i=0; i < S_points;++i){
      for(size_t j=0; j < g_points; ++j){
        metaparams.S0=S_interval[i];
        metaparams.gamma0=g_interval[j];
        myfile << S_interval[i] << " " << g_interval[j] <<" "<<metaparams.feasible_alpha_max(1e-6) << std::endl;
      }
    }
    myfile.close();

  }catch(error e){
    e.handle();
  }
  return 0;
}
