#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{

    Metaparameters metaparams(argc, argv);

    ntype S_min=0.01, S_max=0.1;
    ntype S_points=100;

    ntype g_min=0., g_max=0.1;
    ntype g_points=100;

    nvector S_interval;
    for(size_t i=0; i < S_points; ++i){
      ntype point=S_min+(S_max-S_min)*i/(S_points-1.);
      S_interval.push_back(point);
    }

    nvector g_interval;
    for(size_t i=0; i < g_points;++i){
      ntype point = g_min+(g_max-g_min)*i/(g_points-1.);
      g_interval.push_back(point);
    }

    std::ofstream myfile=open_external_file_truncate(metaparams.save_path);

    for(size_t i=0; i < S_points;++i){
      for(size_t j=0; j < g_points; ++j){
        metaparams.S0=S_interval[i];
        metaparams.gamma0=g_interval[j];
        myfile << S_interval[i] << " " << g_interval[j] <<" "<<metaparams.feasible_alpha_max(1e-8) << std::endl;
      }
    }
    myfile.close();

  }catch(error e){
    e.handle();
  }
  return 0;
}
