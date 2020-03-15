#include "../../include/CRModel.h"



int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    ntype min_conn=1.1/metaparams.NR;
    if(min_conn > 1./metaparams.NS){
      min_conn = 1.1/metaparams.NS;
    }
    nvector nest_range=linear_interval(0.05, 0.6, 5);

    MonteCarloSolver mcsolv;
    mcsolv.T=10.;
    mcsolv.max_steps=1000000;
    mcsolv.max_fails=1000;
    mcsolv.annealing_freq=1000;
    mcsolv.annealing_const=1.-1e-2;
    mcsolv.display_stride=10000;
    mcsolv.cost_function=quadratic_form_nestedness;

    for(auto nest: nest_range){
      /* target nestedness is nest */
      mcsolv.additional_params=&nest;
      nvector conn_range=linear_interval(min_conn, nest, 5);
      for(auto conn: conn_range){

        std::cout << "Creating gamma matrix with target connectance " << conn ;
        std::cout << " and nestedness " << nest << std::endl;

        nmatrix gamma = optimal_consumption_matrix(metaparams.NR, metaparams.NS, conn, mcsolv);
        std::string mat_name="RandTrix_Nr"+std::to_string(metaparams.NR)+"_Nc"+std::to_string(metaparams.NS);
        mat_name+="_Nest";
        mat_name+=std::to_string(nestedness(gamma))+"_Conn"+std::to_string(connectance(gamma))+".txt";

        std::ofstream myfile=open_external_file_truncate("test_matrices/"+mat_name);
        display_food_matrix(myfile, gamma);
        myfile.close();
        
      }
    }

  }catch(error e){
    e.handle();
  }

  return 0;
}
