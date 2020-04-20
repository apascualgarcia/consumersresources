#include "../../include/CRModel.h"



int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    ntype min_conn=1.1/metaparams.NR;
    if(min_conn > 1./metaparams.NS){
      min_conn = 1.1/metaparams.NS;
    }

    // std::vector<std::array<ntype,2>> nest_conn={{0.15,0.0832}, {0.15,0.128}, {0.15,0.176}, {0.1,0.0832}, {0.25,0.0848}, {0.25,0.1312}, {0.2,0.0848}, {0.2,0.1312}, {0.35,0.0832}, {0.35,0.1264}, {0.35,0.1712}, {0.3,0.0848}, {0.3,0.1296}, {0.3,0.1824}, {0.45,0.08}, {0.45,0.1296}, {0.45,0.1744}, {0.45,0.2336}, {0.4,0.0752}, {0.4,0.128}, {0.4,0.176}, {0.4,0.2176}, {0.55,0.0848}, {0.55,0.1232}, {0.55,0.168}, {0.55,0.216}, {0.55,0.2784}, {0.5,0.0816}, {0.5,0.128}, {0.5,0.1776}, {0.5,0.2208}, {0.5,0.2736}, {0.6,0.0816}, {0.6,0.1344}, {0.6,0.1712}, {0.6,0.2304}, {0.6,0.2768}};
    std::vector<std::array<ntype,2>> nest_conn={{0.4, 0.1}, {0.4, 0.2}, {0.4, 0.3}, {0.4, 0.35}, {0.4, 0.38}, {0.4, 0.4}};
    MonteCarloSolver mcsolv;
    mcsolv.max_steps=1000000;
    mcsolv.max_fails=1000;
    mcsolv.annealing_freq=1000;
    mcsolv.annealing_const=1.-1e-2;
    mcsolv.display_stride=10000;
    mcsolv.cost_function=quadratic_form_nestedness_rank;

    for(auto el : nest_conn){
      mcsolv.T=5.0;
      /* target nestedness is nest */
      ntype conn = el[1], nest=el[0];
      mcsolv.additional_params=&nest;
      std::cout << "Creating gamma matrix with target connectance " << conn ;
      std::cout << " and nestedness " << nest << std::endl;

      nmatrix gamma = optimal_consumption_matrix(metaparams.NR, metaparams.NS, conn, mcsolv);
      gamma=order_matrix_by_row_degree(order_matrix_by_column_degree(gamma));
      std::string mat_name="RandTrix_Nr"+std::to_string(metaparams.NR)+"_Nc"+std::to_string(metaparams.NS);
      mat_name+="_Nest";
      mat_name+=std::to_string(nestedness(gamma))+"_Conn"+std::to_string(connectance(gamma))+".txt";
      std::cout << "Saving matrix in " << "test_matrices/"+mat_name << std::endl;
      std::ofstream myfile=open_external_file_truncate("test_matrices/"+mat_name);
      display_food_matrix(myfile, gamma);
      myfile.close();

    }

  }catch(error e){
    e.handle();
  }

  return 0;
}
