#include "CRModel.h"

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <math.h>
#include <string>
#include <iomanip>
#include <ctime>
#include <stdlib.h>


/*** ALL THE CONSTRUCTORS for the classes ****/
Dynamical_variables::Dynamical_variables(nvector* resources_, nvector* consumers_){
  resources=resources_;
  consumers=consumers_;
  return;
}

/**** ALL THE GET FUNCTIONS for the classes *****/

nvector* Dynamical_variables::get_resources() const {
  return resources;
}
nvector* Dynamical_variables::get_consumers() const {
  return consumers;
}


void write_av_number_extinctions_delta_interval(Metaparameters* m, const nvector& deltas, unsigned int Nsimul){
  unsigned int Npoints = deltas.size();
  nvector av_extinctions = nvector(Npoints, 0.);
  nvector std_extinctions = nvector(Npoints, 0.);

  std::ofstream myfile;
  myfile.open(m->save_path,std::ios_base::app);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open "<<m->save_path <<" for saving the simulation "<< std::endl;
  }else{
    save_success=true;
    for(size_t i = 0 ; i < deltas.size(); ++i){

      Extinction_statistics ext_stat = compute_average_extinction(m, deltas[i], Nsimul);

      myfile << "# "  << m->p << " " << m->epsilon << " " << m->foodmatrixpath << " " << m->verbose << " ";
      myfile << m->energy_constraint << " " << m->budget_constraint << " " << m->nb_attempts<<" " << m->seed_number <<" ";
      myfile << m->save_path << " ";
      time_t now=time(0);
      myfile << ctime(&now) ;
      myfile << deltas[i] << " " << ext_stat.extinct.mean_ << " " << ext_stat.extinct.std_deviation_ << std::endl;
    }
  }
  myfile.close();
  return;
}
void write_prob_greater_than_one_delta_interval(Metaparameters* m, const nvector& deltas, unsigned int Nsimul){
  unsigned int Npoints = deltas.size();
  nvector probability_ext = nvector(Npoints, 0.);

  std::ofstream myfile;
  myfile.open(m->save_path,std::ios_base::app);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open "<<m->save_path <<" for saving the simulation "<< std::endl;
  }else{
    save_success=true;
    for(size_t i = 0 ; i < deltas.size(); ++i){

      double proba = probability_of_extinction_greather_than_one(m, deltas[i], Nsimul);

      myfile << "# "  << m->p << " " << m->epsilon << " " << m->foodmatrixpath << " " << m->verbose << " ";
      myfile << m->energy_constraint << " " << m->budget_constraint << " " << m->nb_attempts<<" " << m->seed_number <<" ";
      myfile << m->save_path << " ";
      time_t now=time(0);
      myfile << ctime(&now) ;
      myfile << deltas[i] << " " << proba << std::endl;
    }
  }
  myfile.close();
  return;

}
