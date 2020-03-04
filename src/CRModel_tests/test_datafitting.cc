#include "../../include/CRModel.h"
#include <array>
int main(){
  try{
    nvector x_data;
    nvector y_data;

    std::array<double, 10> x_array = {0.0157463,0.0174171,0.0190878,0.0207586,0.0224293,0.0241001,0.0257708,0.0274416,0.0291123,0.0307831};

    std::array<double, 10> y_array ={-0.44,-0.38,-0.18,-0.26,-0.22,1.11022e-16,0.02,0.22,0.1,0.32};

    for(size_t i = 0; i < 10; ++i){
      x_data.push_back(x_array[i]);
      y_data.push_back(y_array[i]);
    }

    NUMBER_OF_FITTING_PARAMETERS = 2;

    gsl_vector* params = gsl_vector_alloc(NUMBER_OF_FITTING_PARAMETERS);
    gsl_vector* error = gsl_vector_alloc(NUMBER_OF_FITTING_PARAMETERS);

    fitting_parameters fit_params = {params, error};

    fit_points_with_function(x_data, y_data, fit_params, fitmode(sigmoidal));
    std::cout << "[";
    for(size_t i = 0; i < NUMBER_OF_FITTING_PARAMETERS-1; ++i){
      std::cout << gsl_vector_get(params, i) << "+/-"<< gsl_vector_get(error, i) << ",";
    }
    std::cout <<  gsl_vector_get(params, NUMBER_OF_FITTING_PARAMETERS-1) << "+/-"<< gsl_vector_get(error, NUMBER_OF_FITTING_PARAMETERS-1) << "]" << std::endl;
  }catch(error e){
    e.handle();
  }
  return 0;
}
