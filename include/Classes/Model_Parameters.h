#ifndef MODEL_PARAMETERS_H
#define MODEL_PARAMETERS_H

class Model_parameters{
private:
  Parameter_set  params;
public:
  Model_parameters();
  Model_parameters(const Model_parameters&);
  ~Model_parameters();

  // this function allows to give ANOTHER relationship between the parameters
  // (e.g. if tau is for instance zero or given by already existing parameters)
  Parameter_set* get_parameters();

  Parameter_set get_parameter_set() const;


  void display(std::ostream& ) const;

  // all the set functions
  void set_sigma(const nmatrix&) ;
  void set_alpha(const nmatrix&) ;
  void set_gamma(const nmatrix&) ;
  void set_tau(const nmatrix&) ;
  void set_l(const nvector&) ;
  void set_m(const nvector&) ;
  void set_d(const nvector&) ;
  void set_NR(const unsigned int&) ;
  void set_NS(const unsigned int&) ;
};

#endif
