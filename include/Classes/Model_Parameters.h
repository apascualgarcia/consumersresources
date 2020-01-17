#ifndef MODEL_PARAMETERS_H
#define MODEL_PARAMETERS_H

class Model_parameters{
private:
  Parameter_set*  params;
public:
  Model_parameters();
  ~Model_parameters();
  // this function allows to give ANOTHER relationship between the parameters
  // (e.g. if tau is for instance zero or given by already existing parameters)
  Parameter_set* get_parameters() const;

  // all the set functions
  void set_sigma(const nmatrix&) const;
  void set_alpha(const nmatrix&) const;
  void set_gamma(const nmatrix&) const;
  void set_tau(const nmatrix&) const;
  void set_l(const nvector&) const;
  void set_m(const nvector&) const;
  void set_d(const nvector&) const;
  void set_NR(const unsigned int&) const;
  void set_NS(const unsigned int&) const;
};

#endif
