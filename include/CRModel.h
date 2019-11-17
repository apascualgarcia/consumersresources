#ifndef CRMODEL_H
#define CRMODEL_H
#include<vector>
#include<string>
#include<complex>
#include<random>
#include <gsl/gsl_vector.h>


typedef long double ntype;
typedef std::complex<ntype> nctype;
typedef std::vector<ntype> nvector;
typedef std::vector<nctype> ncvector;
typedef std::vector<nvector> nmatrix;
typedef std::vector<nmatrix> ntensor;
typedef nmatrix foodmatrix;

// tau0 : tau = 0; taualpha : tau = alpha
enum taumode{tau0,taualpha};
enum gammamode{random_val, nested, antinested};
enum alphamode{random_structure, no_release_when_eat};

const unsigned int print_precision=4;
const ntype EIGENSOLVER_PRECISION = 1e-15;

struct Parameter_set{
  nmatrix alpha;
  nmatrix gamma;
  nmatrix tau;
  nmatrix sigma;
  nvector l;
  nvector m;
  nvector d;
  unsigned int NR;
  unsigned int NS;
};
struct Metaparameters{
  ntype gamma0;
  ntype alpha0;
  ntype sigma0;
  ntype p;
  ntype R0;
  ntype S0;
  ntype l0;
  ntype epsilon;
  unsigned int NR;
  unsigned int NS;
  gammamode gamma_mode;
  taumode tau_mode;
  alphamode alpha_mode;
  std::string foodmatrixpath;
  bool verbose;
  bool energy_constraint;
  bool budget_constraint;
  unsigned int nb_attempts;
  unsigned int seed_number;
  ntype tf;
  std::string save_path;
  ntype perturb_eq;
  ntype perturb_parameters;

  Metaparameters(int argc, char *argv[]);
};

struct Extinction{
  ntype t_eq;
  ntype extinct;
  nvector new_Req;
  nvector new_Seq;
};

struct Solver_Parameters{
  // these dictate the values we will take for compute_average_extinction
  Metaparameters* metaparameters;
  unsigned int Nsimul;
};

struct Delta_critical{
  ntype delta_crit;
  ntype delta_low;
  ntype delta_high;
  ntype accuracy;
};

class Dynamical_variables{
private :
  nvector* resources;
  nvector* consumers;
public:
  Dynamical_variables(nvector* resources_, nvector* consumers_);
  nvector* get_resources() const;
  nvector* get_consumers() const;
};

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
class CRModel{
private:
  Metaparameters* metaparameters;
  ntensor* eq_vals; // we allow the possibility of multiple equilibria (hence vector of vector)
  Model_parameters* model_param;

public:
  CRModel(const foodmatrix&, Metaparameters&);
  CRModel(Model_parameters*);
  nvector equations_of_evolution(const Dynamical_variables&) const; // returns the value of the RHS of the equations of evolution
  nmatrix jacobian_at_equilibrium() const;
  ncvector eigenvalues_at_equilibrium() const;
  nmatrix jacobian(const Dynamical_variables&) const; // returns the jacobian for the given dynamical variables
  void save(std::ostream&) const; // outputs the model to the external file
  std::ostream& display(std::ostream&) const;
  bool energy_constraint() const;
  bool constraints_fulfilled(const Metaparameters& m) const;
  bool positive_parameters() const;
  bool dynamically_stable() const;
  void save_simulation() const;
  void save_jacobian_at_equilibrium(std::string) const;
  void write_time_evolution(const Dynamical_variables&, ntype) const;
  void write_time_evolution_from_equilibrium() const;
  void write_time_evolution_until_equilibrium(const Dynamical_variables &, ntype, ntype) const;
  void write_death_rates(std::string) const;
  nmatrix time_evolution(const Dynamical_variables&, ntype) const ;
  Dynamical_variables perturb_equilibrium() const;
  void perturb_parameters() const;
  void perturb_parameters(const ntype &) const;
  Extinction evolve_until_equilibrium(ntype) const;
  void save_new_equilibrium(const Extinction&) const;
};

nmatrix build_sigma(const Metaparameters&);
nvector build_resources(const Metaparameters&);
nvector build_consumers(const Metaparameters&);
nmatrix build_gamma(const foodmatrix&, const Metaparameters&);
nmatrix build_alpha(const Parameter_set*, Metaparameters&, const nvector&, unsigned int);
nmatrix build_tau(Parameter_set*, Metaparameters&, unsigned int);
nvector build_l(const Metaparameters& );

foodmatrix load_food_matrix(const Metaparameters&);
ntype norm(const nvector&);
ntype mean(const nmatrix &);
nmatrix random_uniform_matrix(const unsigned int&, const unsigned int&, const ntype&);
void rescale_mean(nmatrix&, const ntype&);

gammamode string_to_gamma_mode(std::string);
taumode string_to_tau_mode(std::string);
alphamode string_to_alpha_mode(std::string);

std::ostream& display_matrix_w_name(std::ostream&, std::string, const nmatrix&);
std::ostream& display_vector_w_name(std::ostream&, std::string, const nvector&);

bool non_neg_elements(const nmatrix&);
bool non_neg_elements(const nvector&);

void initialize_random_engine(const Metaparameters&);
void print_rand_number();

int ode_equations_of_evolution(double, const double[], double[], void*);

std::ostream& operator<<(std::ostream&, const Metaparameters&);
std::ostream& operator<<(std::ostream&, const nctype&);
std::ostream& operator<<(std::ostream&, const ncvector&);
std::ostream& operator<<(std::ostream&, const nvector&);
std::ostream& operator<<(std::ostream&, const nmatrix&);
std::ostream& operator<<(std::ostream&, const ntensor&);
std::ostream& operator<<(std::ostream&, const Parameter_set&);
std::ostream& operator<<(std::ostream&, const Model_parameters&);
std::ostream& operator<<(std::ostream&, const CRModel&);

// for a given set of metaparameters, computes the average extinction of the system
Extinction compute_average_extinction(Metaparameters*, const ntype &, unsigned int);
/* computes the critical delta (i.e. the delta of structural stability for which
   we have an average of 1.0 +/- accuracy number of extinctions
*/
double compute_critical_Delta(Metaparameters, ntype);
double function_av_extinct_solver(double, void*);
double average_number_of_extinctions(double , void*);
double estimate_delta_crit_from_interval(const nvector&, const nvector&);
int function_to_fit(const gsl_vector* , void* , gsl_vector*);

/* everything for the root solving here */
double solve_for_delta_with_fit(const gsl_vector* fit_parameters);

/* everything for the curve fitting here */
#define NUMBER_OF_FITTING_PARAMETERS 4;
double choice_of_fitting_function(double x, void* params);

#endif
