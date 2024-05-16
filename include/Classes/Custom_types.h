#ifndef CUSTOM_TYPES_H
#define CUSTOM_TYPES_H

#include<iostream>
#include<complex>
#include<vector>
#include <gsl/gsl_vector.h>
#include <fstream>

typedef long double ntype;

typedef std::complex<ntype> nctype;
typedef std::vector<ntype> nvector;
typedef std::vector<nctype> ncvector;
typedef std::vector<nvector> nmatrix;
typedef std::vector<nmatrix> ntensor;
typedef nmatrix foodmatrix;

typedef int(*func_equ_evol)(double, const double[], double[], void *);

// tau0 : tau = 0; taualpha : tau = alpha
enum taumode{tau0,taualpha};
enum gammamode{random_val, nested, antinested};
// when loading an external matrix, please use optimal_matrix mode
enum alphamode{fully_connected, random_structure, no_release_when_eat, one_release, optimal_matrix, binary_metamatrix, user_input,
// alphamodes with v2 suffix refer to the generation procedure proposed by APG in Apr 2024
fully_connected_v2,random_structure_v2, optimal_matrix_v2};
/* eqmode tells you when you stop your time evolution algorithm */
enum eqmode{oneextinct, convergence};
/* alpha value: input value or overriden by critical value */
enum alphavalue{input, critical};

/*
perturbmode tells which type of perturbation is chosen :
    - perturb_l : l_mu is perturbed by an amount delta
    - remove_res : Delta*NR resources are removed i.e. (Delta*NR) l_mu's are set to zero
*/
enum perturbmode{perturb_l, remove_l};


/*  when using polynomial please specify the degree manually otherwise,
    there will be a runtime error */
enum fitmode{sigmoidal, polynomial, sigmoidal_erf, linear};
enum stabilitymode{dynamical, structural};

enum systemstability {stable, marginal, unstable};
enum CRModelType{full, effective};

/* different ways of building the system, do we choose l? m? */
enum buildingmode{use_l, use_m};

/*  three different types to run the MC algorithm :
      - A_only : alpha only is changed and is unconstrained
      - both_modified : BOTH alpha and gamma are modified by the algorithm. alpha is modified without constraints
                        but gamma must have full rank AND a given connectance   */
enum MCmode{A_only, both_modified};


/*  writemode is used in the general time evolution of the system. It tells you whether you should, and if so Where
    write the time evolution of the system */
struct writemode{
  bool write;
  std::ostream& write_path;

  writemode(bool, std::ostream&);
  writemode();
};

struct statistics{
  ntype mean_;
  ntype std_deviation_;
  ntype median_;
  statistics(const nvector&, const unsigned int ddof=0);
  statistics(const statistics&);
  statistics();
};

struct nvector_statistics{
  nvector means;
  nvector std_deviations;
  nvector medians;
};

struct stability_metrics{
  statistics resilience;
  statistics angle_between_equilibria;
  statistics distance_between_equilibria;
  statistics extinctions;
};

struct interval{
  ntype begin;
  ntype end;

  interval(ntype a, ntype b){
    if(a < b){
      this->begin = a;
      this->end = b;
    }
    else{
      this->begin = b;
      this->end = a;
    }
  };
};

struct Extinction{
  ntype t_eq;
  ntype extinct;
  nvector new_Req;
  nvector new_Seq;
  nvector old_Req;
  nvector old_Seq;
};

struct Extinction_statistics{
  statistics t_eq;
  statistics extinct;
  nvector_statistics new_Req;
  nvector_statistics new_Seq;
};

struct stability{
  ntype stable;
  ntype marginally_stable;
  ntype unstable;
};

struct error{
  std::string message;
  unsigned int category;
  error(std::string a, unsigned int cat=0):message(a), category(cat){};
  void handle();
};

struct MCS_Running_Parameters{
    /* constant values determined at the initialization */
    const unsigned int max_steps = 2000000;
    const unsigned int max_fails = 300;
    const unsigned int annealing_freq = 1000;
    const double annealing_const = 0.99;
    const unsigned int display_stride = 10000;

    /* constant values used for the determination of energy convergence */
    const unsigned int convergence_Naverage = 50;
    const unsigned int required_convergence = 100;
    const double eps = 1e-2;

    /* variables changing while the algorithm is running */
    nvector last_changed_elements;
    ntype current_temperature = 5;
    unsigned int current_fails = 0;
    unsigned int current_step = 0;
    unsigned int current_convergence = 0;
    bool changed_matrix = false;


    /* functions */
    void display(std::ostream&);
};




#endif
