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
enum alphamode{random_structure, no_release_when_eat, one_release};
/* eqmode tells you when you stop your time evolution algorithm */
enum eqmode{oneextinct, convergence};
/*  when using polynomial please specify the degree manually otherwise,
    there will be a runtime error */
enum fitmode{sigmoidal, polynomial, sigmoidal_erf};
enum stabilitymode{dynamical, structural};

enum systemstability {stable, marginal, unstable};
enum CRModelType{full, effective};

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
  statistics(const nvector&);
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


#endif
