#ifndef CUSTOM_TYPES_H
#define CUSTOM_TYPES_H

#include<iostream>
#include<complex>
#include<vector>


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
enum eqmode{oneextinct, convergence};
/*  when using polynomial please specify the degree manually otherwise,
    there will be a runtime error */
enum fitmode{sigmoidal, polynomial};
enum stabilitymode{dynamical, structural};

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


#endif
