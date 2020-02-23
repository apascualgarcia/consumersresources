#ifndef MATRICES_BUILDING_H
#define MATRICES_BUILDING_H

#include "../Classes/Metaparameters.h"
#include "../Classes/Custom_types.h"
#include "../Classes/ParameterSet.h"

nmatrix build_sigma(const Metaparameters&);
nvector build_resources(const Metaparameters&);
nvector build_consumers(const Metaparameters&);
nmatrix build_gamma(const foodmatrix&, const Metaparameters&);
nmatrix build_alpha(const Parameter_set*, Metaparameters&, const nvector&, unsigned int);
nmatrix build_tau(Parameter_set*, Metaparameters&, unsigned int);
nvector build_l(const Metaparameters& );
nvector build_m(const Metaparameters& ); 

nmatrix build_sigma_Butler(const Metaparameters&);

#endif
