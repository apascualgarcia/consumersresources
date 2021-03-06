#ifndef INPUT_READING_H
#define INPUT_READING_H

#include "../Classes/Custom_types.h"
#include <string>

gammamode string_to_gamma_mode(std::string);
taumode string_to_tau_mode(std::string);
alphamode string_to_alpha_mode(std::string);
eqmode string_to_eq_mode(std::string);
buildingmode string_to_building_mode(std::string);
MCmode string_to_mcmode(std::string);
perturbmode string_to_perturbmode(std::string);
alphavalue string_to_alpha_value(std::string);


#endif
