#ifndef OUTPUT_WRITING_H
#define OUTPUT_WRITING_H

#include "../Classes/Metaparameters.h"
#include "../Classes/Custom_types.h"
#include "../Classes/ParameterSet.h"
#include "../Classes/Model_Parameters.h"
#include "../Classes/Consumer_Resource_Model.h"
#include <iostream>

std::ostream& display_matrix_w_name(std::ostream&, std::string, const nmatrix&);
std::ostream& display_vector_w_name(std::ostream&, std::string, const nvector&);

std::ostream& operator<<(std::ostream&, const Metaparameters&);
std::ostream& operator<<(std::ostream&, const nctype&);
std::ostream& operator<<(std::ostream&, const ncvector&);
std::ostream& operator<<(std::ostream&, const nvector&);
std::ostream& operator<<(std::ostream&, const nmatrix&);
std::ostream& operator<<(std::ostream&, const ntensor&);
std::ostream& operator<<(std::ostream&, const Parameter_set&);
std::ostream& operator<<(std::ostream&, const Model_parameters&);
std::ostream& operator<<(std::ostream&, const CRModel&);
std::ostream& operator<<(std::ostream&, const stability_metrics&);
std::ostream& operator<<(std::ostream&, const stabilitymode&);
std::ostream& operator<<(std::ostream&, const statistics&);
std::ostream& operator<<(std::ostream&, const interval&);
std::ostream& operator<<(std::ostream&, const fitmode &);

#endif
