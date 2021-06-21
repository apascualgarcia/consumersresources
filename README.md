# Consumers resources model
We here detail how the software we developed to numerically study consumers-resources models can be installed and used.
## Installation
### Pre-requisites
In order for every script written here, the following tools need to be installed on your system. All of these tools can be easily installed on your system e.g. with a package manager like **apt-get** on Ubuntu or **brew** on macOS. We will clearly write down the commands as an example for Ubuntu using apt-get.

* **A working C++ compiler** : the part of the program that produces data is written in C++. Hence you need a working C++ compiler. We recommend gcc although other compilers might work too. The appropriate command on Ubuntu using apt-get is
```
sudo apt-get install gcc
```
The two following C++ libraries are also required so that the program can be built.

* **GSL library** : This library (the GNU Scientific Library) can be installed through the following command on Ubuntu
```
sudo apt-get install libgsl-dev
```
Make sure that at least the 2.4 version is installed on your system.

* **Eigen 3 library** : This library can be installed through the following command on Ubuntu
```
sudo apt-get install eigen3
```
Make sure that at least the 3.3 version is installed on your system.

* **Cmake** : this command is crucial to write the makefile that will generate the executables. The appropriate command on Ubuntu using apt-get is
```
sudo apt-get install cmake
```
It is important that the version is at least 3.0. This can be checked using the command
```
cmake --version
```

* **Python 3** : the data analysis scripts are mostly written in Python 3. If Python 3 is not already installed on your system, you can use the following command to get it from the terminal :
```
sudo apt-get install python3
```

* **Python libraries** : Four Python libraries need to be installed on your Python 3 path in order to make the programs work. These are **os**, **numpy**, **matplotlib** and **glob**. If these are not already present on your system, you can install them e.g. with pip.
```
sudo apt-get install python3-pip
```
will install pip on your system. You can then type
```
pip3 install [libname]
```
where [libname] is the name of the wanted library.

* **ts** : this tool is needed in order to print timestamps in the console output when a script is running. It can be found e.g. in the **moreutils** package which can be installed with
```
sudo apt-get install moreutils
```
### Configuration
Once you made sure that the required tools are installed and you pulled the depository on your system, simply run
```
./configure
```
from the main folder to set up everything.

# Examples of different tasks that this package offers
## Generation of matrices that minimize a new energy
Let's say we would like to generate syntrophy matrices that minimize a new quadratic form  _E'(A,G,m)_ that has not been implemented in the code yet. The first step is to include _E'(A,G,m)_ in the code, i.e. implement it as a C++ function in the **src/CRModel_source/optimize_matrix.cc** file (it can e.g. be impedended at the end of the file). That C++ function, let's call it for instance `new_qf` should have the following structure:
```
ntype new_qf(const nmatrix& alpha, const nmatrix& gamma, void* additional_params){
   //definition of the new quadratic form.
}
```
The first argument must be the syntrophy matrix _A_, while the second argument is the competition matrix _G_, even if it ends up not being used in the quadratic form. Finally, the third argument can be a pointer to any additional parameters that would be needed, e.g. metaparameters.

Please note that the declaration of the function (so, without its complete definition) must also be added to the corresponding header file **include/Functions/optimize_matrix.h**.

Once this has been done, the remaining steps are simply to change the line 20 of the file **src/CRModel_targets/optimize_matrices.cc** so that the function called energy_function returns `new_qf` instead of whatever it was set to before. Finally, one needs to recompile everything (from the main directory) and build the needed dependancies:
```
make -C build
```
Finally, the matrices can be generated using the command:
```
build/optimize_matrices LOCATION_OF_THE_CONFIGURATION_FILE path_to_food_matrix=LOCATION_OF_THE_CONSUMPTION_MATRIX_LIST
```
A configuration file is always needed, in case for instance metaparameters are used in the energy function. The _A_ matrices corresponding to each _G_ in the provided matrix list will be stored in the location of the `path_to_syntrophy_matrix` variable, which can also be specified as an additional argument of the script.

### Scripts: what they do and how to use them
* **compute_critical_Delta_matrices** : this one is actually fairly simple. It computes the critical delta of a given set of matrices (given by default but which can be changed if needed) for every configuration of metaparameters specified. This means the total number of critical delta computed will be #matrices x #configuration. If you don't change the default matrix list and the default set of metaparameters, then you can simply run the script with the command
```
main_scripts/compute_delta_critical_matrices CORES
```
where **CORES** is an integer specifying the number of cores you want to allocate for this task. If you need to change the list of matrices or the set of configurations, you can do it manually in the script with a standard text editor, e.g. **nano** :
```
nano main_scripts/compute_delta_critical_matrices
```
Then you can specify the name of your set of matrices by changing the value of **MATRIX_LIST**. Be careful not to include the file extension. Note that the file extension must be '.in' and the file must be in the **config** folder, e.g.
```
MATRIX_LIST="other_matrix_list"
```
for a list of matrices that would be config/other_matrix_list.in.
Similarly, for an other set of metaparameters configurations
```
METAPARAMS_LIST="other_metaparams_list"
```
for the file config/other_metaparams_list.in

* **build/optimize_matrices**:
This script allows to create matrices which respect the needed energy conditions.
More specifically, it transforms each matrix of a set (the list of matrices
needs to be given as the path_to_food_matrix value of the input configuration file)
into a form which minimizes a given cost function energy_function (which can be)
changed on line 18). New energy functions can be added on the **optimize_matrix.cc** file
from the CRModel_source folder (the command ``./configure`` will need to be run again after the new energy form is added).

Typical usage (from main folder):
```
build/optimize_matrices PATH_TO_CONFIG_FILE path_to_food_matrix=PATH_OF_MATRIX_LIST
```
<<<<<<< HEAD
## Feasibility and dynamical stability
=======
# Feasibility and dynamical stability
>>>>>>> c4ff93926874c2d80c784c1ab553734f609a1039

The workload needed to compute the feasibility and dynamical stability data we are interested in is generally separated in two distinct steps. First, a text file containing the exact commands we would like to run is generated, either by hand or *as strongly advised* through the means of another script, and placed in the ``commands`` folder. The commands listed on the target text file, which we can call ``target.txt``, may then be executed with the command:

```
main_scripts/run_commands LOCATION_OF_THE_TARGET_FILE NUMBER_OF_CORES
```

The variable ``LOCATION_OF_THE_TARGET_FILE`` has a self-explicit name and would be for instance ```commands/target.txt```. The variable ```NUMBER_OF_CORES``` is an integer which specifies on how many cores the simulations should be run (**warning**: once the simulations are started, it is very painful to delete them by hand, it is therefore really important to think well before launching them).

## Feasibility

How feasibility works is explained in the main Thesis. For a matrix consumption matrix _G_ and a set of metaparemeters _M_, the most interesting metric is the percentage of feasible systems -- which we also simply call feasibility -- denoted _F_.
The file **find_feas_var_synt** located in the *gen_comm_files* folder allows to generate a command file which when executed will compute the feasibility of each configuration thrown as an input. The different variables that can be chosen are:
* `MATRIX_LIST`
* `ALPHA_MODE`
* `SAVENAME`
* `OPTIMAL_MATRIX_FOLDER`
* `ALPHA_VALS`
* `ADDITIONAL_MODIF`
* `CONFIGURATION_FILE`
* `COMMAND_FILE`
Note that all of them should be initialized in order for the command to work (even `ADDITIONAL_MODIF`. In case there are no additional modifications, just set it to the empty string `""`). `MATRIX_LIST` should be set to the name of the file which contains the list of all matrices for which _F_ . `ALPHA_MODE` is a list of strings which should be set to all the desired alpha modes. `SAVENAME` is the name under which the results should be saved. `OPTIMAL_MATRIX_FOLDER` is the location of the folder which contains the optimal syntrophy matrix, if needed. If no syntrophy matrix is used (e.g. if alpha=0) its value is not relevant. `ALPHA_VALS` is a string list which contains all the different alpha0 values for which we compute _F_. `CONFIGURATION_FILE` is the location of the configuration file which contains all the metaparameters _M_ (note that the value of alpha0 on that file will be overriden by `ALPHA_VALS`). Finally `COMMAND_FILE` is the name under which the generated command file should be saved in the */commands* folder.

### Usage example

To make it more clear, let's look at an example.
```
MATRIX_LIST="test_matrices"
ALPHA_MODE="optimal_matrix"
SAVENAME="run_001"
OPTIMAL_MATRIX_FOLDER="optimal_LRI_corrected_NR25_NS25"
ALPHA_VALS="0 0.5 1"
ADDITIONAL_MODIF="verbose-level=1"

CONFIGURATION_FILE="configuration"
COMMAND_FILE="test_feasibility"
```
Feasibility for all matrices in the **matrix_list/test_matrices.in** file will be computed. Alpha is here in "optimal" mode, and the optimal matrices are located in the **optimal_matrices/syntrophy/optimal_LRI_corrected_NR25_NS25** folder. Please note that it is important to put any new folder created in the right place and to use the right extensions, otherwise this script *will not* work without modifications. _F_ will be computed for three different values of alpha0, i.e. 0, 0.5 and 1. Additionally we ask that the operations should be "mildly" displayed on the terminal, i.e. a `verbose-level` equal to 1. The file **config/configuration.in** will be loaded as a configuration file when the commands will be executed. Finally, the command file generated with this script will be **commands/test_feasibility**.


## Dynamical stability

How dynamical stability works is also explained in the main Thesis. For a matrix consumption matrix _G_ and a set of metaparemeters _M_ *which we know is fully feasible* the main point of interest is the percentage of dynamically stable -- or dynamical stability -- systems for that configuration.

The file **find_dyn_stab** located in the *gen_comm_files* folder allows to generate a command file which when executed will compute the dynamical stability of each configuration thrown as an input. The different variables that can be chosen are:

* `MATRIX_LIST`
* `ALPHA_MODE`
* `SAVENAME`
* `OPTIMAL_MATRIX_FOLDER`
* `ALPHA_VALS`
* `ADDITIONAL_MODIF`
* `CONFIGURATION_FILE`
* `COMMAND_FILE`

The script and variables function exactly the same way as for the feasibility case.
<<<<<<< HEAD
=======

# Examples of different tasks
## Generation of matrices that minimize a new energy
Let's say we would like to generate syntrophy matrices that minimize a new quadratic form  _E'(A,G,m)_ that has not been implemented in the code yet. The first step is to include _E'(A,G,m)_ in the code, i.e. implement it as a C++ function in the **src/CRModel_source/optimize_matrix.cc** file (it can e.g. be impedended at the end of the file). That C++ function, let's call it for instance `new_qf` should have the following structure:
```
ntype new_qf(const nmatrix& alpha, const nmatrix& gamma, void* additional_params){
   //definition of the new quadratic form.
}
```
The first argument must be the syntrophy matrix _A_, while the second argument is the competition matrix _G_, even if it ends up not being used in the quadratic form. Finally, the third argument can be a pointer to any additional parameters that would be needed, e.g. metaparameters.

Please note that the declaration of the function (so, without its complete definition) must also be added to the corresponding header file **include/Functions/optimize_matrix.h**.

Once this has been done, the remaining steps are simply to change the line 20 of the file **src/CRModel_targets/optimize_matrices.cc** so that the function called energy_function returns `new_qf` instead of whatever it was set to before. Finally, one needs to recompile everything (from the main directory) and build the needed dependancies:
```
make -C build
```
Finally, the matrices can be generated using the command:
```
build/optimize_matrices LOCATION_OF_THE_CONFIGURATION_FILE path_to_food_matrix=LOCATION_OF_THE_CONSUMPTION_MATRIX_LIST
```
A configuration file is always needed, in case for instance metaparameters are used in the energy function. The _A_ matrices corresponding to each _G_ in the provided matrix list will be stored in the location of the `path_to_syntrophy_matrix` variable, which can also be specified as an additional argument of the script.
>>>>>>> c4ff93926874c2d80c784c1ab553734f609a1039
