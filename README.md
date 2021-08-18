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
We recommend to install pip3, which will make it easier to install the required packages:
```
sudo python3 -m pip install pip
```

* **Python libraries** : Some Python libraries need to be installed on your Python 3 path in order to make the programs work. To avoid messing up your current setup, we recommend the usage of a virtual environment. The file requirements.txt provides a list of all needed python3 packages, which have to be downloaded within the virtual environment. We provide here a working example on how to set this up using the command **venv** to create a virtual environment named _consumers_resources_:
```
python3 -m venv consumers_resources
```
We then activate the virtual environment so we can work in it without disturbing the rest of our python3 installation:

```
source bin/activate/consumers_resources
```
This new virtual environment does not have any package installed yet so we install all the ones required for this project:
```
pip3 install -r requirements.txt
```

We are now ready to work with any data analysis command in the project. When you have finished working on this project, do not forget to deactivate the virtual environment:

```
deactivate
```

Whenever you start working with this project, don't forget to activate your virtual environment!

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
## General structure of the package
The software we developed is very comprehensive and it is easy to get lost in it. Scripts are provided to perform in an explicit way the most important tasks required by the project. The user has the possibility (and must) provide the number of CPU cores which are reserved for the specific task at hand. The general procedure is the following: the user first generates a "command" file which is then executed. One can generate a command file by hand (simply by stacking on a text file all the commands that should be executed) but this is not recommended. An easier way is to use the provided script configuration files which can be found in the **./gen_comm_files** folder.

It is surely clearer to explain the procedure with an example. Let us say we would like to create a script which will generate the optimal syntrophy matrix of a given list of consumption matrices. The file **./gen_comm_files/optimize_matrices** may be customized according to the comments on the file. The script is then generated by running the command:
```
./gen_comm_files/optimize_matrices
```
The newly created command files is then found in the **commands** folder under the name the user chose, e.g. "comm_to_run.txt".
The script can be executed with the following command:
```
./run_commands commands/comm_to_run.txt CORES
```
where _CORES_ is the number of cores (integer) which should be reserved for the script. In order to avoid cluttering the terminal, every script is run in the background. It uses the _nohup_ command and runs even when the user logs out (e.g when an ssh session is left) and can hence only be killed manually commands similar to _pkill_.

## Tasks that can be accomplished with this package
### Generation of optimized matrices
This package offers the possibility of finding the ecological network _(A, G)_ which minimizes a given energy function _E(A,G,m)_ (_m_ are simply metaparameters). Using a Monte Carlo algorithm, one can find the syntrophy matrix _A_ which minimizes _E(A,G,m)_ for a given consumption matrix _G_ and a set of metaparameters _m_. One also has the possibility to find both _A_ and _G_ for a given _m_. This can be set up using the provided script in **./gen_comm_files/optimize_matrices**. As explained above, the scripts can be run with the commands:
```
./gen_comm_files/optimize_matrices
./run_commands commands/optimize_matrices.txt CORES
```

The energy function _E(A,G,m)_ is by default the one from the main text but may be easily modified if needed. Let us say we would like to generate syntrophy matrices that minimize a new quadratic form  _E'(A,G,m)_ that has not been implemented in the code yet. The first step is to include _E'(A,G,m)_ in the code, i.e. write it as a C++ function in the **src/CRModel_source/optimize_matrix.cc** file (it can e.g. be appended at the end of the file). That C++ function, let's call it for instance `new_qf` must have the following structure:
```
ntype new_qf(const nmatrix& alpha, const nmatrix& gamma, void* additional_params){
   //definition of the new quadratic form.
   //energy form must be returned at the end as a ntype number (aka long double)
}
```
The first argument must be the syntrophy matrix _A_, while the second argument is the competition matrix _G_, even if it ends up not being used in the quadratic form. Finally, the third argument is a pointer to any additional parameters that would be needed, e.g. metaparameters. It must be present in the definition of the function even if no extra parameters are needed. In that case, the *null* pointer may simply be passed.

Please note that the declaration of the function (so, without its complete definition) must also be added to the corresponding header file **include/Functions/optimize_matrix.h**.

Once this has been done, the remaining steps are simply to change the line 20 of the file **src/CRModel_targets/optimize_matrices.cc** so that the function called energy_function returns `new_qf` instead of whatever it was set to before. Finally, one needs to recompile everything (from the main directory) and build the needed dependancies:
```
make -C build
```
The new energy form has then been added and the commands above can be run to generate the new optimized matrices.


### Feasibility and dynamical stability
One script is used to assess the feasibility and local dynamical stability of a given ecological network. We first explain which metrics are used and how they are computed and then we discuss about how this can be performed.

To assess feasibility/dynamical stability of a given ecological network, we estimate the percentage of feasible/dynamically stable within a rectangular cuboid ![cuboid](http://www.sciweavers.org/tex2img.php?eq=%28%5Cgamma_0%2C%20S_0%2C%20%5Calpha_0%29%20%5Cin%20%5B0.01%2C%201%5D%20%5Ctimes%20%5B0.01%2C%201%5D%20%5Ctimes%20%5B0%2C%200.02%5D&bc=White&fc=Black&im=jpg&fs=12&ff=modern&edit=0). For each ecological network, the procedure is the following:
1. Choose a random ![point](http://www.sciweavers.org/tex2img.php?eq=%28%5Cgamma_0%2C%20S_0%2C%20%5Calpha_0%29&bc=White&fc=Black&im=jpg&fs=12&ff=modern&edit=0) point in the metaparameters-cuboid.
1. Generate a parameters set from these metaparameters.
1. Check whether that parameters set is feasible or not.
1. If the parameters set is feasible, check if it is dynamically stable.
1. If it is dynamically stable, compute its rate of return to equilibrium and its effective competition.


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



## Structural stability

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
