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

## Tree structure of the package

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
## Feasibility and dynamical stability

The workload needed to compute the feasibility and dynamical stability data we are interested in is generally separated in two distinct steps. First, a text file containing the exact commands we would like to run is generated, either by hand or *as strongly advised* through the means of another script, and placed in the ``commands`` folder. The commands listed on the target text file, which we can call ``target.txt``, may then be executed with the command:

```main_scripts/run_commands LOCATION_OF_THE_TARGET_FILE NUMBER_OF_CORES```

The variable ``LOCATION_OF_THE_TARGET_FILE`` has a self-explicit name and would be for instance ```commands/target.txt```. The variable ```NUMBER_OF_CORES``` is an integer which specifies on how many cores the simulations should be run (**warning**: once the simulations are started, it is very painful to delete them by hand, it is therefore really important to think well before launching them).

### Feasiblity 





