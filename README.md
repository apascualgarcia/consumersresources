# Consumer resources model
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
where [libname] is the name of the wanted libraries.

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
