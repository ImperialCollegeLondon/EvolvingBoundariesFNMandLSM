# A novel formulation for the explicit discretisation of evolving boundaries with application to topology optimisation 

[![DOI](https://badgen.net/badge/DOI/10.1016%2Fj.cma.2020.113077/cyan)](https://doi.org/10.1016/j.cma.2020.113077) [![DOI](https://zenodo.org/badge/232163541.svg)](https://zenodo.org/badge/latestdoi/232163541)
---
## Authors  
Rui O. S. S. da Costa ([r.costa18@imperial.ac.uk](mailto:r.costa18@imperial.ac.uk))  
Silvestre T. Pinho ([silvestre.pinho@imperial.ac.uk](mailto:silvestre.pinho@imperial.ac.uk))

Department of Aeronautics  
Imperial College London  
South Kensington Campus  
SW7 2AZ, London  
United Kingdom

<br>

## Description
In the folder **_Code_** you will be able to find the implementation  of both, baseline and proposed, methods from the paper. All results files and images of the paper can be found in the folder **_Results_**.

The proposed method is capable of explicit modelling evolving boundaries, by coupling the _floating node method_ and the _level set method_, in an accurate and efficient way. You can find more details on the theory, implementation and results in the [manuscript](https://authors.elsevier.com/a/1b8WWAQEIt0mR).

<br>

## How to use?
This code is made available under the GPL license enclosed with the software. Further to the legal implications of the license, to use this code, please reference:  
* [the publication](https://doi.org/10.1016/j.cma.2020.113077) [![DOI](https://badgen.net/badge/DOI/10.1016%2Fj.cma.2020.113077/cyan)](https://doi.org/10.1016/j.cma.2020.113077)

* [the correct software release](https://zenodo.org/badge/latestdoi/232163541) [![DOI](https://zenodo.org/badge/232163541.svg)](https://zenodo.org/badge/latestdoi/232163541)


### Opening the files
Open Matlab (tested on versions 2018 and above) and navigate to the folder containing the source files.

### Understanding the file names
There are essentially three file names: files with the prefix "fem_" belong to the FEM analysis; files with the prefix "ls_" belong to the Level-set analysis; files with no prefix are shared by both. 

### Running the code
You can run the code in two ways:
1. _runscript.m_

Select the plotting options, output options and the desired test case and run the code directly from this file.

<br>

2. _main.m_

Define the plotting and output options (s_plot) in the Matlab command window with the following commands:
```
s_plot.mesh = 0; % initial mesh
s_plot.ls = 1; % ls field
s_plot.fnm = 0; % mesh partition
s_plot.dens = 1; % meshed domain
s_plot.dof = 0; % displacement solution
s_plot.str = 0; % stress and strain fields
s_plot.vel = 0; % velocity field
s_plot.save = 1; % save .mat file for each iteration
```
in which the values can be changed according to preference.

Run the code with the following command:
```
% Cantilever beam example
[obj, vol, file] = main(1, 2, 1, 100, 50, 7, 5, 0.5, 1, 10, s_plot);
```
in wich the input data can be modified accordingly.
