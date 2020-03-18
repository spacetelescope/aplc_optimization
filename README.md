# Installation

#### Clone aplc-optimization repository 

* `git clone https://github.com/spacetelescope/aplc-optimization.git`

#### Create a new conda environment and install dependencies

* `cd aplc_optimzation`
* `conda env create --file environment.yml`

#### Obtain a license for and download Gurobi 

* [Register](https://pages.gurobi.com/registration) an academic account with Gurobi.
* Download the latest version of the [Gurobi Optimizer](https://www.gurobi.com/downloads/gurobi-optimizer-eula/).
* Request an academic Gurobi [license](https://www.gurobi.com/downloads/end-user-license-agreement-academic/).
* In the directory where Gurobi Optimizer is installed, run ‘`grbgetkey`’ using the argument provided (e.g. `grbgetkey ae36ac20-16e6-acd2-f242-4da6e765fa0a`).
    * The ‘`grbgetkey`’ program will prompt you to store the license key on your machine, as well as validate your eligibility by confirming your academic domain (e.g., any ‘.edu’ address).















