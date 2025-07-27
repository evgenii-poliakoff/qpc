# QPC

This repository holds the python packages for the generalized light cone approach to the problem of out-of-equilibrium dynamics of quantum point contact  

## Getting started

Clone the repository:

```
git clone https://github.com/evgenii-poliakoff/qpc.git
cd qpc
```

Install the requirements (only once on a give host):

```
./scripts/install_requirements.sh
```

If the package was already built, then clean the project:

```
./scripts/clean.sh
```

Build the package:

```
./scripts/build.sh
```

Install the package:

```
./scripts/install.sh
```

Then cd to the notebooks subfolder and explore the jupyter notebooks:

```
cd notebooks
```

### Virtual python environment

#### When using qpc package in your python programs

The package is installed in the virtual python environment in the subfolder ./conda-env. Before using the package in your python programs you need to activate the virtual environment: 

```
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda activate ./conda-env
```

#### When using qpc package in VS Code jupyter notebooks

When running jupyter notebooks from VS Code, the virtual environment is activated by clicking on "Select Kernel", "Select Another Kernel ...", "Jupyter Kernel ...", "qpc-env (Python 3.XX.XX) conda-env/bin/python3". 

#### When using qpc package in JupyterLab notebooks

When running jupyter notebooks from JupyterLab web interface, the virtual environment is activated by clicking on "Python 3 (ipykernel)", "Select kernel", and choosing "qpc-env" in the dropdown list, then "Select".


## Working with jupyter notebooks

### Install jupyter lab

In order to install jupyter lab into the virtual environment of the package, execute the following commands:

```
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda activate ./conda-env
conda install jupyterlab
```

(run_jupyter_cmd)=
### Run jupyter notebook on the server as a command

Now you are ready to start the notebooks on the server side. Go to the directory where you have downloaded the notebooks (incliding this one, which you are currently viewing), and run:

```
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda activate ./conda-env
jupyter lab --no-browser --port 8889
```
    
This command will produce some output. You need the last line. It will look like:

    http://127.0.0.1:8889/?token=0aa445431a56428f9a71e2e7382b6a4924a08afcd3055971

Here the last part, something like 0aa445431a56428f9a71e2e7382b6a4924a08afcd3055971, is a token (password) you will need to enter
to login to your notebook though your web browser on a local computer.


(start_jupyter_service) = 
### Run jupyter notebook on the server as a background service (for lenghty computations)


If you want your notebook to run in background on the server even if you logout (e.g. to perform some lenghty Monte Carlo computation), then start the notebook with the command (on the SERVER):

```
eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
conda activate ./conda-env
jupyter lab --no-browser --port 8889 & disown
```
    
You can close the ssh connection - and the jupyter will continue to work.

(stop_jupyter) = 

### How to stop the jupyter notebook service

In order to stop the service (e.g. you want to change the environment variables for your python or other development environment), find the PID by running the command

```
pgreg jupyter
```
   
It will print the PID, e.g. 1234

Then terminate it by the command

```
kill 1234
```

### Login to jupyter notebook from local computer

Now we assume you can login via ssh to your server, 
and suppose your account looks like

john.smith@server

Then on your local machine you start ssh tunnel by running the following command on your LOCAL machine:

```
ssh -N -f -L localhost:8887:localhost:8889 john.smith@server
```
    
(reconnect_jupyter) = 
#### Reconnecting

The above command establishes connection to the server. You need to rerun it if the connection is lost. 


#### Web browser interface

Now open your local web browser, and type the address:

```
localhost:8887
```

then enter the token (see {ref}`run_jupyter_cmd`).
