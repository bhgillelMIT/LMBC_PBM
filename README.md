## Overview

  This repo contains a population balance model (PBM) being developed at MIT to simulate a liquid metal bubble column (LMBC) reactor. The model  resolves the performance of the reactor,
  accounting for effects that degrade performance. It considers a population of bubbles and accounts for their rise, expansion, coalescence, breakage, heat transfer, a reaction, and 
  liquid circulation. The model is intended to provide an intermediate level of fidelity, balancing computational cost with accuracy.

  While the model is setup for our specific case, it has been written more generally to solve population balanace models for other cases. However, it will require familiarity with MATLAB 
  to implement other applications.


## How to Use 

  The model was developed in MATLAB R2025b. The repo is self-contained, only referencing folders and files within it, or built into this version of MATLAB and a few of its toolboxes. 
  These include the Statistics and Machine Learning toolbox, and the Curve Fitting toolbox. The code should run so long as all the dependencies are installed. 

  The main inputs for the model are specified in the function PBM_inputs(), which outputs a structure containing them. They can either be changed within this file, or the file can be 
  called from a separate function allowing them to be modified for specific cases. Additional, more detailed inputs, are contained in the main function, PBM_v3(). Explanations for each
  input are included in comments beside them. Care is required when changing these.
