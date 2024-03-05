# Functions `shiny` app

## Description

This repository contains the R functions for the shiny-app CompARE (https://www.grbio.eu/compare/). CompARE is a web-platform inspired to provide help on issues relating to trials with composite endpoints.

[**CompARE**](https://www.grbio.eu/compare/) comprises two different
Shiny apps: one devoted to time-to-event endpoints, the other to binary
endpoints.

  - **CompARE for Time-to-event endpoints**: [Time-to-event
    app](https://www.grbio.eu/compare/CompARETimeToEvent/).
  - **CompARE for Binary endpoints**: [Binary
    app](https://www.grbio.eu/compare/CompAREBinary/).


Their user-friendly interface allows one to input the main parameters
included in the trial -such as the treatment effect on the components of
the composite endpoint, and the frequencies of occurrence- and the app
provides sample size calculations among others.

## Getting Started

If you are a newcomer, we recommend starting with the tutorial
vignettes. These vignettes provide an introduction to CompARE:

  - **Time-to-event Tutorial**: Guide document of CompARE for
    Time-to-event endpoints [Time-to-event
    Tutorial](https://www.grbio.eu/compare/CompARETimeToEvent/help_Tutorial.html).
  - **Binary Tutorial**: Guide document of CompARE for Binary endpoints
    [Binary
    Tutorial](https://www.grbio.eu/compare/CompAREBinary/Help-Tutorial.html).


## Contents of this repository

- _**Functions_Binary**_ contains the functions referred to binary composite endpoints.
- _**Functions_T2E**_ contains the functions referred to time-to-event composite endpoints.
- _**BCG_CaseStudies**_ contains the code to reproduce the figures and tables shown in the app based on parameters from real studies.

Visit https://www.grbio.eu/compare/ to learn more.

## Publications and related repositories

The repository https://github.com/MartaBofillRoig/CompARE  contains the source files of the papers:

- Selection of composite binary endpoints in clinical trials. Biom J. 2018 Mar; 60(2):246-261. (https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201600229).
	
- Sample size derivation for composite binary endpoints. Stat in Med. DOI: 10.1002/sim.8092. (https://arxiv.org/abs/1807.01305).
