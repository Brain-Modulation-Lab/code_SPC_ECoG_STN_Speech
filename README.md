# Sample code for: Spike-phase coupling of subthalamic neurons to posterior opercular cortex predicts speech sound accuracy

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/license/mit)
[![Generic badge](https://img.shields.io/badge/release-1.0.0-green.svg)](https://github.com/rutishauserlab/paper_SPC_ECoG_STN_Speech/releases/tag/v1.0.0)
[![Generic badge](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.10494534-orange.svg)](https://doi.org/10.5281/zenodo.12610957)
## Introduction
This repository contains the code and the preprocessed dataset for the preprint [Vissani et al 2024](https://doi.org/10.1101/2023.10.18.562969) "Spike-phase coupling of subthalamic neurons to posterior opercular cortex predicts speech sound accuracy". The full raw dataset is hosted in the [Data Archive BRAIN Initiative](https://dabi.loni.usc.edu/dsi/1U01NS098969) and is available upon request.

Abstract of the paper:
>Speech provides a rich context for understanding how cortical interactions with the basal ganglia contribute to unique human behaviors, but opportunities for direct intracranial recordings across cortical-basal ganglia networks are rare. We recorded electrocorticographic signals in the cortex synchronously with single units in the basal ganglia during awake neurosurgeries where subjects spoke syllable repetitions. We discovered that individual STN neurons have transient (200ms) spike-phase coupling (SPC) events with multiple cortical regions. The spike timing of STN neurons was coordinated with the phase of theta-alpha oscillations in the posterior supramarginal and superior temporal gyrus during speech planning and production. Speech sound errors occurred when this STN-cortical interaction was delayed. Our results suggest that the STN supports mechanisms of speech planning and auditory-sensorimotor integration during speech production that are required to achieve high fidelity of the phonological and articulatory representation of the target phoneme. These findings establish a framework for understanding cortical-basal ganglia interaction in other human behaviors, and additionally indicate that firing-rate based models are insufficient for explaining basal ganglia circuit behavior.

<p align="center">
  <img width="700" height="400" src="https://github.com/Brain-Modulation-Lab/code_SPC_ECoG_STN_Speech/blob/main/image/Figure1.png">
</p>

## Installation (Code)

This repository can be downloaded by entering the following commands:
```
git clone https://github.com/Brain-Modulation-Lab/code_SPC_ECoG_STN_Speech.git
cd code_SPC_ECoG_STN_Speech
```

## Minimal Dataset 

The minimum dataset `SPC_ECoG_STN_Speech` required to run the repository can be downloaded in [Zenodo](https://doi.org/10.5281/zenodo.12610957):
* To run `spc_demo.m` you need to copy the file `SPC_ECoG_STN_Speech/demos/intracranial-data-examples/intracranial_data.mat` in `code_SPC_ECoG_STN_Speech/demos/intracranial-data-examples`.
* To run the figures of the manuscript you need to copy the content of the folder `SPC_ECoG_STN_Speech/data` in `code_SPC_ECoG_STN_Speech/data`.
  
## External dependencies

The code depends on these repositories:

* [fieldtrip](https://www.fieldtriptoolbox.org/): toolbox to analyze electrophysiological data
* [bml](https://github.com/Brain-Modulation-Lab/bml): fieldtrip wrapper developed by the BrainModulation Lab.

You need to manually download and include (only the main folder!) them in your MATLAB dependencies.
After that just run these commands in MATLAB to manage dependencies:
```
bml_defaults
ft_defaults
```

## MATLAB Analysis

* DEMO Spike-phase coupling: the script `spc_demo.m` is designed to illustrate the computation of the time-resolved spike-phase coupling between 2 exemplary ECoG channels and 2 exemplary neurons contained in intracranial_data.mat.
* MAIN ANALYSIS: `Figure_02.m, Figure_03.m, Figure_04.m and Figure_05.m` are designed to analyze the minimal dataset and to reproduce  figures & metrics noted in the manuscript.
  - `Figure_02.m`: This figure describes the properties of the spike-phase coupling interaction between Subthalamic Nucleus and Cortex during speech production. To reduce computational effort, please make sure that the folder permutation_avgmaps is in your data folder. 
  - `Figure_03.m`: This figure describes the spatial distribution of the spike-phase coupling on the Cortex.
  - `Figure_04.m`: This figure describes the spatial distribution of the spike-phase coupling on the Subthalamic Nucleus.
  - `Figure_05.m`: This figure illustrates the error analysis contained in the manuscript. To reduce computational effort, please make sure that the folder permutation_avgmaps is in your data folder. 

If you encounter any problems, please report them as issues in the repository or send an [email](mailto:mvissani@mgh.harvard.edu).
This repository has been tested successfully in MATLAB versions 2022a and 2023a.

## Contributors
* [Matteo Vissani](mailto:mvissani@mgh.harvard.edu)

>Citation: [INSERT DOI HERE]


## Funding
This work was funded by the National Institute of Health (BRAIN Initiative), through grants U01NS098969, U01NS117836 and R01NS110424.

# License
This project is covered under the **MIT License**.
