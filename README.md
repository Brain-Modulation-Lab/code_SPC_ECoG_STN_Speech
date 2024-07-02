# Sample code for: Spike-phase coupling of subthalamic neurons to posterior opercular cortex predicts speech sound accuracy

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/license/mit)
[![Generic badge](https://img.shields.io/badge/release-1.0.0-green.svg)](https://github.com/rutishauserlab/paper_SPC_ECoG_STN_Speech/releases/tag/v1.0.0)
[![Generic badge](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.10494534-orange.svg)](https://doi.org/10.5281/zenodo.12610957)
## Introduction
This repository contains the code and the preprocessed dataset for the preprint [Vissani et al 2024](https://doi.org/10.1101/2023.10.18.562969) "Spike-phase coupling of subthalamic neurons to posterior opercular cortex predicts speech sound accuracy". The full raw dataset is hosted in the [Data Archive BRAIN Initiative](https://dabi.loni.usc.edu/dsi/1U01NS098969) and is available upon request.
This repository can be used to calculate time-resolved spike-phase coupling in event-based paradigms adapting the width-specific window approach proposed by [Fischer et al](https://elifesciences.org/articles/51956).

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
* To run the figures of the manuscript you need to copy the content of the folder `SPC_ECoG_STN_Speech/data` in `code_SPC_ECoG_STN_Speech/data`:
  - `DB_main_analysis.mat` contains the results of the spike-phase coupling computation in the main analysis. It contains the description of the ECoG and neuron pairs and the features of the identified spike-phase coupling events.
  - `DB_error_analysis.mat` contains the results of the spike-phase coupling computation in the error analysis.
  - `DISTAL_atlas.mat` contains the meshes for the visualization of the Subthalamic Nucleus as depicted by the [DISTAL](https://doi.org/10.1016/j.neuroimage.2017.05.015) atlas.
  - `Cortex_MNI.mat` contains the meshes for the visualization of the Cortex. ROIs are parcellated using the [Destrieux](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation) atlas.
  - `permutation_avgmaps` is a folder that contains pre-computed permutation tests. Using this set of permutations can accelerate the run of the code.
  - `tSPC_density_ECoG.txt` contains the spatial density of the spike-phase coupling on the Cortex. Alternatively, the `.node` version can be used to plot the results using [SurfIce](https://www.nitrc.org/projects/surfice/).
  - `tSPC_density_STN.csv` contains the spatial density of the spike-phase coupling on the Subthalamic Nucleus.

## External dependencies

The code depends on these repositories:

* [fieldtrip](https://www.fieldtriptoolbox.org/): toolbox to analyze electrophysiological data
* [bml](https://github.com/Brain-Modulation-Lab/bml): fieldtrip wrapper developed by the [BrainModulation Lab](https://www.brainmodulationlab.org/).

You need to manually download and include (only the main folder!) them in your MATLAB dependencies.
After that just run these commands in MATLAB to manage dependencies:
```
bml_defaults
ft_defaults
```
The external folder in the repo contains other libraries:
* [Circstat](https://github.com/circstat/circstat-matlab): toolbox to analyze circular data
* [RainCloud](https://github.com/RainCloudPlots/RainCloudPlots): toolbox to visualize data distributions.
* [PERMUTOOLS](https://github.com/mickcrosse/PERMUTOOLS/tree/master): toolbox to implement permutation-based statistics.

## MATLAB Analysis

* DEMO SPIKE-PHASE COUPLING: the script `spc_demo.m` is designed to illustrate the computation of the time-resolved spike-phase coupling between 2 exemplary ECoG channels and 2 exemplary neurons contained in intracranial_data.mat.
  - The script `calc_spike_PLV_all.m` is the core function that calculates the spike-phase coupling metric. The pipeline inherits the approach proposed by [Fischer et al](https://elifesciences.org/articles/51956).
  - The script `set_configs.m` allows to toggle different parameters for the spike-phase coupling computation. Some of these parameters include:
    - `cfg.plv.NUM_PERMS = 500 (recommended)`: number of permutation maps to normalize the spike-phase coupling metric (500 recommended but the computational effort increases a lot, use 80 to have an initial guess)
    - `cfg.locked = true (recommended)`: if true it removes the event-locked component before the computation (true is recommended)
    - `cfg.powertrim = true (recommended)`: if true it removes the most-extreme power events, i.e., pauses (< 10th percentile) and bursts (> 90th percentile) (true is recommended)
    - `cfg.plv.MIN_TRIALS = 10 (recommended at least)`: number of minimum set of required trials

* MAIN ANALYSIS: `Figure_02.m, Figure_03.m, Figure_04.m and Figure_05.m` are designed to analyze the minimal dataset and to reproduce  figures & metrics noted in the manuscript.
  - `Figure_02.m`: This figure describes the properties of the spike-phase coupling interaction between Subthalamic Nucleus and Cortex during speech production. To reduce computational effort, please make sure that the folder permutation_avgmaps is in your data folder. (~20 min)
  - `Figure_03.m`: This figure describes the spatial distribution of the spike-phase coupling on the Cortex. (~10 min)
  - `Figure_04.m`: This figure describes the spatial distribution of the spike-phase coupling on the Subthalamic Nucleus. (~10 min)
  - `Figure_05.m`: This figure illustrates the error analysis contained in the manuscript. To reduce computational effort, please make sure that the folder permutation_avgmaps is in your data folder. (~10 min)

If you encounter any problems, please report them as issues in the repository or send an [email](mailto:mvissani@mgh.harvard.edu).
This repository has been tested successfully in MATLAB versions 2022a and 2023a on MacOS.

## Contributors
* [Matteo Vissani](mailto:mvissani@mgh.harvard.edu)

>Citation: [INSERT DOI HERE]

## Full Dataset request
The full raw dataset is hosted in the [Data Archive BRAIN Initiative](https://dabi.loni.usc.edu/dsi/1U01NS098969) and is available upon request.

## Funding
This work was funded by the National Institute of Health (BRAIN Initiative), through grants U01NS098969, U01NS117836 and R01NS110424.

# License
This project is covered under the **MIT License**.
