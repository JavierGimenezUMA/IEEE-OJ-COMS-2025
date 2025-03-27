# Datasets and Scripts for "Spectral Shaping Method for OFDM Combining Time-shifted Active Interference Cancellation and Adaptive Symbol Transition"

**This repository contains datasets and MATLAB scripts to generate the curves presented in the figures of our paper. Only the results corresponding to the method we propose can be generated with these materials.**

📄 **"Spectral Shaping Method for OFDM Combining Time-shifted Active Interference Cancellation and Adaptive Symbol Transition"**, submitted to **IEEE Open Journal of the Communications Society**, 2025

📌 DOI: 

## Overview
This repository contains:
- 💾 **Datasets** (/data/): One dataset per curve is provided, which are grouped in folders according to the Fig. they belong to. Each dataset contains all the parameters required to generate one curve. That includes (please, refer to the paper for the symbol definitions):
  - OFDM signal parameters: $N$, $N_{gi}$, $\beta$, $L$, $g(n)$
  - Parameters of the proposed pulse: $\mathcal{D}$, $\mathcal{D}^h$, $\mathcal{C}$, $M_a$, $M_d$, $\boldsymbol{\epsilon_{k,i}}$, $\boldsymbol{\delta}_{k,i}$
  - The optimal set of coefficients to generate the pulse: $\boldsymbol{\gamma}_{k,i}$
- ⚙️ **Scripts** (/scripts/): There is one MATLAB script per curve, which are grouped in folders according to the Fig. they belong to. Each script loads one of the datasets in and generates the curve.

📜 These materials are provided for the sake of reproducibility of the results in our paper. Readers are encouraged to utilize them under the terms of the journal and this respository's license.

 ## Requirements
 To run the scripts, it is required:
 
 ✔️ MATLAB R2022b (or later)


 ## Usage Instructions
 1. The repository can be cloned or downloaded as a zip.
 2. Open MATLAB and navigate to main/
 3. Add /scripts/ and /data/ to the workspace's path
 4. Run the script
 5. The corresponding curve(s) should be plotted
 

