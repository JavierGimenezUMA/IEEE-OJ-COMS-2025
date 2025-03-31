# DATASETS

In this folder you can find the datasets that are employed by the MATLAB scripts in '/scripts/'. These are organized accordingly to the Fig. or Table of our paper they belong to.
For convenience, we detail here how the most important variables relate to the parameters used in our manuscript. To that end, we have arranged them in two sets:

### OFDM_param
Parameters that define the OFDM signal.
- N: size of the DFT ($N$)
- Ngi: size of the guard interval ($N_{\textrm{GI}}$)
- L: symbol period ($N_{\textrm{s}}$)
- b: size of each of the tapered transitions of the shaping pulse ($\beta$)
- g: time-domain samples of the shaping pulse ($g(n)$)

Parameters that appear just in the scripts
- Ms: frequency-domain oversampling factor
- us, ue: time-domain pulses of length $\beta$ that are non-zero just in the starting and ending region of the pulse, respectively

### GenPulse_param
Parameters that define the proposed pulse
- Pd: vector containing the indices of all the data carriers that use the proposed pulse (all the indices in $\mathcal{D}^\textrm{h}$)
- Pd_todas: vector containing the indices of all the data carriers (all the indices in $\mathcal{D}$)
- Paic: cell array with one entry per index in Pd. Each entry contains the indices of the CC employed by the corresponding proposed pulse
- Ptk: cell array with one entry per index in Pd. Each entry contains the indices of the terms used in the harmonic design of the transition pulses
- NexTemp_a/NexTemp_d: Number of advanced/delayed terms the proposed pulse has ($M_{\textrm{a}}/M_{\textrm{d}}$)
- NexTemp: Number of advanced and delayed terms the proposed pulse has when it is defined with Hermitian symmetry ($M_{\textrm{a}}=M_{\textrm{d}}$)
- R_aic: Cell array with one entry per index in Pd. Each entry contains the optimal coefficients $\alpha_{k,i}$ for one proposed pulse
- R_tk: Cell array with one entry per index in Pd. Each entry contains the optimal coefficients $\zeta_{k,i}$ for one proposed pulse ($\xi_{k,i}$ if transition pulses are harmonically designed)
- cell_coef_lim: Cell array with one entry per index in Pd. Each entry contains the bounds used to limit the coefficients' magnitude in the optimization procedure.

The way the coefficients' bounds are arranged in the entries of cell_coef_lim depends on whether the proposed pulse has Hermitian symmetry or not. When it is designed with Hermitian symmetry, it is:

![alt text](https://github.com/JavierGimenezUMA/IEEE-OJ-COMS-2025/blob/main/data/Bounds_vector_Hermitian_pulse.PNG "Content of one of the entries of cell_coef_lim when the pulse has Hermitian symmetry")

where 

![alt text](https://github.com/JavierGimenezUMA/IEEE-OJ-COMS-2025/blob/main/data/Epsilon_vector.PNG "Definition for the epsilon vector of bounds")

and

![alt text](https://github.com/JavierGimenezUMA/IEEE-OJ-COMS-2025/blob/main/data/Delta_vector.PNG "Definition for the delta vector of bounds")

depending on whether the transition pulses are regularly designed (first definition) or harmonically designed (second definition)

When the proposed pulse is not designed with Hermitian symmetry, the boundsa are arranged as follows
