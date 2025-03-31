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
Parameters that define the generalized pulse
- Pd: vector containing the indices of all the data carriers that use the generalized pulse (all the indices in $\mathcal{D}^\textrm{h}$)
- Pd_todas: vector containing the indices of all the data carriers (all the indices in $\mathcal{D}$)
- Paic: cell array with one entry per index in Pd. Each entry contains the indices of the CC employed by the corresponding generalized pulse
- Ptk: cell array with one entry per index in Pd. Each entry contains the indices of the terms used in the harmonic design of the transition pulses
- NexTemp_a/NexTemp_d: Number of advanced/delayed terms the generalized pulse has ($M_{\textrm{a}}/M_{\textrm{d}}$)
- NexTemp: Number of advanced and delayed terms the generalized pulse has when it is defined with Hermitian symmetry ($M_{\textrm{a}}=M_{\textrm{d}}$)
- R_aic: Cell array with one entry per index in Pd. Each entry contains the optimal coefficients $\alpha_{k,i}$ for one generalized pulse
- R_tk: Cell array with one entry per index in Pd. Each entry contains the optimal coefficients $\zeta_{k,i}$ for one generalized pulse ($\xi_{k,i}$ if transition pulses are harmonically designed)
- cell_coef_lim: Cell array with one entry per index in Pd. Each entry contains the bounds used to limit the coefficients' magnitude in the optimization procedure.
