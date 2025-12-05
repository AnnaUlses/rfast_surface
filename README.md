# rfast_flux_grey
Modified version of the open source radiative transfer code created by Tyler Robinson and Arnaud Salvador: [Robinson & Salvador 2023](https://doi.org/10.3847/PSJ/acac9a)  <br/>

Code has been modified to run forward models and retrievals using a grey surface albedo. Retrievals are configured for flux ratio (Fp/Fs) spectra. 

**Description of changeable files** <br/>
```rfast_inputs.scr``` is the inputs file with user-chosen parameters <br/>
```rfast_user_model.py``` script with models for the forward model, surface albedo parametrisation can be modified through ```surfalb``` function <br/>
```rfast_genspec_alb.py``` is the forward modeling script <br/>
```rfast_noise.py``` generates noisy spectra with user specified noise model <br/>
```rfast_retrieve_pp.py``` main script for MCMC analysis - outputs a .h5 file with the chain data <br/>
```rfast_analyze_pp.py``` first analysis script, generates corner plots, individual posteriors for land fractions, and fitted spectra <br/>
```spectra_gen_pp.py``` second analysis script for plotting confidence intervals of fits <br/>

**How to run** 
1. ```rfast_genspec_alb.py rfast_inputs.scr```
2. ```rfast_noise.py rfast_inputs.scr```
3. ```rfast_retrieve_pp.py rfast_inputs.scr```
4. ```rfast_analyze_pp.py rfast_inputs.scr```
5. ```spectra_gen_pp.py rfast_inputs.scr``` <br/>

This code relies on a co-located opacities file called ```hires_opacities```. This and installation instructions can be found in the main [rfast](https://github.com/hablabx/rfast) github. <br/>

**Publications that use this code** <br/>
Recent paper submitted to Astrobiology: [https://arxiv.org/pdf/2507.14771](url)

