# rfast_modified
Modified version of the open source radiative transfer code from [rfast](https://github.com/hablabx/rfast) (**not my code**)
Code has been modified to run forward models and retrievals using a user-defined albedo model (defined in user_models). Currently set-up for a two parameter parametrisation, further parameters will require small modifications in all scripts.

**Description of changeable files**
'''rfast_inputs.scr''' is the inputs file with user-chosen parameters <br/>
'''rfast_user_model.py''' script with models for the forward model, surface albedo parametrisation can be modified through '''surfalb''' <br/>
'''rfast_genspec_alb.py''' is the forward modeling script <br/>
'''rfast_noise.py''' generates noisy spectra with user specified noise model <br/>
'''rfast_retrieve_pp.py''' main script for MCMC analysis - outputs a .h5 file with the chain data <br/>
'''rfast_analyze_pp.py''' first analysis script, generates corner plots, individual posteriors for land fractions, and fitted spectra <br/>
'''spectra_gen_pp.py''' second analysis script for plotting confidence intervals of fits <br/>

**How to run**
1. '''rfast_genspec_alb.py rfast_inputs.scr'''
2. '''rfast_noise.py rfast_inputs.scr'''
3. '''rfast_retrieve_pp.py rfast_inputs.scr'''
4. '''rfast_analyze_pp.py rfast_inputs.scr'''
5. '''spectra_gen_pp.py rfast_inputs.scr'''
