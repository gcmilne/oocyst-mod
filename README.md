## Estimating the longevity of *Toxoplasma gondii* oocyst-specific antibodies

### Description of model
This project presents a serocatalytic model that tracks the change in two *T. gondii*-specific antibodies in human sera over age. The first, non-stage-specific IgG, is acquired following exposure to *T. gondii*, irrespective of the infecting stage. The second antibody, anti-*T. gondii* embyrogenesis-related protein (TgERP), is specific to the oocyst stage and hence is only present in oocyst-infected individuals. Individuals seroconvert from seronegative to seropositive at a rate determined by the age-specific force of infection (the functional form of which is determined by the values of three fitted parameters). For oocyst-specific antibodies, the force of infection is multiplied by the proportion of infections that are oocyst-acquired (as opposed to infections acquired via other parasite stages) to derive the oocyst-specific force of infection. The model presented here estimates the duration of anti-TgERP IgG antibodies in human sera (1/rate of seroreversion for anti-TgERP IgG) by fitting to previously-published human serologcal data collected in Brazil.

<img align="centre" alt="Fig1" width="800px" src="https://user-images.githubusercontent.com/65221055/148261797-ceac35b8-1dfa-41e5-8d9d-be5334bdc5e9.png" />

### Associated publication
**Milne, G.**, Webster, J.P and Walker, M., 2020. Toward Improving Interventions Against Toxoplasmosis by Identifying Routes of Transmission Using Sporozoite-specific Serological Tools. *Clinical Infectious Diseases*, 71(10): [https://doi.org/10.1093/cid/ciaa428](https://academic.oup.com/cid/article/71/10/e686/5819395).

### Description of R scripts

#### model.R
The serocatalytic model, based on ordinary differential equations, which is solved to derive age-specific seroprevalences for both antibody types.

#### pars.R
Serocatalytic model parameters, including those that remain fixed, and those that are estimated during the fitting procedure. This script includes a toggle parameter ```P$data_choice``` which allows the user to change the dataset being fit when sourcing in the data (found in data.R). 

#### funcs.R
Provides all the user-defined functions (excluding the serocatalytic model, found in model.R), to help with various aspects of the workflow, including model fitting, data creation and cleaning.  

#### data.R
Loads in the data, calculates midpoints of age groups (using a function defined in funcs.R), and calculates age-specific seroprevalence with 95% confidence intervals. The loaded dataset can be changed by altering the ```P$data_choice``` toggle parameter in pars.R.

#### sensitivity.R
Performs two distinct sensitivity analyses to estimate the effect of changing several parameters on estimates of the duration of anti-TgERP IgG, and parameters defining the age-specific force of infection. Sensitivity analysis 1: varying the assumed TgERP immunoassay performance (sensitivity and specificity). Sensitivity analysis 2: varying the proportion of infections that are oocyst-acquired.

#### fitmod.R
Fits the serocatalytic model to the seroprevalence data using a maximum likelihood approach, with parameters fixed to appropriate values determined by sensitivity analyses (in sensitivity.R). Also calculates uncertainty in seroprevalence estimates using a bootstrapping approach. The script output is a dataframe containing the age-specific modelled observed seroprevalence (with 95% CIs), which can be plotted in plots.R. 

#### plots.R
Makes ggplots showing the fit of the modelled observed seroprevalence to the data, for both datasets and antibodies. Also provides ggplots of anti-TgERP IgG duration estimates for different assumed sensitivity and specificity values of the TgERP diagnostic (using data created with analysis in sensitivity.R).

### A note
There is a small (negligible) mismatch between the fitted parameter estimates for the Vieira et al. (2015) dataset when comparing the results from the paper to those obtained from the model  presented here. This is because the highest age band for these data, which has a midpoint of 54.5 years, was accidentally recorded as 55.5 years in the analyses performed prior to publication. This mistake has now been corrected in the code presented here, and makes no difference to the results nor conclusions.
