# mlMS
Developmental R-Package to fit various multilevel Bayesian Markov switching models and some nested models. 
All Markov switching models are fitted through a combination of the well known forward-algorithm and various samplers from the Nimble package. 
The most complex model that can be fitted is the multilevel Bayesian Markov switching vector autoregressive (ML-MS-VAR) model. 
Other models that can be fitted, can be seen as special cases of the ML-MS-VAR model. 
The implemented ML-MS-VAR model features random effects in all model parts: the intercepts, the autoregressive effects, 
the residual variances (through nested wishart distributions) and the transition probabilities. See "test_mlMSVAR.R" for an example of how to simulate data
and fit the model. 


Also check out some of the main inspirations for this project:

Crayen, C., Eid, M., Lischetzke, T., Courvoisier, D. S., & Vermunt, J. K. (2012). 
Exploring Dynamics in Mood Regulation—Mixture Latent Markov Modeling of Ambulatory Assessment Data. 
Psychosomatic Medicine, 74(4), 366–376. https://doi.org/10.1097/PSY.0b013e31825474cb

https://github.com/SachaEpskamp/mlVAR

Altman, R. M. (2007). Mixed Hidden Markov Models: An Extension of the Hidden Markov Model to the Longitudinal Data Setting. 
Journal of the American Statistical Association, 102(477), 201–210.

de Haan-Rietdijk, S., Kuppens, P., Bergeman, C. S., Sheeber, L. B., Allen, N. B., & Hamaker, E. L. (2017). 
On the Use of Mixed Markov Models for Intensive Longitudinal Data. 
Multivariate Behavioral Research, 52(6), 747–767. https://doi.org/10.1080/00273171.2017.1370364

Turek, D., de Valpine, P., & Paciorek, C. J. (2016). 
Efficient Markov chain Monte Carlo sampling for hierarchical hidden Markov models. 
Environmental and Ecological Statistics, 23(4), 549–564. https://doi.org/10.1007/s10651-016-0353-z

https://github.com/emmekeaarts/mHMMbayes

Bringmann, L. F., Ferrer, E., Hamaker, E. L., Borsboom, D., & Tuerlinckx,
F. (n.d.). Modeling nonstationary emotion dynamics in dyads using a time-
varying vector-autoregressive model. , 53 (3), 293–314. Retrieved 2022-03-30, from
https://doi.org/10.1080/00273171.2018.1439722 (Publisher: Routledge eprint:
https://doi.org/10.1080/00273171.2018.1439722)

https://github.com/qureshlatif/QSLpersonal/blob/master/R/RunNimbleParallel.R

Asparouhov, T., Hamaker, E. L., & Muthén, B. (2018). Dynamic Structural Equation Models. Structural Equation Modeling: A Multidisciplinary Journal, 25(3), 359–388. https://doi.org/10.1080/10705511.2017.1406803


