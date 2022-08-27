# POPAN model
A modified version of POPAN model with transience and per capita recruitment parameters

# Description
Our model was a modification from a standard POPAN model by Schwarz and Arnason (1996). This modified POPAN model was specifically written to examine demographic parameters of reef manta rays (_Mobula alfredi_) in two marine protected areas (MPAs): Dampier Strait MPA and South East Misool MPA in Raja Ampat, eastern Indonesia.

The parameters estimated using this modified POPAN model include expected population sizes, detection probabilities, apparent survival probabilities, and per capita recruitment rates. These parameters are separately estimated using 11 year of reef manta ray sighting data collected from the two MPAs. 

We considered each year to be an occasion, and so the data required by our model is a capture history for each detected individual, indicating which years in which each individual was detected.
The codes in this repository allows us to conduct some procedures, including goodness-of-fit testing, model fitting that incorporate transience and/or per capita recruitment rate, model selection, and model averaging.

# Usage
Start the analysis by running “manta-analysis.r” and this will source other files including “main.r”, “internals.r”, “internals-t.r”, and “popan-ll.cpp”. The manta ray capture history data “misool_dampier_caphist.RData” will also be loaded.

A brief overview of functions to fit a POPAN model can be found in “main.r”. All procedures for undertaking the main analyses using the modified POPAN model can be found in “manta-analysis.r.”, including goodness-of-fit (GOF) tests, model fitting and model averaging, plotting the results, and checking model summary.
