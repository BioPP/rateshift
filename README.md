# Bio++ RateShift

A reimplementation of the RateShift test by Pupko and Galtier:
> A covarion-based method for detecting molecular adaptation: application to the evolution of primate mitochondrial genomes.
> Pupko T, Galtier N.
> Proc Biol Sci. 2002 Jul 7;269(1498):1313-6.

This implementation compares two sets of branches, arbitrarily refered to as 'foreground' and 'background'. The program compares for each site a model where the two sets evolve with a distinct evolutionary rate (estimated by the program) to a model where the two sets evolve under the same rate. The output file contains, for each site:
* the common rate
* the foreground and background rates
* the AIC of the one-rate and two-rate models
* the p-value of the likelihood ratio test between the two models

BppRateShift follows the BppO syntax, and all sequence and tree format supported by Bio++ (see the BppSuite manual for a complete reference). Input branches are specified as a list of node ids. Such nodes can be output by programms such as bppML, or visualized using the BppPhyView.
