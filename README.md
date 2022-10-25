## Introduction
Decision-making is a fundamental part of life. Which shirt to buy? What food to eat? Where to move or go to school? Some decisions are harder than others, based on varied environmental stimuli (e.g. others' opinions or your own personal research), whereas some may be so trivial that you don't even consider them to be decisions! (One might say you're "biased" against all other possibilities.)

[The Physics of Optimal Decision Making: A Formal Analysis of Models of Performance in Two-Alternative Forced-Choice Tasks](https://psycnet.apa.org/doiLanding?doi=10.1037%2F0033-295X.113.4.700) (Bogacz et al 2006) is a seminal paper in decision theory with regards to two-alternative forced-choice (TAFC) tasks, in which an agent accumulates stochastic evidence towards one of two choices, namely a "positive" choice and a "negative" choice. Each choice has an associated threshold value -- if the agent's evidence crosses that value, the agent decides on that particular option.

[Decision Dynamics in Groups with Interacting Members](https://arxiv.org/abs/1906.04377) (Caginalp & Doiron 2017) introduces coupled agents to TAFCs, and they investigate how often agent decisions coincide. Crucially, each agent accumulates evidence independently until the one agent makes a decision, at which time the other agent receives an instantaneous "kick" in the direction of the first agent's decision.

My research takes this a step further and lets each agent receive a reward for a "correct" decision. One question we can then ask is, how can agents maximize their reward rate? And if they are coupled together, how is their reward rate influenced by their cognitive biases (represented by the respective values of their decision thresholds) and the extent to which they pay attention to others' choices (represented by the size of the "kick" they get after another agent's decision).

## Repository Guide
The above-mentioned papers by Bogacz and Caginalp contain the required prerequisite knowledge for navigating my writeups and code. I include three other PDFs in this repo:

 - ***optimal_strats_asymm_thresh.pdf*** is the primary reference for the derivations and calculations involved in my code. 
 - ***backward_forward_FPE_methods.pdf*** is a short writeup on two equivalent methods to compute a single agent's survival probability, based on the forward and backward Fokker-Planck equations. 
 - ***math_bio_11-15-21.pdf*** contains the slides I presented to the CU Boulder Mathematical Biology Seminar in November 2021.

And below is a summary of the code in this repo:

- ***parameters.m*** controls all parameters relevant to each agent: thresholds, kick magnitudes, rewards at each threshold, etc.
- Files that begin with ***"diffusionTrial"*** run Monte Carlo simulations for agent decisions. Files that begin with ***"passageTimes"*** use these MC trials to accumulate simulated first passage times.
- ***c.m***, ***dcdx.m***, and ***intc_x.m*** all contain explicit series expansions for survival probabilities and related quantities.
- Files that begin with ***RR*** calculate expected reward rates in different ways.
    - NOTE: Any file whose name has the word ***"simple"*** refers to a "simpler" type of task in which the second agent to decide is NOT allowed to continue diffusing after the first agent's decision, but instead must IMMEDIATELY decide at the threshold to which it is closest.
- ***caginalpFigReplication.m*** replicates the first few relevant figures in the above paper by Caginalp.
- ***twoAgentRR.m*** is the go-to file for figures related to two-agent reward rates.
