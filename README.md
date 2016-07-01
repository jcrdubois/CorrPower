# CorrPower
Matlab functions to calculate power for a correlation analysis

The function analCorrPower.m uses parametric statistics to compute power given an expected correlation and a criterion for significance.

To plot curves of power vs. sample size given an expected correlation rho, use power=plotAnalCorrPower(rho); this will compute power for criteria alpha = 0.05, 0.01, and 0.001, both one-sided and two-sided. It calls analCorrPower.m

See each function's help and code for further information.
