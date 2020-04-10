# CorrPower
Matlab & Python functions to calculate power for a correlation analysis

Python
------
See example notebook `Python/plot_pearson_analytical_power.ipynb`

Matlab
------
The function `pearson_analytical_power.m` uses parametric statistics to compute power given an expected correlation and a criterion for significance.

To plot curves of power vs. sample size given an expected correlation rho, use 
>>> power = plot_pearson_analytical_power(rho); 
This will compute power for criteria alpha = 0.05, 0.01, and 0.001, both one-sided and two-sided. It calls `pearson_analytical_power.m`

