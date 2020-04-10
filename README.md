# CorrPower
Matlab & Python functions to calculate power for a correlation analysis
The calculation are based on parametric statistics to compute power given an expected correlation and a criterion for 
significance.

Python
------
See example notebook `Python/plot_pearson_analytical_power.ipynb`
Requires: numpy, scipy, matplotlib

Matlab
------
To plot curves of power vs. sample size given an expected correlation rho, use 
>>> power = plot_pearson_analytical_power(rho); 

Enjoy!

![example_output](https://github.com/jcrdubois/CorrPower/blob/master/Python/pearson_analytical_power.png)
