import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

def compute_power(
    r:np.ndarray,
    alpha:float, 
    Ns:np.ndarray,
    tail:str='right'  
) -> np.ndarray:
    r"""Compute analytical statistical power for a Pearson correlation, 
    given an effect size `r`, a significance level `$\alpha$`, and the sample size `N`. 
    
    Parameters
    ----------
    r : np.ndarray
        expected Pearson correlation (effect size)
    alpha : float
        threshold p-value for rejecting the null hypothesis (Type I error rate)
    Ns : np.ndarray
        sample size(s)
    tail : str, optional
        tail for the statistical test, by default 'right'
    
    Returns
    -------
    np.ndarray
        the power (`$\beta$`) of the experiment

    Example
    -------
    >>> from pearson_analytical import compute_power
    >>> power = compute_power(r=0.4, alpha=0.05, N=100, tail='right')
    >>> print(f'The power of this experiment is: {power:.2f}')
    """
    C = 0.5 * np.log( (1 + r) / ( 1 - r ) )

    if tail in ['left', 'right']:
        z_alpha = norm.ppf(1-alpha)
    elif tail == 'both':
        z_alpha = norm.ppf(1-alpha/2)
    else:
        raise RuntimeError(f'{tail} is not a valid value for tail')

    z_beta = C * np.sqrt(Ns - 2) - z_alpha

    return norm.cdf(z_beta)


def plot_power(
    r: np.ndarray = .4, 
    target_power: float = .95, 
    alphas: np.ndarray = [0.05, 0.01, 0.001],
    tail: str = 'right', 
    Nmax: int = 200, 
    ax = None
) -> np.ndarray:
    """make a plot of power vs. sample size for a Pearson correlation
    using the analytical power calculation. 
    Show results for different criteria and tails.
    
    Parameters
    ----------
    r : np.ndarray, optional
        expected effect size, by default .4
    target_power : float, optional
        targeted power level, by default .95
    alphas : np.ndarray, optional
        significance criteria to plot, by default [0.05, 0.01, 0.001]
    tail : str, optional
        tail of the statistical test, by default 'right'
    Nmax : int, optional
        maximum sample size to plot, by default 200
    ax : plt.Axes, optional
        axis to plot into, by default None (spawn new figure)
    
    Returns
    -------
    np.ndarray
        power, size [len(alphas), Nmax - 5]

    Example
    -------
    >>> from pearson_analytical import plot_power
    >>> power = plot_power(r=0.4)
    """    
    if isinstance(alphas, float):
        alphas = np.array([alphas])

    Ns = np.arange(5, Nmax+1).astype(int)
    # DO THE WORK

    power = np.zeros((len(alphas),len(Ns)))
    for i, alpha in enumerate(alphas):
        power[i, :] = compute_power(r=r, alpha=alpha, Ns=Ns, tail=tail)

    target_N = np.zeros(len(alphas)).astype(int)
    for i, alpha in enumerate(alphas):
        indmini = np.argmin(np.abs(power[i, :] - target_power))
        target_N[i] = Ns[indmini]

    # plot
    if ax is None:
        _, ax = plt.subplots(figsize=(10,10))

    for i, alpha in enumerate(alphas):
        ax.plot(Ns, power[i, :], label=r'$\alpha$=' + f'{alpha:.4f} -- N={target_N[i]:d}')
        ax.axvline(x=target_N[i], color=ax.get_lines()[-1].get_color(), ls=':')

    ax.legend()
    ax.set_ylabel(r'analytical power ($\beta$)')
    ax.set_xlabel('Sample size')
    ax.set_title(f'Expected effect size: r = {r:.3f}, tail = {tail}')

    # plot target power
    ax.axhline(y=target_power, color='k', ls=':')

    return power