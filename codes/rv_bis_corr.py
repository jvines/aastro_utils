"""rv_bis_corr.

Author: Jose Vines

Calculate and plot RV vs BIS correlation
"""

import scipy as sp
import statsmodels.api as sm
from scipy.stats import pearsonr
import scipy.stats as st
import matplotlib.pyplot as plt


def rv_bis_corr(data, confidence=0.05, name='last'):
    """Calculate RV vs BIS correlation and plot it.

    Parameters
    ----------
    data : dict
        A dictionary containing the datasets. Each dataset must be an array
        of size (5, m) in the following order: t, x, y, xerr, yerr.

    confidence : float
        The confidence level.

    name : str, optional
        Target name for saving the plot.

    """
    # Linear Model fitting
    tlow = sp.inf
    tup = -sp.inf
    x = sp.array([])
    y = sp.array([])
    for key in data.keys():
        tl = data[key][:, 0].min()
        tu = data[key][:, 0].max()
        if tl < tlow:
            tlow = tl
        if tu > tup:
            tup = tu
        x = sp.concatenate((x, data[key][:, 1]))
        y = sp.concatenate((y, data[key][:, 3]))

    r, p_val = pearsonr(x, y)

    X = sm.add_constant(x)
    model = sm.OLS(y, X)
    fitted = model.fit()

    error_kwargs = {'lw': .75, 'zorder': 0}

    # Confidence interval calculation
    y_hat = fitted.predict(X)
    y_err = y - y_hat
    x_mean = X.T[1].mean()
    n = len(x)
    dof = n - fitted.df_model - 1  # Degrees of freedom

    # 2 tailed t-stat calculation
    t = st.t.ppf(1 - confidence / 2, df=dof)
    s_err = sp.sum(sp.power(y_err, 2))

    markers = ['o', 'v', '^', '>', '<', '8', 's', 'p', 'H', 'D', '*', 'd']

    f, ax = plt.subplots(figsize=(20, 10))
    ims = []
    for i, key in enumerate(data.keys()):
        x = data[key][:, 1]
        y = data[key][:, 3]
        xerr = data[key][:, 2]
        yerr = data[key][:, 4]
        ti = data[key][:, 0]
        im = ax.scatter(
            x, y, marker=markers[i], edgecolors='k', c=ti, cmap='cool_r', s=180
        )
        ims.append(im)
        ax.errorbar(
            x, y, xerr=xerr, yerr=yerr, marker=None,
            linestyle='', ecolor='k', **error_kwargs
        )
    for im in ims:
        im.set_clim(tlow, tup)

    xmin, xmax = ax.get_xlim()

    x_pred = sp.linspace(xmin, xmax, 1000)
    x_pred2 = sm.add_constant(x_pred)
    y_pred = fitted.predict(x_pred2)

    conf = t * sp.sqrt((s_err / (n - 2)) *
                       (1. / n + (sp.power((x_pred - x_mean), 2) /
                                  ((sp.sum(sp.power(x_pred, 2))) - n *
                                   (sp.power(x_mean, 2))))))

    upper = y_pred + abs(conf)
    lower = y_pred - abs(conf)
    cb = f.colorbar(ims[-1], pad=.005)
    lab = 'Pearson r: {:.3f}'.format(r)
    ax.plot(x_pred, y_pred, '-', color='midnightblue', linewidth=2, label=lab)
    ax.fill_between(x_pred, lower, upper, color='#888888', alpha=0.4)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel('RV (km s$^{-1}$)', fontsize=30)
    ax.set_ylabel('Bisector Velocity Span (km s$^{-1}$)', fontsize=30)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(28)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(28)
    cb.ax.tick_params(labelsize=28)
    cb.set_label('JD - 2450000', rotation=270, labelpad=25, fontsize=30)
    fname = '{}_bisector_rv.pdf'.format(name)
    plt.legend(loc=0, prop={'size': 28})
    plt.savefig(fname, bbox_inches='tight')
    return fitted
