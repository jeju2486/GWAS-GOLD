#!/usr/bin/env python3
def get_options():
    import argparse

    description = 'Draw a QQ-plot from pyseer lrt-pvalue results'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('table',
                        help='Pyseer output (must contain an "lrt-pvalue" column)')
    parser.add_argument('plot_file',
                        nargs='?',
                        default='qq_plot.png',
                        help='Filename for the QQ-plot (default: qq_plot.png)')

    return parser.parse_args()


def main():
    options = get_options()

    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    import matplotlib.pyplot as plt

    # Read in p-values
    m = pd.read_csv(options.table, usecols=['lrt-pvalue'], sep='\t')['lrt-pvalue']
    # Replace zero p-values with a very small number
    m[m == 0] = 1e-300

    # Convert to -log10
    y = -np.log10(m)
    x = -np.log10(np.random.uniform(0, 1, m.shape[0]))

    # Create figure
    plt.figure(figsize=(4, 3.75))
    ax = plt.subplot(111)

    fig = sm.qqplot_2samples(x, y,
                             xlabel='Expected $-\\log_{10}(p)$',
                             ylabel='Observed $-\\log_{10}(p)$',
                             line='45',
                             ax=ax)
    ax = fig.axes[0]
    points = ax.lines[0]
    ref_line = ax.lines[1]

    # Style the points: black edges, white fill, slightly transparent
    points.set_marker('o')
    points.set_markerfacecolor('white')
    points.set_markeredgecolor('black')
    points.set_linestyle('None')
    points.set_alpha(0.4)

    # Style the reference line: red, full opacity
    ref_line.set_color('red')
    ref_line.set_alpha(1.0)

    # Adjust x and y limits (as in original code)
    ax.set_xlim(-0.5, x.max() + 0.5)
    ax.set_ylim(-0.5, y.max() + 0.5)

    plt.tight_layout()
    plt.savefig(options.plot_file, dpi=150)


if __name__ == "__main__":
    main()
