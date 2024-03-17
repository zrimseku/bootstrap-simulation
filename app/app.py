import itertools
import time
from pathlib import Path
import pandas as pd
import numpy as np
import scipy
import seaborn as sns
from matplotlib import pyplot as plt

from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

sns.set_theme()

# Load RDS file
rds_file_path = Path(__file__).parent / "data_rds/"
r_data1 = robjects.r['readRDS'](rds_file_path / "aggregated_one_sided.rds")
r_data2 = robjects.r['readRDS'](rds_file_path / "aggregated_two_sided.rds")

# Convert R data to Python object
pandas2ri.activate()
py_data1 = pandas2ri.ri2py(r_data1)
py_data2 = pandas2ri.ri2py(r_data2)
pandas2ri.deactivate()

# Convert Python object to DataFrame
df1 = pd.DataFrame(py_data1)
df2 = pd.DataFrame(py_data2)

df = pd.read_csv(Path(__file__).parent / "coverage.csv", na_values="NA")
df = df[df['B'] == 1000]
selectable_cols = ['Confidence level', 'Statistic', 'Distribution']
sel_cols_x = ['Statistic', 'Distribution', "No X grid"]
sel_cols_y = selectable_cols + ["No Y grid"]            # Confidence level can't be on x axis because of scales
methods = df["method"].unique().tolist()
alphas1 = df["alpha"].unique().tolist()
alphas1.sort()
alphas2 = [round(b - a, 5) for a, b in zip(alphas1[:len(alphas1)//2], alphas1[len(alphas1)//2:][::-1])]
statistics = df["statistic"].unique().tolist()
distributions = df["dgp"].unique().tolist()


app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_action_button("render_btn", "Render plot with parameters selected below"),
        # ui.input_slider("width_slider", "Width slider", min=1, max=3000, value=1000),
        ui.input_radio_buttons("sided", "Confidence intervals", {"1": "One-sided", "2": "Two-sided"}),
        ui.input_selectize(
            "xgrid", "X grid", sel_cols_x, selected="No X grid"
        ),
        ui.panel_conditional(
            "input.xgrid === 'Confidence level'", ui.input_checkbox_group(
                "alphas_x", "Confidence levels", alphas1, selected=alphas1
            )),
        ui.panel_conditional(
            "input.xgrid === 'Statistic'", ui.input_checkbox_group(
                "statistics_x", "Statistics", statistics, selected=statistics
            )),
        ui.panel_conditional(
            "input.xgrid === 'Distribution'", ui.input_checkbox_group(
                "distributions_x", "Distributions", distributions, selected=distributions
            )),
        ui.input_selectize(
            "ygrid", "Y grid", sel_cols_y, selected="No Y grid"
        ),
        ui.panel_conditional(
            "input.ygrid === 'Confidence level'", ui.input_checkbox_group(
                "alphas_y", "Confidence levels", alphas1, selected=alphas1
            )),
        ui.panel_conditional(
            "input.ygrid === 'Statistic'", ui.input_checkbox_group(
                "statistics_y", "Statistics", statistics, selected=statistics
            )),
        ui.panel_conditional(
            "input.ygrid === 'Distribution'", ui.input_checkbox_group(
                "distributions_y", "Distributions", distributions, selected=distributions
            )),
        ui.panel_conditional(
            "input.xgrid != 'Confidence level' & input.ygrid != 'Confidence level'", ui.input_selectize(
            "alpha", "Confidence level", alphas1, selected=0.95
            )),
        ui.panel_conditional(
            "input.xgrid != 'Statistic' & input.ygrid != 'Statistic'", ui.input_selectize(
                "statistic", "Statistic", statistics, selected=statistics[0]
            )),
        ui.panel_conditional(
            "input.xgrid != 'Distribution' & input.ygrid != 'Distribution'", ui.input_selectize(
            "distribution", "Distribution", distributions, selected=distributions[0]
            )),
        ui.input_checkbox_group(
            "methods", "Filter by methods", methods, selected=methods
        ),
        ui.hr(),
        # ui.input_switch("show_margins", "Show marginal plots", value=True),
    ),
    ui.output_ui("plot")
    # ui.card(
    # ui.page_fluid(
    #     ui.output_plot("coverages", height="100%"),
    # ),
)


def compare_cov_dis_grid(df=None, comparing='coverage', filter_by={'alpha': [0.95]}, x='n', row='statistic', col='dgp',
                         hue='method', save_add=None, title=None, ci=95, scale='linear', folder_add='', subfolder='',
                         set_ylim=False, colors=None):
    if df is None:
        df = pd.read_csv(f'results{folder_add}/{comparing}.csv')

    for key in filter_by.keys():
        df = df[df[key].isin(filter_by[key])]

    if colors is None:
        nm = df['method'].nunique()
        if nm > 10:
            cols = plt.cm.tab20(np.linspace(0.05, 0.95, df['method'].nunique()))
        else:
            cols = plt.cm.tab10(np.linspace(0.05, 0.95, df['method'].nunique()))
        colors = {m: c for (m, c) in zip(df['method'].unique(), cols)}

    g = sns.FacetGrid(df, row=row, col=col, margin_titles=True, sharex=True, sharey='row', palette=colors, aspect=2,
                      height=3, legend_out=True)
    g.figure.set_size_inches(10, 20)
    if comparing == 'coverage':
        g.map_dataframe(plot_coverage_bars, colors=colors, ci=ci, scale=scale, set_ylim=set_ylim,
                        order=df[hue].unique(), hue=hue, x=x)
    else:
        g.map(sns.boxplot, x, comparing, hue, hue_order=df[hue].unique(), fliersize=0,
              whis=[(100 - ci) / 2, 50 + ci / 2],
              palette=colors)
        ylim = np.nanquantile(df['distance'], (0.01, 0.99))
        g.set(ylim=ylim)

        for axs in g.axes:
            for ax in axs:
                ax.axhline(0, linestyle='--', color='gray')

    # g.add_legend(title='method')
    # sns.move_legend(g, 'center left', bbox_to_anchor=(.9, 0.5), borderaxespad=1)

    plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
    # plt.tight_layout()


    if title is not None:
        g.fig.subplots_adjust(top=0.95)
        g.fig.suptitle(title, fontsize=16)


def plot_coverage_bars(data, **kwargs):
    colors = kwargs['colors']
    ci = kwargs['ci']
    scale = kwargs['scale']

    if 'cov_kind' in kwargs:                # for the possibility of plotting variance coverage with it
        cov_kind = kwargs['cov_kind']
    else:
        cov_kind = 'coverage'

    data['ci'] = np.sqrt(data[cov_kind] * (1 - data[cov_kind]) / data['repetitions'])
    if ci != 'se':
        data['ci'] *= scipy.stats.norm.ppf(0.5 + ci / 200)
    data['low'] = data[cov_kind] - data['ci']

    n_levels = len(kwargs['order'])
    group_width = 0.8
    bar_width = group_width / n_levels
    offsets = np.linspace(0, group_width - bar_width, n_levels)
    offsets -= offsets.mean()

    bar_pos = np.arange(data[kwargs['x']].nunique())
    for i, method in enumerate(kwargs['order']):
        data_m = data[data[kwargs['hue']] == method]
        offset = bar_pos + offsets[i]
        if data_m['ci'].shape[0] == 0:
            continue
        plt.bar(offset, data_m['ci'], bar_width, bottom=data_m[cov_kind], label=method, color=colors[method],
                ec=colors[method])
        plt.bar(offset, data_m['ci'], bar_width, bottom=data_m['low'], color=colors[method], ec=colors[method])

    for p in bar_pos[:-1]:
        plt.axvline(p + 0.5, ls=':', alpha=0.2)

    a = data['alpha'].values[0]
    if a > 0.9:
        if scale == 'logit':
            ylim = (0.8, 0.99)
        else:
            ylim = (0.8, 1)
    elif a < 0.1:
        ylim = (0, 0.2)
    else:
        ylim = (a - 0.1, a + 0.1)

    ax = plt.gca()

    ax.set_yscale(scale)

    if kwargs['set_ylim']:
        ax.set(ylim=ylim)

    ax.axhline(a, linestyle='--', color='gray')
    plt.yticks(list(plt.yticks()[0]) + [a])

    ax.set_xlabel(kwargs['x'])
    ax.set_ylabel(cov_kind)
    plt.xticks(bar_pos, sorted(data[kwargs['x']].unique()))


def server(input: Inputs):

    @reactive.effect()
    def update_choices():
        """Updates possible choices based on selected values."""
        # setting correct alphas for one or two-sided intervals
        if input.sided() == '1':
            ui.update_selectize('alpha', choices=alphas1)
            selected_alphas_x = input.alphas_x()
            ui.update_checkbox_group('alphas_x', choices=alphas1,
                                     selected=[str(a) for a in alphas1 if str(a) in selected_alphas_x])
            ui.update_checkbox_group('alphas_y', choices=alphas1,
                                     selected=[str(a) for a in alphas1 if str(a) in input.alphas_y()])
        else:
            ui.update_selectize('alpha', choices=alphas2)
            ui.update_checkbox_group('alphas_x', choices=alphas2,
                                     selected=[str(a) for a in alphas2 if str(a) in input.alphas_x()])
            ui.update_checkbox_group('alphas_y', choices=alphas2,
                                     selected=[str(a) for a in alphas2 if str(a) in input.alphas_y()])

        bootstrap_methods = ['percentile', 'basic', 'bca', 'bc', 'standard', 'smoothed', 'double', 'studentized']
        other_methods = {'mean': ['wilcoxon', 'ttest'], 'std': ['chi_sq'],
                         'median': ['wilcoxon', 'ci_quant_param', 'ci_quant_nonparam', 'maritz-jarrett'],
                         'percentile_5': ['ci_quant_param', 'ci_quant_nonparam', 'maritz-jarrett'],
                         'percentile_95': ['ci_quant_param', 'ci_quant_nonparam', 'maritz-jarrett'],
                         'corr': ['ci_corr_pearson']}

        # setting correct methods for each statistic
        if 'Statistic' not in [input.xgrid(), input.ygrid()]:

            possible_methods = bootstrap_methods + other_methods[input.statistic()]

            # setting correct distributions for selected statistic
            possible_dist = [d for d in distributions if d[:9] != 'DGPBiNorm'] if input.statistic() != 'corr' else \
                [d for d in distributions if d[:9] == 'DGPBiNorm']
            selected_dist_possible = input.distribution() in possible_dist
            ui.update_selectize('distribution', choices=possible_dist,
                                selected=input.distribution() if selected_dist_possible else possible_dist[0])

        else:
            stats = input.statistics_x() if 'Statistic' == input.xgrid() else input.statistics_y()
            possible_methods = bootstrap_methods + list(itertools.chain(*[other_methods[s] for s in stats]))

        selected_methods = input.methods()
        ui.update_checkbox_group('methods', choices=possible_methods,
                                 selected=[m for m in possible_methods if m in selected_methods])

        # not possible to select the same dimension on the X and Y axis of the grid
        if input.xgrid() != 'No X grid':
            ui.update_selectize('ygrid', choices=[c for c in sel_cols_y if c != input.xgrid()],
                                selected=input.ygrid())
        else:
            ui.update_selectize('ygrid', choices=sel_cols_y, selected=input.ygrid())

        if input.ygrid() != 'No Y grid':
            ui.update_selectize('xgrid', choices=[c for c in sel_cols_x if c != input.ygrid()],
                                selected=input.xgrid())
        else:
            ui.update_selectize('xgrid', choices=sel_cols_x, selected=input.xgrid())


    @reactive.Calc()
    def filtered_df() -> pd.DataFrame:
        """Returns a Pandas data frame that includes only the desired rows"""
        # This calculation "req"uires that at least one species is selected
        req(len(input.methods()) > 0)

        fil_df = df.copy()
        if 'Distribution' not in [input.xgrid(), input.ygrid()]:
            fil_df = fil_df[(fil_df['dgp'] == input.distribution())]
        if 'Statistic' not in [input.xgrid(), input.ygrid()]:
            fil_df = fil_df[(fil_df['statistic'] == input.statistic())]

        if input.sided() == '1':
            if 'Confidence level' not in [input.xgrid(), input.ygrid()]:
                fil_df = fil_df[(fil_df['alpha'] == float(input.alpha()))]

        else:
            al = round((1 - float(input.alpha())) / 2, 5)
            au = 1 - al
            if 'Confidence level' not in [input.xgrid(), input.ygrid()]:
                fil_df = fil_df[fil_df['alpha'] == au]

            fil_df2 = df.copy()
            if 'Distribution' not in [input.xgrid(), input.ygrid()]:
                fil_df2 = fil_df2[(fil_df2['dgp'] == input.distribution())]
            if 'Statistic' not in [input.xgrid(), input.ygrid()]:
                fil_df2 = fil_df2[(fil_df2['statistic'] == input.statistic())]
            if 'Confidence level' not in [input.xgrid(), input.ygrid()]:
                fil_df2 = fil_df2[fil_df2['alpha'] == al]

            # merge dataframes to be able to subtract coverages of same experiments
            fil_df = pd.merge(fil_df, fil_df2, on=['method', 'dgp', 'statistic', 'n', 'B', 'repetitions'],
                              suffixes=('_au', '_al'))
            fil_df['coverage'] = fil_df['coverage_au'] - fil_df['coverage_al']
            fil_df['alpha'] = fil_df['alpha_au'] - fil_df['alpha_al']

        return fil_df

    @render.ui
    @reactive.event(input.render_btn)
    def plot():
        if input.xgrid() == 'No X grid':
            w = 1000
        else:
            if input.xgrid() == 'Statistic':
                w = len(input.statistics_x()) * 500
            elif input.xgrid() == 'Confidence level':
                w = len(input.alphas_x()) * 500
            else:
                w = len(input.distributions_x()) * 500

        if input.ygrid() == 'No Y grid':
            h = 500
        else:
            if input.ygrid() == 'Statistic':
                h = len(input.statistics_y()) * 250
            elif input.ygrid() == 'Confidence level':
                h = len(input.alphas_y()) * 250
            else:
                h = len(input.distributions_y()) * 250

        return ui.div(
            render.plot(_fn=coverages, width=w, height=h)
        )

    def coverages():
        with reactive.isolate():
            """Generates a plot for Shiny to display to the user"""
            # req(update_choices()) doesn't help (+ needs Effect to Calc)
            # ne dela -> vse kar caka in se updejta gre na stack in se pol izrise. Pogruntat rabim kako lahko izrisem samo
            # tazadn plot, ga klicem sele ko se do konca poupdajtajo zadeve

            # The plotting function to use depends on whether we are plotting grids
            no_grids = input.xgrid() == "No X grid" and input.ygrid() == "No Y grid"

            current_df = filtered_df()
            current_methods = [m for m in input.methods() if m in current_df['method'].unique()]

            if no_grids:

                nm = len(current_methods)
                if nm > 10:
                    cols = plt.cm.tab20(np.linspace(0.05, 0.95, nm))
                else:
                    cols = plt.cm.tab10(np.linspace(0.05, 0.95, nm))
                colors = {m: c for (m, c) in zip(current_methods, cols)}

                plt.figure(figsize=(10, 12))

                plot_coverage_bars(data=current_df, colors=colors, ci=95, scale='linear', set_ylim=True,
                                   order=current_methods, hue='method', x='n')
                # plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda value, _: f'{value:.2f}'))

                handles, labels = plt.gca().get_legend_handles_labels()
                plt.legend(handles, labels, loc='center left', title="Method", bbox_to_anchor=(1, 0.5))

            else:
                row_col_dict = dict(zip(selectable_cols + ["No X grid", "No Y grid"],
                                        ['alpha', 'statistic', 'dgp', None, None]))

                filters = {'method': input.methods()}
                if 'Confidence level' == input.xgrid():
                    filters['alpha'] = [float(a) for a in input.alphas_x()]
                elif 'Confidence level' == input.ygrid():
                    filters['alpha'] = [float(a) for a in input.alphas_y()]
                else:
                    filters['alpha'] = [float(input.alpha())]
                if 'Statistic' in [input.xgrid(), input.ygrid()]:
                    filters['statistic'] = input.statistics_x() if 'Statistic' == input.xgrid() else input.statistics_y()
                else:
                    filters['statistic'] = [input.statistic()]
                if 'Distribution' in [input.xgrid(), input.ygrid()]:
                    filters['dgp'] = input.distributions_x() if 'Distribution' == input.xgrid() else input.distributions_y()
                else:
                    filters['dgp'] = [input.distribution()]

                compare_cov_dis_grid(df=current_df, comparing='coverage', filter_by=filters, x='n',
                                     row=row_col_dict[input.ygrid()], col=row_col_dict[input.xgrid()],
                                     hue='method', save_add=None, title=None, ci=95, scale='linear',
                                     set_ylim=True, colors=None)

                # plt.legend([], [], frameon=False)
                # handles, labels = plt.gca().get_legend_handles_labels()
                # plt.legend(handles, labels, loc='center left', title="Method", bbox_to_anchor=(1, 0.5))


app = App(app_ui, server)
# app.run()
