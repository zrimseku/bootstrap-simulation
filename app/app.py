import itertools
from pathlib import Path
import pandas as pd
import numpy as np
import scipy
import seaborn as sns
from matplotlib import pyplot as plt

from shiny import App, Inputs, reactive, render, req, ui

sns.set_theme()

df1 = pd.read_csv(Path(__file__).parent / "aggregated_one_sided.csv")
df2 = pd.read_csv(Path(__file__).parent / "aggregated_two_sided.csv")

df1 = df1[df1['B'] == 1000]
df2 = df2[df2['B'] == 1000]
repetitions = 10000
df1['repetitions'] = repetitions
df2['repetitions'] = repetitions
selectable_cols = ['Confidence level', 'Statistic', 'Distribution']
sel_cols_x = ['Statistic', 'Distribution', "No X grid"]
sel_cols_y = selectable_cols + ["No Y grid"]            # Confidence level can't be on x axis because of scales
methods = df1["method"].unique().tolist()
alphas1 = df1["nominal_coverage"].unique().tolist()
alphas2 = df2["nominal_coverage"].unique().tolist()
statistics = df1["functional"].unique().tolist()
distributions = df1["dgp"].unique().tolist()

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_action_button("render_btn", "Render plot with parameters selected below"),
        ui.panel_conditional(
            "input.sided === '1'", ui.input_radio_buttons("criteria", "Criteria to draw",
                                                          {"coverage": "Coverage",
                                                           "abs_dist_to_exact": "Distance from exact"})),
        ui.panel_conditional(
            "input.sided === '2'", ui.input_radio_buttons("criteria2", "Criteria to draw",
                                                          {"coverage": "Coverage", "length": "Length"})),
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
            "alpha", "Confidence level", alphas1, selected='0.05'
            )),
        ui.panel_conditional(
            "input.xgrid != 'Statistic' & input.ygrid != 'Statistic'", ui.input_selectize(
                "statistic", "Statistic", statistics, selected='mean'
            )),
        ui.panel_conditional(
            "input.xgrid != 'Distribution' & input.ygrid != 'Distribution'", ui.input_selectize(
            "distribution", "Distribution", distributions, selected='DGPNorm_0_1'
            )),
        ui.input_checkbox_group(
            "methods", "Filter by methods", methods, selected=methods
        ),
        ui.hr(),
    ),
    ui.output_ui("plot")

)


def compare_cov_dis_grid(df, comparing='coverage', filter_by={'nominal_coverage': [0.95]}, x='n', row='functional',
                         col='dgp', hue='method', title=None, ci=95, scale='linear', set_ylim=False, colors=None):

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

    g.map_dataframe(plot_coverage_bars, colors=colors, ci=ci, scale=scale, set_ylim=set_ylim,
                    order=df[hue].unique(), hue=hue, x=x, comparing=comparing)


    plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))

    if title is not None:
        g.fig.subplots_adjust(top=0.95)
        g.fig.suptitle(title, fontsize=16)


def plot_coverage_bars(data, **kwargs):
    colors = kwargs['colors']
    ci = kwargs['ci']
    scale = kwargs['scale']
    comparing = kwargs['comparing']

    data['ci'] = data[comparing + '_sd'] / np.sqrt(data['repetitions']) * scipy.stats.norm.ppf(0.5 + ci / 200)
    data['low'] = data[comparing] - data['ci']

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
        plt.bar(offset, data_m['ci'], bar_width, bottom=data_m[comparing], label=method, color=colors[method],
                ec=colors[method])
        plt.bar(offset, data_m['ci'], bar_width, bottom=data_m['low'], color=colors[method], ec=colors[method])

    for p in bar_pos[:-1]:
        plt.axvline(p + 0.5, ls=':', alpha=0.2)

    a = data['nominal_coverage'].values[0]

    ax = plt.gca()

    ax.set_yscale(scale)

    if kwargs['set_ylim']:
        if comparing == 'abs_dist_to_exact':
            data_bt = data[data['method'] == 'B-t']
            data_not_bt = data[data['method'] != 'B-t']
            if max(data_bt[comparing] + data_bt['ci']) > (max(data_not_bt[comparing] + data_not_bt['ci']) * 50):
                ylim = [min(-0.01, min(data_not_bt['low'])), max(data_not_bt[comparing] + data_not_bt['ci']) * 1.05]
            else:
                ylim = [min(-0.01, min(data['low'])), max(data[comparing] + data['ci']) * 1.05]
            ax.set(ylim=ylim)
        elif comparing == 'coverage':
            if a > 0.9:
                if scale == 'logit':
                    ylim = (0.8, 0.99)
                else:
                    ylim = (0.8, 1)
            elif a < 0.1:
                ylim = (0, 0.2)
            else:
                ylim = (a - 0.1, a + 0.1)
            ax.set(ylim=ylim)

    if comparing == 'coverage':
        ax.axhline(a, linestyle='--', color='gray')
        plt.yticks(list(plt.yticks()[0]) + [a])

    ax.set_xlabel(kwargs['x'])
    ax.set_ylabel(comparing)
    plt.xticks(bar_pos, sorted(data[kwargs['x']].unique()))


def server(input: Inputs):

    @reactive.effect()
    def update_choices():
        """Updates possible choices based on selected values."""
        # setting correct nominal_coverages for one or two-sided intervals
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

        bootstrap_methods = ['PB', 'BB', 'BCa', 'BC', 'B-n', 'SB', 'DB', 'B-t']
        other_methods = {'mean': ['wilcoxon', 't-test'], 'std': ['chi-sq'],
                         'median': ['wilcoxon', 'q-par', 'q-nonpar', 'm-j'],
                         'Q(0.05)': ['q-par', 'q-nonpar', 'm-j'],
                         'Q(0.95)': ['q-par', 'q-nonpar', 'm-j'],
                         'corr': ['fisher']}

        # setting correct methods for each functional
        if 'Statistic' not in [input.xgrid(), input.ygrid()]:

            # Check if current distribution is Bernoulli
            current_dist = input.distribution()
            if current_dist.startswith('DGPBernoulli'):
                # For Bernoulli distributions, only show t-test, agresti_coull, clopper_pearson for mean
                if input.statistic() == 'mean':
                    possible_methods = bootstrap_methods + ['t-test', 'a-c', 'c-p']
                else:
                    possible_methods = bootstrap_methods + other_methods[input.statistic()]
            else:
                possible_methods = bootstrap_methods + other_methods[input.statistic()]

            # setting correct distributions for selected statistic
            if input.statistic() == 'corr':
                possible_dist = [d for d in distributions if d[:9] == 'DGPBiNorm']
            elif input.statistic() == 'mean':
                possible_dist = [d for d in distributions if d[:9] != 'DGPBiNorm']
            else:
                possible_dist = [d for d in distributions if d[:9] != 'DGPBiNorm' and not d.startswith('DGPBernoulli')]
            
            selected_dist_possible = input.distribution() in possible_dist
            ui.update_selectize('distribution', choices=possible_dist,
                                selected=input.distribution() if selected_dist_possible else possible_dist[0])

        else:
            stats = input.statistics_x() if 'Statistic' == input.xgrid() else input.statistics_y()
            
            # Check if any Bernoulli distributions are selected in the grid
            selected_dists = input.distributions_x() if 'Distribution' == input.xgrid() else input.distributions_y()
            has_bernoulli = any(d.startswith('DGPBernoulli') for d in selected_dists)
            
            # Create a copy of other_methods to modify
            grid_methods = other_methods.copy()
            if has_bernoulli and 'mean' in stats:
                grid_methods['mean'] = ['t-test', 'a-c', 'c-p']
                
            possible_methods = bootstrap_methods + list(itertools.chain(*[grid_methods[s] for s in stats]))

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
        fil_dfs = []

        if input.sided() == '1':
            fil_df = df1.copy()
        else:
            fil_df = df2.copy()
        if 'Distribution' not in [input.xgrid(), input.ygrid()]:
            fil_df = fil_df[(fil_df['dgp'] == input.distribution())]
        if 'Statistic' not in [input.xgrid(), input.ygrid()]:
            fil_df = fil_df[(fil_df['functional'] == input.statistic())]

        # if input.sided() == '1':
        if 'Confidence level' not in [input.xgrid(), input.ygrid()]:
            fil_df = fil_df[(fil_df['nominal_coverage'] == float(input.alpha()))]

        return fil_df

    @render.ui
    @reactive.event(input.render_btn)
    def plot():
        if input.xgrid() == 'No X grid':
            w = 1000
        else:
            if input.xgrid() == 'Statistic':
                if input.distribution()[:12] == 'DGPBernoulli':
                    w = 1000
                else:
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
            # The plotting function to use depends on whether we are plotting grids
            no_grids = input.xgrid() == "No X grid" and input.ygrid() == "No Y grid"

            current_df = filtered_df()
            current_methods = [m for m in input.methods() if m in current_df['method'].unique()]
            comparing = input.criteria() if input.sided() == '1' else input.criteria2()

            if no_grids:

                nm = len(current_methods)
                if nm > 10:
                    cols = plt.cm.tab20(np.linspace(0.05, 0.95, nm))
                else:
                    cols = plt.cm.tab10(np.linspace(0.05, 0.95, nm))
                colors = {m: c for (m, c) in zip(current_methods, cols)}

                plt.figure(figsize=(10, 12))

                plot_coverage_bars(data=current_df, colors=colors, ci=95, scale='linear', set_ylim=True,
                                   order=current_methods, hue='method', x='n', comparing=comparing)

                handles, labels = plt.gca().get_legend_handles_labels()
                plt.legend(handles, labels, loc='center left', title="Method", bbox_to_anchor=(1, 0.5))

            else:
                row_col_dict = dict(zip(selectable_cols + ["No X grid", "No Y grid"],
                                        ['nominal_coverage', 'functional', 'dgp', None, None]))

                filters = {'method': input.methods()}
                if 'Confidence level' == input.xgrid():
                    filters['nominal_coverage'] = [float(a) for a in input.alphas_x()]
                elif 'Confidence level' == input.ygrid():
                    filters['nominal_coverage'] = [float(a) for a in input.alphas_y()]
                else:
                    filters['nominal_coverage'] = [float(input.alpha())]
                if 'Statistic' in [input.xgrid(), input.ygrid()]:
                    filters['functional'] = input.statistics_x() if 'Statistic' == input.xgrid() else input.statistics_y()
                else:
                    filters['functional'] = [input.statistic()]
                if 'Distribution' in [input.xgrid(), input.ygrid()]:
                    filters['dgp'] = input.distributions_x() if 'Distribution' == input.xgrid() else input.distributions_y()
                else:
                    filters['dgp'] = [input.distribution()]

                compare_cov_dis_grid(df=current_df, comparing=comparing, filter_by=filters, x='n',
                                     row=row_col_dict[input.ygrid()], col=row_col_dict[input.xgrid()],
                                     hue='method', title=None, ci=95, scale='linear',
                                     set_ylim=True, colors=None)


app = App(app_ui, server)
# app.run()
# shinylive export app docs