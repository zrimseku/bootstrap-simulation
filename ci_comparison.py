import itertools
import os
import time
from itertools import repeat
from multiprocessing import Pool

import numpy
import numpy as np
import pandas as pd
import scipy
from tqdm import tqdm
from scipy.stats import norm
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt

import bootstrap_ci as boot
from generators import DGP, DGPNorm, DGPExp, DGPBeta, DGPBiNorm, DGPLogNorm, DGPLaplace, DGPRandEff

# TODO set correct R folder
os.environ['R_HOME'] = "C:/Users/ursau/anaconda3/envs/bootstrap/Lib/R"

import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()


class CompareIntervals:

    def __init__(self, statistic: callable, methods: list[str], data_generator: DGP, n: int, b: int,
                 alphas: list[float], quantile_type='median_unbiased', use_jit: bool = True,
                 sampling: str = 'nonparametric', sampling_parameters=None, group_sizes: list = None):
        self.statistic = statistic
        self.methods = methods
        self.methods_hierarchical = []
        self.dgp = data_generator
        self.n = n
        self.b = b
        self.alphas = np.array(alphas)  # we are interested in one sided intervals, two sided can be computed from them
        self.quantile_type = quantile_type
        self.computed_intervals = {m: defaultdict(list) for m in methods}       # will add all computed intervals
        self.inverse_cdf = {}
        self.times = defaultdict(list)
        self.coverages = {}
        self.distances_from_exact = {}
        self.lengths = {}
        self.use_jit = use_jit
        self.sampling = sampling

        # used for hierarchical sampling
        if sampling == 'hierarchical':
            self.sampling_parameters = sampling_parameters
            if group_sizes is None:
                max_groups_l1 = np.random.randint(2, min(5, int(self.n / 2) + 1))
                self.group_sizes = [max_groups_l1, 2 * self.n / max_groups_l1]
            else:
                self.group_sizes = group_sizes
            self.gt_variance = sum(np.array(self.dgp.stds) ** 2)     # ground truth variance
            self.variances = {m + sstr: {'mean': [], '95CI_l': [], '95CI_h': []}
                              for sstr in ['_' + ''.join([str(s) for s in strategy])
                                           for strategy in sampling_parameters['strategies']]
                              for m in sampling_parameters['sampling_methods']}     # sampling variances

    def compute_bootstrap_intervals(self, data: np.array):
        # initialize and sample so we will have the same bootstrap samples for all bootstrap methods
        t = time.time()
        bts = boot.Bootstrap(data, self.statistic, self.use_jit)
        bts.sample(self.b)
        bts.evaluate_statistic()
        ts = time.time() - t            # time needed for sampling (will add it to all methods)
        for method in self.methods:
            if method not in bts.implemented_methods:
                # skip non-bootstrap methods during calculation of bootstrap CI
                continue
            t = time.time()
            cis = bts.ci(coverages=self.alphas, side='one', method=method, quantile_type=self.quantile_type)
            for a, ci in zip(self.alphas, cis):
                self.computed_intervals[method][a].append(ci)
            self.times[method].append(time.time() - t + ts)
        return bts

    def compute_bootstrap_intervals_hierarchical(self, data: np.array, group_indices: list):
        btss = []
        for sampling_method in self.sampling_parameters['sampling_methods']:
            for strategy in self.sampling_parameters['strategies']:
                strategy_bool = [bool(i) for i in strategy]
                strategy_str = '_' + ''.join([str(s) for s in strategy]) if sampling_method == 'cases' else ''
                if any(strategy):
                    if sampling_method == 'random-effect':
                        # only run random-effect bootstrap once
                        continue
                else:
                    if sampling_method == 'cases':
                        # skip if we sample without replacement on all levels, as this would just be original set
                        continue
                t = time.time()
                bts = boot.Bootstrap(data, self.statistic, self.use_jit, group_indices=group_indices)
                bts.sample(self.b, sampling='hierarchical',
                           sampling_args={'method': sampling_method, 'strategy': strategy_bool})
                bts.evaluate_statistic(sampling='hierarchical',  sampling_args={'method': sampling_method,
                                                                                'strategy': strategy_bool})
                ts = time.time() - t            # time needed for sampling (will add it to all methods)

                # calculation of variance ratio
                var_bts = compute_hierarchical_var(bts)
                self.variances[sampling_method + strategy_str]['mean'].append(np.mean(var_bts))
                self.variances[sampling_method + strategy_str]['95CI_l'].append(np.quantile(var_bts, 0.025,
                                                                                            method='median_unbiased'))
                self.variances[sampling_method + strategy_str]['95CI_h'].append(np.quantile(var_bts, 0.975,
                                                                                            method='median_unbiased'))

                for method in self.methods:
                    if method not in bts.implemented_methods:
                        # skip non-bootstrap methods during calculation of bootstrap CI
                        continue
                    method_str = sampling_method + strategy_str + '_' + method
                    if method_str not in self.methods_hierarchical:
                        self.methods_hierarchical.append(method_str)
                    t = time.time()
                    cis = bts.ci(coverages=self.alphas, side='one', method=method, quantile_type=self.quantile_type,
                                 sampling='hierarchical')
                    if method_str not in self.computed_intervals:
                        self.computed_intervals[method_str] = defaultdict(list)
                    for a, ci in zip(self.alphas, cis):
                        self.computed_intervals[method_str][a].append(ci)
                    self.times[method_str].append(time.time() - t + ts)
                btss.append(bts)
        return btss

    def compute_non_bootstrap_intervals(self, data: np.array):
        """
        Computes CI with non-bootstrap methods, that can be applied to the statistic in use.
        :param data: array containing one sample
        """
        ci = defaultdict(list)
        new_methods = {'mean': ['wilcoxon', 'ttest'],
                       'median': ['wilcoxon', 'ci_quant_param', 'ci_quant_nonparam', 'maritz-jarrett'],
                       'std': ['chi_sq'], 'percentile': ['ci_quant_param', 'ci_quant_nonparam', 'maritz-jarrett'],
                       'corr': ['ci_corr_pearson']}
        if self.statistic.__name__[:10] not in ['mean', 'median', 'std', 'percentile', 'corr']:
            print(f'No known non-bootstrap methods to use for statistic {self.statistic.__name__}.')
            new_methods[self.statistic.__name__] = []

        for method in new_methods[self.statistic.__name__[:10]]:
            if method == 'ttest':
                t = time.time()
                stat = self.statistic(data)
                se = np.std(data) / np.sqrt(self.n)
                ci[method] = scipy.stats.t.ppf(self.alphas, df=self.n - 1, loc=stat, scale=se)
                self.times[method].append(time.time() - t)

            elif method == 'wilcoxon':
                t = time.time()
                for a in self.alphas:
                    f = robjects.r('''f <- function(data, a) {
                                        res <- wilcox.test(data, conf.int = T, conf.level = a, alternative='less')
                                        b <- lapply(res, attributes)$conf.int$conf.level[1]
                                        list(ci=res$conf.int[2], cl=b)
                                        }
                                   ''')
                    r_f = robjects.globalenv['f']
                    ci_r, achieved_a = r_f(data, a)

                    if abs(achieved_a[0] - a) > a/2:        # R criteria for saying that conf level is not achievable
                        ci_r[0] = np.nan
                    ci[method].append(ci_r[0])
                self.times[method].append(time.time() - t)

            elif method in ['ci_quant_param', 'ci_quant_nonparam']:
                t = time.time()
                # 0.5 for median and 5 or 95 for percentiles:
                quant = 0.5 if self.statistic.__name__ == 'median' else int(self.statistic.__name__.split('_')[1])/100

                if method == 'ci_quant_param':
                    m = np.mean(data)
                    s = np.std(data)
                    z = (np.quantile(data, quant, method=self.quantile_type) - m) / s
                    nc = -z * np.sqrt(self.n)
                    ci[method] = m - scipy.stats.nct.ppf(1 - self.alphas, nc=nc, df=self.n - 1) * s / np.sqrt(self.n)

                elif method == 'ci_quant_nonparam':
                    t = time.time()
                    sorted_data = np.array(sorted(data) + [np.nan])
                    ci[method] = sorted_data[scipy.stats.binom.ppf(self.alphas, self.n, quant).astype(int)]

                self.times[method].append(time.time() - t)

            elif method == 'maritz-jarrett':
                t = time.time()
                quant = 0.5 if self.statistic.__name__ == 'median' else int(self.statistic.__name__.split('_')[1])/100
                ci[method] = [scipy.stats.mstats.mquantiles_cimj(data, prob=quant, alpha=abs(2*a-1))[int(a > 0.5)][0]
                              for a in self.alphas]
                self.times[method].append(time.time() - t)

            elif method == 'chi_sq':
                t = time.time()
                s = np.std(data)
                qchisq = scipy.stats.chi2.ppf(1 - self.alphas, self.n - 1)
                ci[method] = np.sqrt((self.n - 1) * s ** 2 / qchisq)
                self.times[method].append(time.time() - t)

            elif method == 'ci_corr_pearson':
                t = time.time()
                res = scipy.stats.pearsonr(data[:, 0], data[:, 1], alternative='less')
                ci[method] = [res.confidence_interval(a).high for a in self.alphas]
                self.times[method].append(time.time() - t)

        for m in ci.keys():
            if m not in self.methods:
                self.methods.append(m)
            if m not in self.computed_intervals:
                self.computed_intervals[m] = defaultdict(list, {a: [c] for a, c in zip(self.alphas, ci[m])})
            else:
                for i in range(len(self.alphas)):
                    self.computed_intervals[m][self.alphas[i]].append(ci[m][i])

    def exact_interval_simulation(self, repetitions: int):
        true_val = self.dgp.get_true_value(self.statistic.__name__)

        if np.size(true_val) == 1:
            stat_values = np.empty(repetitions)
        else:
            stat_values = np.empty((repetitions, np.size(true_val)))

        if self.sampling == 'hierarchical':
            data = self.dgp.sample(sample_size=self.n, nr_samples=repetitions, group_sizes=self.group_sizes)
        else:
            data = self.dgp.sample(sample_size=self.n, nr_samples=repetitions)

        for r in range(repetitions):
            stat_values[r] = self.statistic(data[r])

        distribution = stat_values - true_val

        # inverse cdf that is used for exact interval calculation
        self.inverse_cdf = {a: np.quantile(distribution, 1 - a, method=self.quantile_type) for a in self.alphas}

    def compare_intervals(self, repetitions, length=None):
        true_statistic_value = self.dgp.get_true_value(self.statistic.__name__)
        self.exact_interval_simulation(100000)

        stat_original = []

        if self.sampling == 'hierarchical':
            data = self.dgp.sample(sample_size=self.n, nr_samples=repetitions, group_sizes=self.group_sizes)
            group_indices = self.dgp.group_indices
        else:
            data = self.dgp.sample(sample_size=self.n, nr_samples=repetitions)

        slow_methods = []                   # hack for speeding up in double or studentized method
        # calculation with different bootstrap methods
        for r in range(repetitions):
            if self.sampling != 'hierarchical':     # only compare strategies of hierarchical sampling among themselves
                # calculation with non-bootstrap methods
                self.compute_non_bootstrap_intervals(data[r])
            # calculation with different bootstrap methods
            if r == 10000:
                slow_methods = [m for m in self.methods if m in ['double', 'studentized']]
                print(f"Leaving out methods {slow_methods} during repetitions over 10000.")
                self.methods = [m for m in self.methods if m not in ['double', 'studentized']]

            if self.sampling == 'hierarchical':
                btss = self.compute_bootstrap_intervals_hierarchical(data[r], group_indices)
                stat_original.append(btss[0].original_statistic_value)      # original doesn't depend on sampling
            else:
                bts = self.compute_bootstrap_intervals(data[r])
                stat_original.append(bts.original_statistic_value)
        stat_original = np.array(stat_original)
        self.methods = self.methods + slow_methods

        # exact intervals
        self.computed_intervals['exact'] = {a: stat_original - self.inverse_cdf[a] for a in self.alphas}

        # compute coverages for all methods and alphas
        methods_to_compare = self.methods_hierarchical if self.sampling == 'hierarchical' else self.methods
        self.coverages = {method: {alpha: np.mean(np.array(self.computed_intervals[method][alpha][-repetitions:]) >
                                                  true_statistic_value) for alpha in self.alphas}
                          for method in methods_to_compare}

        self.distances_from_exact = {method: {alpha: np.array(self.computed_intervals[method][alpha][-repetitions:]) -
                                                     self.computed_intervals['exact'][alpha] for alpha in self.alphas}
                                     for method in methods_to_compare}

        if length is not None:
            low_alpha, high_alpha = [round((1 - length) / 2, 5), round((length + 1) / 2, 5)]
            if low_alpha not in self.alphas or high_alpha not in self.alphas:
                raise ValueError(f"Length of {length} CI can't be calculated, because we don't have calculations for "
                                 f"the corresponding alphas.")
        else:
            low_alpha, high_alpha = [min(self.alphas), max(self.alphas)]
            print(f'Calculating lengths of {round(high_alpha - low_alpha, 5)} CI, '
                  f'from {low_alpha} to {high_alpha} quantiles.')

        self.lengths = {method: np.array(self.computed_intervals[method][high_alpha][-repetitions:]) -
                                np.array(self.computed_intervals[method][low_alpha][-repetitions:])
                        for method in methods_to_compare}

        times_stats = {method: {'mean': np.mean(self.times[method]), 'std': np.std(self.times[method])}
                       for method in methods_to_compare}

        return self.coverages, times_stats

    def draw_intervals(self, alphas_to_draw: list[float], show=False):
        # only implemented for nonparametric sampling, we don't need it for hierarchical
        data = self.dgp.sample(sample_size=self.n)
        self.alphas = np.union1d(self.alphas, alphas_to_draw)

        # compute bootstrap intervals
        bts = self.compute_bootstrap_intervals(data)

        # compute non-bootstrap intervals
        self.compute_non_bootstrap_intervals(data)

        # exact intervals calculation
        if not np.array([a in self.inverse_cdf for a in alphas_to_draw]).all():
            self.exact_interval_simulation(10000)
        exact_intervals = [self.statistic(data) - self.inverse_cdf[a] for a in alphas_to_draw]

        # plotting
        colors = iter(plt.cm.jet(np.linspace(0.05, 0.95, len(self.methods))))
        plt.hist(bts.statistic_values, bins=30, label='statistic')
        if 'smoothed' in self.methods:
            if np.nan in bts.statistic_values_noise or np.inf in bts.statistic_values_noise:
                print('skipped drawing of smoothed values because of nan values.')
            else:
                plt.hist(bts.statistic_values_noise, bins=30, label='smoothed stat.', alpha=0.3)
        for method in self.methods:
            col = next(colors)
            for alpha in alphas_to_draw:
                if alpha == alphas_to_draw[0]:  # label only the first line of a method to avoid duplicates in legend
                    plt.axvline(self.computed_intervals[method][alpha][-1], linestyle='--', label=method, color=col,
                                alpha=0.75)
                else:
                    plt.axvline(self.computed_intervals[method][alpha][-1], linestyle='-.', color=col, alpha=0.75)

        # draw exact intervals
        for e in exact_intervals:
            if e == exact_intervals[0]:
                plt.axvline(e, linestyle=':', label='exact', color='black', alpha=0.75)
            else:
                plt.axvline(e, linestyle=':', color='black', alpha=0.75)

        ci = round((alphas_to_draw[1] - alphas_to_draw[0]) * 100)
        plt.title(f'{ci}CI for {self.statistic.__name__} of {type(self.dgp).__name__} (n = {self.n}, B = {self.b})')
        plt.legend()
        if show:
            plt.show()
        else:
            plt.savefig(f'images/{ci}CI_{self.statistic.__name__}_{type(self.dgp).__name__}_n{self.n}_B{self.b}.png')
            plt.close()

    def combine_results(self, repetitions, length=0.9):
        res = self.compare_intervals(repetitions, length)

        # building long dataframes
        methods_to_use = self.methods_hierarchical if self.sampling == 'hierarchical' else self.methods
        cov_methods = list(methods_to_use)*len(self.alphas)
        cov_alphas = np.repeat(self.alphas, len(methods_to_use))
        coverage_df = pd.DataFrame({'method': cov_methods, 'alpha': cov_alphas,
                                    'coverage': [self.coverages[m][a] for m, a in zip(cov_methods, cov_alphas)]})

        if self.sampling == 'hierarchical':
            var_coverage = [np.mean((np.array(self.variances['_'.join(m.split('_')[:2])]['95CI_l'][-repetitions:]) <
                                     self.gt_variance) & (self.gt_variance <
                                     np.array(self.variances['_'.join(m.split('_')[:2])]['95CI_h'][-repetitions:])))
                            for m in cov_methods]
            coverage_df['var_coverage'] = var_coverage

        df_distance = pd.DataFrame({'alpha': np.repeat(self.alphas, repetitions)} |
                                   {m: np.concatenate([self.distances_from_exact[m][a][-repetitions:]
                                                       for a in self.alphas]) for m in methods_to_use})

        df_length = pd.DataFrame({m: v[-repetitions:] for m, v in zip(self.lengths.keys(), self.lengths.values())})
        df_times = pd.DataFrame({m: v[-repetitions:] for m, v in zip(self.times.keys(), self.times.values())})

        # times are for all alphas combined

        true_val = self.dgp.get_true_value(self.statistic.__name__)

        if self.sampling == 'hierarchical':
            # saving long table for hierarchical
            df_intervals = pd.DataFrame([{'method': m, 'alpha': a, 'predicted': v, 'true_value': true_val}
                                         for m in self.computed_intervals.keys() for a in
                                         self.computed_intervals[m].keys()
                                         for v in self.computed_intervals[m][a]])

            df_intervals['gt_variance'] = self.gt_variance
            # np.nan to skip methods for which variances are not calculated (i.e. exact)
            nn = numpy.empty(repetitions)
            nn[:] = np.nan
            df_intervals['mean_var'] = np.concatenate([self.variances['_'.join(m.split('_')[:2])]['mean']
                                                      if m in self.methods_hierarchical else nn
                                                      for m in self.computed_intervals.keys()
                                                      for _ in self.computed_intervals[m].keys()])
            df_intervals['95CI_var_l'] = np.concatenate([self.variances['_'.join(m.split('_')[:2])]['95CI_l']
                                                        if m in self.methods_hierarchical else nn
                                                        for m in self.computed_intervals.keys()
                                                        for _ in self.computed_intervals[m].keys()])
            df_intervals['95CI_var_h'] = np.concatenate([self.variances['_'.join(m.split('_')[:2])]['95CI_h']
                                                        if m in self.methods_hierarchical else nn
                                                        for m in self.computed_intervals.keys()
                                                        for _ in self.computed_intervals[m].keys()])

        else:
            # saving wide table for nonparametric sampling because of memory issues
            df_intervals = pd.DataFrame([{'alpha': a, 'true_value': true_val} |
                                         {m: self.computed_intervals[m][a][i] for m in self.computed_intervals.keys()}
                                         for a in self.alphas for i in range(repetitions)])

        return res, coverage_df, df_length, df_times, df_distance, df_intervals


def compute_hierarchical_var(bts):
    values = bts.original_sample[bts.bootstrap_indices]

    mean = np.mean(values, axis=1, keepdims=True)
    errors = values - mean

    def flatten(indices):
        if isinstance(indices[0], int):
            return indices
        else:
            return flatten(list(itertools.chain.from_iterable(indices)))

    indices = bts.group_indices.copy()
    random_effects = []
    group_sizes = []

    # calculation of stds without sigma_e
    while isinstance(indices[0], list):
        lvl_indices = [flatten(g) for g in indices]
        lvl_predictors = [np.mean(errors[:, ind], axis=1) for ind in lvl_indices]  # groups x b array
        lvl_std = np.std(np.array(lvl_predictors), axis=0)
        group_sizes.append([len(flatten(g)) for g in indices])

        for pred, ind in zip(lvl_predictors, lvl_indices):
            errors[:, ind] -= np.expand_dims(pred, axis=1)  # leave only residuals of next level in errors

        random_effects.append(lvl_std)
        indices = list(itertools.chain.from_iterable(indices))

    # calculate E[Var]
    evar = 0
    for gs, re in zip(group_sizes, random_effects):
        evar += (bts.n**2 - sum(np.array(gs)**2)) / bts.n**2 * re**2

    # add sigma_e part
    evar += (bts.n - 1) / bts.n * np.std(errors, axis=1) ** 2

    evar2 = np.var(values, axis=1)

    return evar             # mean/95CI ??


def run_comparison(dgps, statistics, Bs, methods, alphas, repetitions, ns=None, alphas_to_draw=[0.05, 0.95],
                   leaves=None, branches=None,
                   length=0.9, append=True, nr_processes=24, dont_repeat=False, sampling='nonparametric'):
    if ns is None:
        ns = [1]        # just to hold space in for loop
    names = ['coverage', 'length', 'times', 'distance', 'intervals']
    if sampling == 'nonparametric':
        all_methods = ['percentile', 'basic', 'bca', 'bc', 'standard',  'smoothed', 'double', 'studentized', 'ttest',
                       'wilcoxon', 'ci_quant_param', 'ci_quant_nonparam', 'maritz-jarrett', 'chi_sq', 'ci_corr_pearson',
                       'double-std']
    else:
        sampling_methods = np.concatenate([['cases' + ''.join(map(str, strat))
                                            for strat in itertools.product([0, 1], repeat=n_lvls)]
                                           for n_lvls in np.unique([len(dgp.stds) for dgp in dgps])])
        ci_methods = ['percentile', 'bca']
        all_methods = [s + '_' + c for s in sampling_methods for c in ci_methods]

    if sampling == 'hierarchical':
        cols = {'coverage': ['method', 'alpha', 'coverage', 'dgp', 'statistic', 'n_leaves', 'n_branches', 'n', 'B',
                             'repetitions', 'balance', 'std', 'levels', 'var_coverage'],
                'length': ['CI', 'dgp', 'statistic', 'n_leaves', 'n_branches', 'n', 'B', 'repetitions', 'balance',
                           'std', 'levels'] + all_methods,
                'distance': ['alpha', 'dgp', 'statistic', 'n_leaves', 'n_branches', 'n', 'B',
                             'repetitions', 'balance', 'std', 'levels'] + all_methods,
                'times': ['dgp', 'statistic', 'n_leaves', 'n_branches', 'n', 'B', 'repetitions', 'balance', 'std',
                          'levels'] + all_methods,
                'intervals': ['method', 'alpha', 'predicted', 'true_value', 'dgp', 'statistic', 'n_leaves',
                              'n_branches', 'n', 'B', 'repetitions', 'balance', 'std', 'levels',
                              'gt_variance', 'mean_var', '95CI_var_l', '95CI_var_h']}

    else:
        cols = {'coverage': ['method', 'alpha', 'coverage', 'dgp', 'statistic', 'n', 'B', 'repetitions'],
                'length': ['CI', 'dgp', 'statistic', 'n', 'B', 'repetitions'] + all_methods,
                'distance': ['alpha', 'dgp', 'statistic', 'n', 'B', 'repetitions'] + all_methods,
                'times': ['dgp', 'statistic', 'n', 'B', 'repetitions'] + all_methods,
                'intervals': ['alpha', 'dgp', 'statistic', 'n', 'B', 'repetitions', 'true_value',
                              'exact'] + all_methods}

    folder = 'results_hierarchical' if sampling == 'hierarchical' else 'results'

    if not append:
        # need to use this (append=False) for running first time to set header!!
        print(f'Will delete all results in folder {folder} - ARE YOU SURE???')
        time.sleep(6)
        for name in names:
            pd.DataFrame(columns=cols[name]).to_csv(f'{folder}/' + name + '.csv', index=False)

    if dont_repeat:
        cov = pd.read_csv(f'{folder}/coverage.csv')

    params = []
    nr_skipped = 0
    for dgp in dgps:
        for statistic in statistics:
            if (statistic.__name__ == 'corr' and type(dgp).__name__ != 'DGPBiNorm') or \
                    (type(dgp).__name__ == 'DGPBiNorm' and statistic.__name__ != 'corr') or \
                    (type(dgp).__name__ == 'DGPCategorical' and statistic.__name__ == 'std'):
                continue

            for n in ns:

                for B in Bs:

                    if dont_repeat:
                        if sampling == 'nonparametric':
                            same_pars = cov[(cov['dgp'] == dgp.describe()) & (cov['statistic'] == statistic.__name__) &
                                            (cov['n'] == n) & (cov['B'] == B) & (cov['repetitions'] == repetitions)]

                            if same_pars.shape[0] > 0:
                                # TODO: check if we have calculations for all methods and alphas if needed, group_sizes
                                nr_skipped += 1
                                continue
                            else:
                                print('Adding: ', dgp.describe(), statistic.__name__, n, B, repetitions)

                    if sampling == 'hierarchical':
                        n_lvls = len(dgp.stds)
                        sampling_parameters = {'strategies': list(itertools.product([0, 1], repeat=n_lvls)),
                                               'sampling_methods': ['cases']}

                        for n_leaves in leaves:
                            for n_branches in branches:
                                if n_leaves == 32 and n_lvls == 4:
                                    continue
                                same_pars = cov[
                                    (cov['dgp'] == dgp.describe()) & (cov['statistic'] == statistic.__name__) &
                                    (cov['n_leaves'] == n_leaves) & (cov['n_branches'] == n_branches) & (cov['B'] == B)
                                    & (cov['repetitions'] == repetitions) & (cov['balance'] == 'balanced')
                                    & (cov['std'] == dgp.stds[0]) & (cov['levels'] == n_lvls)]
                                if same_pars.shape[0] > 0:
                                    nr_skipped += 1
                                    continue
                                else:
                                    print('Adding: ', dgp.describe(), statistic.__name__, n_leaves, n_branches, B,
                                          repetitions, 'balanced', dgp.stds[0], n_lvls)

                                group_sizes = [n_leaves for _ in range(n_branches)]
                                n = n_leaves * n_branches

                                for _ in range(n_lvls - 2):
                                    group_sizes = [group_sizes for _ in range(n_branches)]
                                    n *= n_branches

                                params.append(((statistic, methods.copy(), dgp, n, B, alphas.copy()),
                                               (group_sizes, sampling_parameters, n_leaves, n_branches, 'balanced')))

                    else:
                        params.append((statistic, methods.copy(), dgp, n, B, alphas.copy()))

    if dont_repeat:
        print(f'Skipped {nr_skipped} combinations, because we already have results.')

    pool = Pool(processes=nr_processes)
    for dfs in tqdm(pool.imap_unordered(multiprocess_run_function, zip(params, repeat(repetitions), repeat(length),
                                                                       repeat(alphas_to_draw), repeat(sampling)),
                                        chunksize=1), total=len(params)):
        for i in range(len(dfs)):

            if repetitions >= 10000 and names[i] == 'distance':
                continue    # don't save distances to save space when making many repetitions

            dfs[i] = pd.concat([pd.DataFrame(columns=cols[names[i]]), dfs[i]])      # setting right order of columns

            dfs[i].to_csv(f'{folder}/{names[i]}.csv', header=False, mode='a', index=False)


def multiprocess_run_function(param_tuple):
    pars_tup, repetitions, length, alphas_to_draw, sampling = param_tuple
    if sampling == 'hierarchical':
        pars, samp_pars = pars_tup
        group_sizes, sampling_parameters, n_leaves, n_branches, balance = samp_pars
    else:
        pars = pars_tup
        group_sizes = sampling_parameters = None

    statistic, methods, dgp, n, B, alphas = pars

    use_jit = (repetitions >= 100)
    comparison = CompareIntervals(*pars, use_jit=use_jit, sampling=sampling, group_sizes=group_sizes,
                                  sampling_parameters=sampling_parameters)
    _, coverage_df, df_length, df_times, df_distance, df_intervals = comparison.combine_results(repetitions=repetitions,
                                                                                                length=length)
    dfs = [coverage_df, df_length, df_times, df_distance, df_intervals]

    df_length['CI'] = length
    for i in range(len(dfs)):
        if sampling == 'hierarchical':
            dfs[i]['n_leaves'] = n_leaves
            dfs[i]['n_branches'] = n_branches
            dfs[i]['balance'] = balance
            dfs[i]['std'] = dgp.stds[0]         # TODO change if we'll have different on different levels
            dfs[i]['levels'] = len(dgp.stds)

        dfs[i]['n'] = n
        dfs[i]['dgp'] = dgp.describe()
        dfs[i]['statistic'] = statistic.__name__
        dfs[i]['B'] = B
        dfs[i]['repetitions'] = repetitions

    return dfs


def corr(data):
    c = np.corrcoef(data, rowvar=False)
    return c[0, 1]


def percentile_5(data):
    return np.quantile(data, 0.05, method='median_unbiased')


def percentile_95(data):
    return np.quantile(data, 0.95, method='median_unbiased')


if __name__ == '__main__':

    # non-parametric experiments
    seed = 0
    alphas = [0.025, 0.05, 0.25, 0.75, 0.95, 0.975]

    # methods = ['percentile', 'basic', 'bca', 'bc', 'standard', 'smoothed', 'double', 'studentized']
    methods = ['standard', 'double-std']

    dgps = [DGPNorm(seed, 0, 1), DGPExp(seed, 1), DGPBeta(seed, 1, 1), DGPBeta(seed, 10, 2),
            DGPLaplace(seed, 0, 1), DGPLogNorm(seed, 0, 1),
            DGPBiNorm(seed, np.array([1, 1]), np.array([[2, 0.5], [0.5, 1]]))]
    statistics = [np.mean, np.median, np.std, percentile_5, percentile_95, corr]

    # ns = [4, 8, 16, 32, 64, 128, 256]
    ns = [4, 8, 16, 32]
    Bs = [1000]

    np.random.seed(seed)
    repetitions = 10000
    run_comparison(dgps, statistics, Bs, methods, alphas, repetitions, ns, nr_processes=24, dont_repeat=True,
                   append=True, sampling='nonparametric')

    # # hierarchical experiments
    # leaves = [2, 4, 8, 16, 32]
    # branches = [1, 3, 5, 7]
    # stds = [0.1, 1, 10]
    # levels = [2, 3, 4]
    # dgps = [DGPRandEff(seed, 0, [s for l in range(n_lvl)]) for n_lvl in levels for s in stds]
    #
    # run_comparison(dgps, statistics, Bs, methods, alphas, repetitions, leaves=leaves, branches=branches,
    #                nr_processes=24, dont_repeat=True, append=False, sampling='hierarchical')


