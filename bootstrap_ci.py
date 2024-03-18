import itertools
import numpy as np
from numba import njit
from scipy.stats import norm
import matplotlib.pyplot as plt


class Bootstrap:

    def __init__(self, data: np.array, statistic: callable, use_jit: bool = False, group_indices: list = None):
        self.original_sample = data
        self.original_statistic_value = statistic(data)
        self.statistic = statistic
        self.n = data.shape[0]
        self.b = 0
        self.bootstrap_indices = np.empty(0)
        self.statistic_values = np.empty(0)
        self.statistic_values_noise = np.empty(0)
        self.implemented_methods = ['basic', 'standard', 'percentile', 'bc', 'bca', 'studentized', 'smoothed', 'double',
                                    'double-std']
        self.use_jit = use_jit
        self.group_indices = group_indices      # hierarchical structure needed for hierarchical bootstrap
        self.bootstrap_values = np.empty(0)     # in random-effect bootstrap we must save the values, not indices

    def sample(self, nr_bootstrap_samples: int = 1000, seed: int = None, sampling: str = 'nonparametric',
               sampling_args: dict = None):
        """
        Draws bootstrap samples from original dataset.
        :param nr_bootstrap_samples: how many samples to draw
        :param seed: random seed
        :param sampling: select type of sampling: choose between nonparametric or hierarchical
        :param sampling_args: sampling arguments, used when doing hierarchical sampling. They should include 'method',
        implemented methods are 'cases' and 'random-effect'. For 'cases' sampling 'strategy' also needs to be defined.
        """
        if seed is not None:
            np.random.seed(seed)
        self.b = nr_bootstrap_samples

        if sampling in ['nonparametric', 'non-parametric']:
            self.bootstrap_indices = np.random.choice(range(self.n), size=[nr_bootstrap_samples, self.n])

        elif sampling == 'hierarchical':
            method = sampling_args['method']            # cases / random-effect

            if method == 'cases':
                groups_n = [[self.group_indices.copy()] for _ in range(nr_bootstrap_samples)]
                strategy = sampling_args['strategy']  # with (True) or without (False) replacement on each level
                for s in strategy:
                    if s:       # with replacement
                        gr_indices_n = [[np.random.randint(len(g), size=len(g)) for g in groups] for groups in groups_n]
                        groups_n = [list(itertools.chain.from_iterable([[g[i] for i in ids]
                                                                        for ids, g in zip(gr_indices, groups)]))
                                    for gr_indices, groups in zip(gr_indices_n, groups_n)]
                    else:       # without replacement
                        groups_n = [list(itertools.chain.from_iterable(groups)) for groups in groups_n]

                self.bootstrap_indices = groups_n  # can't be numpy array as lengths can be different

            elif method == 'random-effect':
                # calculation of predictors on each level
                mean = np.mean(self.original_sample)
                errors = self.original_sample - mean

                def flatten(indices):
                    if isinstance(indices[0], int):
                        return indices
                    else:
                        return flatten(list(itertools.chain.from_iterable(indices)))

                indices = self.group_indices.copy()
                predictors = []
                group_sizes = []
                while isinstance(indices[0], list):
                    lvl_indices = [flatten(g) for g in indices]
                    group_sizes.append([len(g) for g in indices])
                    lvl_predictors = [np.mean(errors[ind]) for ind in lvl_indices]

                    for pred, ind in zip(lvl_predictors, lvl_indices):
                        errors[ind] -= pred         # leave only residuals of next level in errors

                    predictors.append(np.array(lvl_predictors))
                    indices = list(itertools.chain.from_iterable(indices))

                predictors.append(errors)

                # sampling of predictors to build bootstrap samples
                sampled_pred = [pred[np.random.randint(len(pred), size=(self.b, len(pred)))] for pred in predictors]

                self.bootstrap_values = sampled_pred[0]
                for pred, size in zip(sampled_pred[1:], group_sizes):
                    self.bootstrap_values = np.repeat(self.bootstrap_values, size, axis=1) + pred

            else:
                raise ValueError(f'Method {method} of hierarchical sampling is not implemented. Choose between cases'
                                 f'and random-effect.')

        else:
            raise ValueError(f'{sampling} sampling is not implemented. Choose between nonparametric and hierarchical.')

    def evaluate_statistic(self, noise: np.array = None, sampling: str = 'nonparametric', sampling_args: dict = None):
        """
        Evaluates statistic on bootstrapped datasets.
        :param noise: smoothed bootstrap calculates noise separately, so it needs to be included at statistic evaluation
        :param sampling: type of sampling, that tells which values to use at evaluation
        :param sampling_args: needed for sampling method of hierarchical sampling
        """
        if np.size(self.original_statistic_value) == 1:
            self.statistic_values = np.zeros(self.b)
        else:
            self.statistic_values = np.zeros((self.b, np.size(self.original_statistic_value)))

        if noise is not None:
            # save statistic values with noise separately, so we don't override the original ones when calling smoothed
            if np.size(self.original_statistic_value) == 1:
                self.statistic_values_noise = np.zeros(self.b)
            else:
                self.statistic_values_noise = np.zeros((self.b, np.size(self.original_statistic_value)))

            statistic_input_noise = self.original_sample[self.bootstrap_indices]
            statistic_input_noise += noise

        for i in range(self.b):
            if sampling == 'nonparametric' or sampling_args['method'] == 'cases':
                self.statistic_values[i] = self.statistic(self.original_sample[self.bootstrap_indices[i]])
            elif sampling == 'hierarchical' and sampling_args['method'] == 'random-effect':
                self.statistic_values[i] = self.statistic(self.bootstrap_values[i, :])
            if noise is not None:
                self.statistic_values_noise[i] = self.statistic(statistic_input_noise[i, :])

    def ci(self, coverages: np.array, side: str = 'one', method: str = 'double', nr_bootstrap_samples: int = None,
           seed: int = None, sampling: str = 'nonparametric', quantile_type: str = 'median_unbiased',
           sampling_args: dict = {'kernel': 'norm', 'width': None}) -> np.array:
        """
        Returns confidence intervals, calculated with selected method.
        :param coverages: array of coverages for which the values need to be computed. In the case of two-sided
                            intervals (side='two') it is a float number.
        :param side: it is possible to choose between 'one' and 'two' sided confidence intervals. One-sided returns
                        the left-sided confidence interval threshold x, representing CI in the shape of (-inf, x).
        :param method: how to construct confidence intervals. It is possible to select from
                       'percentile', 'basic', 'bca', 'bc', 'standard', 'smoothed', 'double' and 'studentized'
        :param nr_bootstrap_samples: number of bootstrap samples
        :param seed: random seed
        :param sampling: type of sampling, 'nonparametric' or 'hierarchical'
        :param sampling_args: additional arguments used with hierarchical sampling
        :param quantile_type: type of quantiles, possible to select from methods used in numpy's quantile function
        :return: array of values for corresponding coverages
        """

        if nr_bootstrap_samples is not None:        # we will sample again, otherwise reuse previously sampled data
            self.sample(nr_bootstrap_samples, seed, sampling)

        elif len(self.bootstrap_indices) == 0:
            raise ValueError("Choose the number of boostrap samples to do if you haven't sampled before confidence "
                             "interval calculation.")

        if len(self.statistic_values) == 0:
            self.evaluate_statistic(sampling=sampling, sampling_args=sampling_args)
        quantiles = []
        if side == 'two':
            quantiles = np.array([(1-coverages)/2, 0.5 + coverages/2])
        elif side == 'one':
            quantiles = np.array(coverages)
        else:
            raise ValueError("Choose between 'one' and 'two'-sided intervals when setting parameter side.")

        if method == 'percentile':
            return np.quantile(self.statistic_values, quantiles, method=quantile_type)

        elif method == 'standard':
            sd = np.std(self.statistic_values)
            return self.original_statistic_value + sd * norm.ppf(quantiles)

        elif method in ['basic', 'reverse-percentile']:
            tails = np.quantile(self.statistic_values, 1 - quantiles, method=quantile_type)

            return 2 * self.original_statistic_value - tails

        elif method[:2] == 'bc':
            bias = np.mean(self.statistic_values < self.original_statistic_value)

            if bias == 1 or bias == 0:
                e = np.empty(quantiles.shape)
                e[:] = np.nan
                return e

            a = 0   # for BC method
            if method == 'bca':
                jackknife_values = [self.statistic(self.original_sample[np.arange(self.n) != i]) for i in range(self.n)]
                jack_dot = np.mean(jackknife_values)
                u = (self.n - 1) * (jack_dot - np.array(jackknife_values))
                if np.sum(u**2) == 0:
                    u += 1e-10                      # hack for u = 0
                a = np.sum(u**3) / (6 * np.sum(u**2) ** 1.5)

            z_alpha = norm.ppf(quantiles)
            corrected = norm.cdf(norm.ppf(bias) + (norm.ppf(bias) + z_alpha) / (1 - a * (norm.ppf(bias) + z_alpha)))

            return np.quantile(self.statistic_values, corrected, method=quantile_type)

        elif method in ['studentized', 'bootstrap-t', 'bootstrapt']:
            standard_errors = self.studentized_error_calculation()
            t_samples = (self.statistic_values - self.original_statistic_value) / standard_errors
            se = np.std(self.statistic_values)

            return self.original_statistic_value - np.nanquantile(t_samples, 1 - quantiles, method=quantile_type) * se

        elif method == 'smoothed':
            input_shape = self.original_sample[self.bootstrap_indices].shape
            if sampling_args['width'] is None:
                # rule of thumb width selection (we can improve it with AMISE/MISE approximation if needed)
                iqr = np.quantile(self.statistic_values, 0.75, method=quantile_type) - \
                      np.quantile(self.statistic_values, 0.25, method=quantile_type)
                h = 0.9 * min(np.std(self.statistic_values), iqr / 1.34) * (self.n ** -0.2)
            else:
                h = sampling_args['width']
            if sampling_args['kernel'] == 'uniform':
                noise = np.random.uniform(-h, h, input_shape)
            elif sampling_args['kernel'] == 'norm':
                noise = np.random.normal(0, h, input_shape)
            else:
                noise = np.zeros(input_shape)
                print(f"Unknown kernel: {sampling_args['kernel']}, using percentile method.")

            self.evaluate_statistic(noise)
            return np.quantile(self.statistic_values_noise, quantiles, method=quantile_type)

        elif method == 'double':
            # get percentiles of original value in inner bootstrap samples
            nested_btsp_values = self.nested_bootstrap(self.b)
            sample_quantiles = np.mean(nested_btsp_values < self.original_statistic_value, axis=1)

            if side == 'two':
                t = abs(0.5 - sample_quantiles)     # change quantiles to one parameter to get symmetric interval
                t_quantile = np.quantile(t, coverages)
                new_quantiles = [0.5 - t_quantile, 0.5 + t_quantile]
            else:
                new_quantiles = np.quantile(sample_quantiles, quantiles, method=quantile_type)

            return np.quantile(self.statistic_values, new_quantiles, method=quantile_type)

        elif method == 'double-std':
            nested_btsp_values = self.nested_bootstrap(self.b)      # functional values B x B_inner
            stds_inner = np.std(nested_btsp_values, axis=1)
            ahs = norm.cdf(self.original_statistic_value, loc=self.statistic_values, scale=stds_inner)
            calibrated_alphas = np.quantile(ahs, quantiles, method=quantile_type)
            return norm.ppf(calibrated_alphas, loc=self.original_statistic_value, scale=np.std(self.statistic_values))

        else:
            raise ValueError(f'This method is not supported, choose between {self.implemented_methods}.')

    def studentized_error_calculation(self):
        """
        Calculation of standard errors of nested bootstrap values for studentized method.
        :return: standard errors
        """
        # NESTED BOOTSTRAP:
        nested_btsp_values = self.nested_bootstrap(self.b)
        standard_errors = np.std(nested_btsp_values, axis=1)
        if 0 in standard_errors:
            standard_errors = [s if s != 0 else np.nan for s in standard_errors]

        return standard_errors

    def nested_bootstrap(self, b_inner):
        """
        Nested or double bootstrap.
        :param b_inner: number of inner bootstrap samples
        :return: array of size self.b x b_inner containing statistic value of inner bootstrap samples
        """
        if self.b >= 500 or b_inner >= 500 or self.n > 100 or self.use_jit:
            stat_njit = {'mean': wrapped_mean, 'median': wrapped_median, 'std': wrapped_std,
                         'percentile_5': wrapped_percentile_5, 'percentile_95': wrapped_percentile_95,
                         'corr': wrapped_corr}[self.statistic.__name__]
            new_values = nested_bootstrap_jit(self.b, b_inner, self.n, self.original_sample, stat_njit)
        else:
            new_values = np.zeros([self.b, b_inner])
            for i in range(self.b):
                new_indices = np.random.choice(self.bootstrap_indices[i, :], size=[b_inner, self.n])
                for j in range(b_inner):
                    new_values[i, j] = self.statistic(self.original_sample[new_indices[j, :]])

        return new_values

    def plot_bootstrap_distribution(self):
        """Draws distribution of statistic values on all bootstrap samples."""
        plt.hist(self.statistic_values, bins=30)
        plt.show()

    def jackknife_after_bootstrap(self, bootstrap_statistics=[np.mean]):
        """
        Jackknife-after-bootstrap diagnostic to analyse if sampling is too dependent on any single point.
        :param bootstrap_statistics: which statistic to use
        """
        jk_indices = {i: [j for j in range(self.b) if i not in self.bootstrap_indices[j, :]]
                      for i in range(self.n)}   # gives indices of samples where the point is not included

        jk_stat_values = {i: self.statistic_values[jk_indices[i]] for i in range(self.n)}   # statistic on those samples

        jk_values = np.array([np.mean(jk_stat_values[i]) for i in range(self.n)])
        jk_influence = (np.mean(jk_values) - jk_values) * (self.n - 1)
        point_order = np.argsort(jk_influence)

        min_val = self.original_statistic_value
        for bs in bootstrap_statistics:
            jk_values_bs = np.array([bs(jk_stat_values[i]) for i in range(self.n)])

            plt.plot(jk_influence, jk_values_bs, 'o', color='black')
            plt.plot(jk_influence[point_order], jk_values_bs[point_order], color='black')
            plt.axhline(bs(self.statistic_values), linestyle='--', color='black')

            min_val = min(min_val, min(jk_values_bs))

        for point in range(self.n):
            if abs(jk_influence[point]) > 1:
                plt.text(jk_influence[point], min_val - (self.original_statistic_value - min_val) / 10, point)

        plt.show()


# NUMBA functions, used to significantly speed up nested bootstrap calculation
@njit()
def nested_bootstrap_jit(b, b_inner, n, original_sample, statistic):
    """
    Nested/double bootstrap but using library njit to speed up the calculations.
    :param b: number of bootstrap samples
    :param b_inner: number of inner bootstrap samples
    :param n: dataset size
    :param original_sample: original dataset
    :param statistic: statistic in which we are interested
    :return: array of size self.b x b_inner containing statistic value of inner bootstrap samples
    """
    bootstrap_indices = np.random.choice(np.arange(n), size=(b, n))
    new_values = np.zeros((b, b_inner))
    for i in range(b):
        new_indices = np.random.choice(bootstrap_indices[i, :], size=(b_inner, n))
        for j in range(b_inner):
            new_values[i, j] = statistic(original_sample[new_indices[j, :]])

    return new_values


# Wrapped functions to be able to use them with njit library
@njit()
def wrapped_median(val):
    return np.median(val)


@njit()
def wrapped_mean(val):
    return np.mean(val)


@njit()
def wrapped_std(val):
    return np.std(val)


@njit()
def wrapped_corr(data):
    c = np.corrcoef(data, rowvar=False)
    return c[0, 1]


@njit()
def wrapped_percentile_5(data):
    # same as np.quantile(data, 0.05, method='median_unbiased'), written for njit use
    return quantile_median_unbiased(data, 0.05)


@njit()
def wrapped_percentile_95(data):
    # same as np.quantile(data, 0.95, method='median_unbiased'), written for njit use
    return quantile_median_unbiased(data, 0.95)


@njit()
def quantile_median_unbiased(data, q):
    n = len(data)
    m = (q+1)/3         # a = b = 1/3, m = a + p*(1 - a - b)
    j = int(q*n + m)
    g = q*n + m - j
    data = np.sort(data)
    return (1 - g) * data[j-1] + g * data[j]

