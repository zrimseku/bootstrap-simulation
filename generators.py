import itertools

import numpy as np
import scipy.stats


class DGP:

    def __init__(self, seed: int, true_statistics: dict = None):
        np.random.seed(seed)
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics = true_statistics

    def sample(self, sample_size: int, nr_samples: int = 1) -> np.array:
        raise NotImplementedError()

    def get_true_value(self, statistic_name):
        if statistic_name not in self.true_statistics:
            raise ValueError(f"True value of {statistic_name} is not known. You should specify it at DGP "
                             f"initialization.")
        return self.true_statistics[statistic_name]

    def describe(self):
        return type(self).__name__


class DGPNorm(DGP):

    def __init__(self, seed: int, loc: float = 0, scale: float = 1, true_statistics: dict = None):
        super(DGPNorm, self).__init__(seed, true_statistics)
        self.loc = loc
        self.scale = scale
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics['mean'] = loc
        self.true_statistics['median'] = loc
        self.true_statistics['std'] = scale
        self.true_statistics['percentile_5'] = loc - 1.645 * scale
        self.true_statistics['percentile_95'] = loc + 1.645 * scale

    def sample(self, sample_size: int, nr_samples: int = 1) -> np.array:
        size = (nr_samples, sample_size) if nr_samples != 1 else sample_size
        return np.random.normal(self.loc, self.scale, size=size)

    def describe(self):
        return type(self).__name__ + '_' + str(self.loc) + '_' + str(self.scale)


class DGPExp(DGP):

    def __init__(self, seed: int, scale: float = 1, true_statistics: dict = None):
        super(DGPExp, self).__init__(seed, true_statistics)
        self.scale = scale                      # 1/lambda
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics['mean'] = scale
        self.true_statistics['median'] = scale * np.log(2)
        self.true_statistics['std'] = scale
        self.true_statistics['percentile_5'] = scale * np.log(20/19)
        self.true_statistics['percentile_95'] = scale * np.log(20)

    def sample(self, sample_size: int, nr_samples: int = 1) -> np.array:
        size = (nr_samples, sample_size) if nr_samples != 1 else sample_size
        return np.random.exponential(self.scale, size=size)

    def describe(self):
        return type(self).__name__ + '_' + str(self.scale)


class DGPBeta(DGP):

    def __init__(self, seed: int, alpha: float = 1, beta: float = 1, true_statistics: dict = None):
        super(DGPBeta, self).__init__(seed, true_statistics)
        self.alpha = alpha
        self.beta = beta
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics['mean'] = alpha / (alpha + beta)
        self.true_statistics['std'] = np.sqrt(alpha * beta / (alpha + beta) ** 2 / (alpha + beta + 1))
        self.true_statistics['percentile_5'] = scipy.stats.beta.ppf(0.05, alpha, beta)
        self.true_statistics['percentile_95'] = scipy.stats.beta.ppf(0.95, alpha, beta)
        self.true_statistics['median'] = scipy.stats.beta.ppf(0.5, alpha, beta)

    def sample(self, sample_size: int, nr_samples: int = 1) -> np.array:
        size = (nr_samples, sample_size) if nr_samples != 1 else sample_size
        return np.random.beta(self.alpha, self.beta, size=size)

    def describe(self):
        return type(self).__name__ + '_' + str(self.alpha) + '_' + str(self.beta)


class DGPLogNorm(DGP):

    def __init__(self, seed: int, mean: float, sigma: float, true_statistics: dict = None):
        super(DGPLogNorm, self).__init__(seed, true_statistics)
        self.mean = mean
        self.sigma = sigma
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics['mean'] = np.exp(mean + (sigma ** 2) / 2)
        self.true_statistics['median'] = np.exp(mean)
        self.true_statistics['std'] = (np.exp(2 * mean + sigma ** 2) * (np.exp(sigma ** 2) - 1)) ** 0.5
        self.true_statistics['percentile_5'] = np.exp(mean - 1.645 * sigma)
        self.true_statistics['percentile_95'] = np.exp(mean + 1.645 * sigma)

    def sample(self, sample_size: int, nr_samples: int = 1) -> np.array:
        size = (nr_samples, sample_size) if nr_samples != 1 else sample_size
        return np.random.lognormal(self.mean, self.sigma, size=size)

    def describe(self):
        return type(self).__name__ + '_' + str(self.mean) + '_' + str(self.sigma)


class DGPLaplace(DGP):

    def __init__(self, seed: int, loc: float, scale: float, true_statistics: dict = None):
        super(DGPLaplace, self).__init__(seed, true_statistics)
        self.loc = loc
        self.scale = scale
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics['mean'] = loc
        self.true_statistics['median'] = loc
        self.true_statistics['std'] = scale * 2**0.5
        self.true_statistics['percentile_5'] = loc + scale * np.log(0.1)
        self.true_statistics['percentile_95'] = loc - scale * np.log(0.1)

    def sample(self, sample_size: int, nr_samples: int = 1) -> np.array:
        size = (nr_samples, sample_size) if nr_samples != 1 else sample_size
        return np.random.laplace(self.loc, self.scale, size=size)

    def describe(self):
        return type(self).__name__ + '_' + str(self.loc) + '_' + str(self.scale)


class DGPBernoulli(DGP):

    def __init__(self, seed: int, p: float, true_statistics: dict = None):
        super(DGPBernoulli, self).__init__(seed, true_statistics)
        self.p = p
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics['mean'] = p
        if p == 0.5:
            self.true_statistics['median'] = 0.5
        elif p < 0.5:
            self.true_statistics['median'] = 0
        else:
            self.true_statistics['median'] = 1
        self.true_statistics['std'] = (p * (1 - p)) ** 0.5
        self.true_statistics['percentile_5'] = float(p >= 0.95)
        self.true_statistics['percentile_95'] = float(p >= 0.05)

    def sample(self, sample_size: int, nr_samples: int = 1) -> np.array:
        size = (nr_samples, sample_size) if nr_samples != 1 else sample_size
        return np.random.binomial(1, self.p, size=size).astype(float)

    def describe(self):
        return type(self).__name__ + '_' + str(self.p)


class DGPCategorical(DGP):

    def __init__(self, seed: int, pvals: np.array, true_statistics: dict = None):
        super(DGPCategorical, self).__init__(seed, true_statistics)
        self.pvals = pvals
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics['mean'] = np.sum(pvals * np.array(range(len(pvals))))
        self.true_statistics['median'] = np.where(np.cumsum(pvals) > 0.5)[0][0]
        self.true_statistics['percentile_5'] = np.where(np.cumsum(pvals) > 0.05)[0][0]
        self.true_statistics['percentile_95'] = np.where(np.cumsum(pvals) > 0.95)[0][0]

    def sample(self, sample_size: int, nr_samples: int = 1) -> np.array:
        size = (nr_samples, sample_size) if nr_samples != 1 else sample_size
        return np.array([np.argmax(c, axis=-1) for c in np.random.multinomial(1, self.pvals, size=size)]).astype(float)

    def describe(self):
        return type(self).__name__ + '_' + str(self.pvals)


class DGPBiNorm(DGP):

    def __init__(self, seed: int, mean: np.array, cov: np.array, true_statistics: dict = None):
        super(DGPBiNorm, self).__init__(seed, true_statistics)
        self.mean = mean        # means of both variables, 1D array of length 2
        self.cov = cov          # covariance matrix, 2D array (2x2)
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics['mean'] = mean
        self.true_statistics['median'] = mean
        self.true_statistics['std'] = np.diag(cov) ** 0.5
        self.true_statistics['corr'] = cov[0, 1] / (cov[0, 0] * cov[1, 1]) ** 0.5   # extend it to multiple??

    def sample(self, sample_size: int, nr_samples: int = 1) -> np.array:
        size = (nr_samples, sample_size) if nr_samples != 1 else sample_size
        return np.random.multivariate_normal(self.mean, self.cov, size=size)

    def describe(self):
        return type(self).__name__ + '-' + '_'.join([str(par) for par in [self.mean[0], self.mean[1], self.cov[0, 0],
                                                                          self.cov[0, 1], self.cov[1, 1]]])


class DGPRandEff(DGP):

    def __init__(self, seed: 0, mean: float, stds: list, true_statistics: dict = None):
        super(DGPRandEff, self).__init__(seed, true_statistics)
        self.group_indices = []
        self.mean = mean
        self.stds = stds
        if true_statistics is None:
            self.true_statistics = {}
        self.true_statistics['mean'] = mean
        self.true_statistics['median'] = mean
        self.true_statistics['std'] = stds
        self.group_sizes = None

    def sample(self, sample_size: int = None, nr_samples: int = 1, group_sizes: list = None,
               max_group_sizes: list = None) -> np.array:
        if group_sizes is None:
            # we will generate groups on random, based on specified max sizes on each level
            if len(self.stds) != len(max_group_sizes):
                raise ValueError(f'Specified standard deviations of the generator imply different number of levels '
                                 f'({len(self.stds) - 1}) than max_group_sizes ({len(max_group_sizes) - 1})!')
            else:
                def get_sizes(max_sizes):
                    if len(max_sizes) == 1:
                        return np.random.randint(1, max_sizes[0])
                    else:
                        return [get_sizes(max_sizes[1:]) for _ in range(np.random.randint(1, max_sizes[0]))]
                group_sizes = get_sizes(max_group_sizes)

        else:
            # only checking if depths match
            depth = 1
            depth_test = group_sizes.copy()
            while isinstance(depth_test, list):
                depth += 1
                depth_test = depth_test[0]
            if len(self.stds) != depth:
                raise ValueError(f'Specified standard deviations of the generator imply different number of levels '
                                 f'({len(self.stds)}) than group_sizes ({depth})!')

        self.group_sizes = group_sizes

        counter = itertools.count()

        def get_indices(sizes, cnt):
            if isinstance(sizes, int):
                return [next(cnt) for _ in range(sizes)]
            else:
                return [get_indices(sizes[j], cnt) if isinstance(sizes[j], int) else
                        get_indices(sizes[j], cnt) for j in range(len(sizes))]

        self.group_indices = get_indices(group_sizes, counter)

        sample_size_true = next(counter)

        if sample_size_true != sample_size:
            if max_group_sizes is None:   # not best way to check this, but group_sizes is not None anymore
                raise ValueError(f'Actual sample size, obtained from group_sizes ({sample_size_true}) is not the same '
                                 f'as specified sample_size ({sample_size}).')
            else:
                # group sizes are random, so they are probably not the same as specified sample_size. No error.
                if sample_size is not None:
                    print(f'Actual sample size, obtained from group_sizes ({sample_size_true}) is not the same as '
                          f'specified sample_size ({sample_size}).')

        size = (nr_samples, sample_size)
        data = np.zeros(size) + self.mean

        def add_variance(indices, data, depth, varr):
            if isinstance(indices, int):
                data[:, indices] += varr
            else:
                for i in range(len(indices)):
                    # self.stds can be lists on each level if we want to set each groups std differently
                    if self.stds[depth] is None:
                        s = np.random.uniform(0, 1, size=nr_samples)
                    elif isinstance(self.stds[depth], int) or isinstance(self.stds[depth], float):
                        s = self.stds[depth]
                    else:
                        if len(self.stds[depth]) != len(indices):
                            raise ValueError(
                                f'If you set each groups std separately, the length of specified stds on this '
                                f'level ({len(self.stds[depth])}) should coincide with the number of groups '
                                f'of the level ({len(indices)}).')
                        s = self.stds[depth][i]
                    add_variance(indices[i], data, depth + 1, varr + np.random.normal(0, s, size=nr_samples))

        add_variance(self.group_indices, data, 0, 0)

        if nr_samples == 1:
            return data[0]      # removes unnecessary parenthesis
        else:
            return data

    def describe(self):
        return type(self).__name__ + '_' + str(self.mean)

