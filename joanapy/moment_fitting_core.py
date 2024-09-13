import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import count
import numpy as np
from scipy.stats import beta, kstest, anderson_ksamp
from joanapy.cramervanmises import testCramerVonMises2Sample
import pandas as pd
from joanapy.utils import write_readme_file
import seaborn as sns


class MOMENT_FITTING:


    def __init__(self, data, resultpath, steps=1000, tolerance=1E-4, components=["falling", "unimodal", "rising"], init_fct=None):
        self.data = data
        self.steps = steps
        self.tolerance = tolerance
        self.components = components
        self.resultpath = resultpath
        self.init_fct = init_fct

    def __get_values(self, x, left, right):
        y = x[np.logical_and(x >= left, x <= right)]
        n = len(y)
        if n == 0:
            m = (left + right) / 2.0
            v = (right - left) / 12.0
        else:
            m = np.mean(y)
            v = np.var(y)
            if v == 0.0:
                v = (right - left) / (12.0 * (n + 1))
        return m, v, n

    def __get_initialization(self, x, ncomponents, limit=0.8):
        # TODO: work with specific components instead of just their number
        points = np.linspace(0.0, 1.0, ncomponents + 2)
        means = np.zeros(ncomponents)
        variances = np.zeros(ncomponents)
        pi = np.zeros(ncomponents)
        # init first component
        means[0], variances[0], pi[0] = self.__get_values(x, points[0], points[1])
        # init intermediate components
        N = ncomponents - 1
        for j in range(1, N):
            means[j], variances[j], pi[j] = self.__get_values(x, points[j], points[j + 2])
        # init last component
        means[N], variances[N], pi[N] = self.__get_values(x, points[N + 1], points[N + 2])

        # compute parameters ab, pi
        ab = [self.__ab_from_mv(m, v) for (m, v) in zip(means, variances)]
        pi = pi / pi.sum()

        # adjust first and last
        if ab[0][0] >= limit:  ab[0] = (limit, ab[0][1])
        if ab[-1][1] >= limit:  ab[-1] = (ab[-1][0], limit)
        return ab, pi

    def __get_initialization_joana(self, x):
        ab = list()
        pi = list()

        # ab.append([1., 100.])
        # ab.append([100., 1.])
        # ab.append([100., 1.])
        # pi.append(0.33)
        # pi.append(0.33)
        # pi.append(1.-0.66)

        mean, variance, temp_pi = self.__get_values(x, 0., 0.1)
        pi.append(temp_pi)
        ab.append(self.__ab_from_mv(mean, variance))
        mean, variance, temp_pi = self.__get_values(x, 0.1, 0.9)
        pi.append(temp_pi)
        ab.append(self.__ab_from_mv(mean, variance))
        mean, variance, temp_pi = self.__get_values(x, 0.9, 1.)
        pi.append(temp_pi)
        ab.append(self.__ab_from_mv(mean, variance))
        
        pi = np.asarray(pi)
        pi = pi / pi.sum()

        return ab, pi


    def __estimate_mixture(self, x, init, steps=1000, tolerance=1E-5, verbose=False):
        """
        estimate a beta mixture model from the given data x
        with the given number of components and component types
        """
        (ab, pi) = init
        n, ncomponents = len(x), len(ab)

        for step in count():
            if step >= steps:
                break
            abold = list(ab)
            piold = pi[:]
            # E-step: compute component memberships for each x
            w = self.__get_weights(x, ab, pi)
            # compute component means and variances and parameters
            for j in range(ncomponents):
                wj = w[:, j]
                pij = np.sum(wj)
                m = np.dot(wj, x) / pij
                v = np.dot(wj, (x - m) ** 2) / pij
                if np.isnan(m) or np.isnan(v):
                    m = 0.5;
                    v = 1 / 12  # uniform
                    ab[j] = (1, 1)  # uniform
                    assert pij == 0.0
                else:
                    assert np.isfinite(m) and np.isfinite(v), (j, m, v, pij)
                    ab[j] = self.__ab_from_mv(m, v)
                pi[j] = pij / n
            delta = self.__get_delta(ab, abold, pi, piold)
            if delta < tolerance:
                break
            if verbose:
                if step % 100 == 0:
                    print((step, delta,
                           ("{:.5f}".format(round(ab[0][0], 5)),
                           "{:.5f}".format(round(ab[0][1], 5)),
                           "{:.5f}".format(round(pi[0], 5))),
                           ("{:.5f}".format(round(ab[1][0], 5)),
                           "{:.5f}".format(round(ab[1][1], 5)),
                           "{:.5f}".format(round(pi[1], 5))),
                           ("{:.5f}".format(round(ab[2][0], 5)),
                           "{:.5f}".format(round(ab[2][1], 5)),
                           "{:.5f}".format(round(pi[2], 5)))))
        usedsteps = step + 1
        return (ab, pi, usedsteps)

    def __estimate_mixture_uniform(self, x, init, steps=1000, tolerance=1E-5, verbose=False):
        """
        estimate a beta mixture model from the given data x
        with the given number of components and component types
        """
        (ab, pi) = init
        n, ncomponents = len(x), len(ab)

        for step in count():
            if step >= steps:
                break
            abold = list(ab)
            piold = pi[:]
            # E-step: compute component memberships for each x
            w = self.__get_weights(x, ab, pi)
            # compute component means and variances and parameters
            for j in range(ncomponents):
                wj = w[:, j]
                pij = np.sum(wj)
                m = np.dot(wj, x) / pij
                v = np.dot(wj, (x - m) ** 2) / pij
                if np.isnan(m) or np.isnan(v) or j == 1:
                    m = 0.5;
                    v = 1 / 12  # uniform
                    ab[j] = (1, 1)  # uniform
                    # assert pij == 0.0
                else:
                    assert np.isfinite(m) and np.isfinite(v), (j, m, v, pij)
                    ab[j] = self.__ab_from_mv(m, v)
                pi[j] = pij / n
            delta = self.__get_delta(ab, abold, pi, piold)
            if delta < tolerance:
                break
            if verbose:
                if step % 100 == 0:
                    print((step, delta,
                           ("{:.5f}".format(round(ab[0][0], 5)),
                           "{:.5f}".format(round(ab[0][1], 5)),
                           "{:.5f}".format(round(pi[0], 5))),
                           ("{:.5f}".format(round(ab[1][0], 5)),
                           "{:.5f}".format(round(ab[1][1], 5)),
                           "{:.5f}".format(round(pi[1], 5))),
                           ("{:.5f}".format(round(ab[2][0], 5)),
                           "{:.5f}".format(round(ab[2][1], 5)),
                           "{:.5f}".format(round(pi[2], 5)))))
        usedsteps = step + 1
        return (ab, pi, usedsteps)

    def __get_weights(self, x, ab, pi):
        """return nsamples X ncomponents matrix with association weights"""
        bpdf = beta.pdf
        n, c = len(x), len(ab)
        y = np.zeros((n, c), dtype=float)
        s = np.zeros((n, 1), dtype=float)
        for (j, p, (a, b)) in zip(count(), pi, ab):
            y[:, j] = p * bpdf(x, a, b)
        s = np.sum(y, 1).reshape((n, 1))
        np.seterr(divide='ignore', invalid='ignore')
        w = y / s  # this may produce inf or nan; this is o.k.!
        # clean up weights w, remove infs, nans, etc.
        wfirst = np.array([1] + [0] * (c - 1), dtype=float)
        wlast = np.array([0] * (c - 1) + [1], dtype=float)
        bad = (~np.isfinite(w)).any(axis=1)
        badfirst = np.logical_and(bad, x < 0.5)
        badlast = np.logical_and(bad, x >= 0.5)
        w[badfirst, :] = wfirst
        w[badlast, :] = wlast
        # now all weights are valid finite values and sum to 1 for each row
        assert np.all(np.isfinite(w)), (w, np.isfinite(w))
        assert np.allclose(np.sum(w, 1), 1.0), np.max(np.abs(np.sum(w, 1) - 1.0))
        return w

    def __get_delta(self, ab, abold, pi, piold):
        epi = max(self.__relerror(p, po) for (p, po) in zip(pi, piold))
        ea = max(self.__relerror(a, ao) for (a, _), (ao, _) in zip(ab, abold))
        eb = max(self.__relerror(b, bo) for (_, b), (_, bo) in zip(ab, abold))
        return max(epi, ea, eb)

    def loglikelihood(self):
        beta_active = beta(self.ab[0][0], self.ab[0][1])
        beta_inactive_1 = beta(self.ab[1][0], self.ab[1][1])
        beta_inactive_2 = beta(self.ab[2][0], self.ab[2][1])
        
        llh = 0
        for x in self.data:
            llh += np.log(self.pi[0]*beta_active.pdf(x) + self.pi[1]*beta_inactive_1.pdf(x) + self.pi[2]*beta_inactive_2.pdf(x))
            
        return llh

    def __relerror(self, x, y):
        if x == y:  return 0.0
        return abs(x - y) / max(abs(x), abs(y))

    def __ab_from_mv(self, m, v):
        """
        estimate beta parameters (a,b) from given mean and variance;
        return (a,b).

        Note, for uniform distribution on [0,1], (m,v)=(0.5,1/12)
        """
        phi = m * (1 - m) / v - 1  # z = 2 for uniform distribution
        return (phi * m, phi * (1 - m))  # a = b = 1 for uniform distribution

    def __estimate(self, x, components, steps=1000, tolerance=1E-4, verbose=False, init='moment', second_comp_uniform=False):
        if self.init_fct is None:
            if init == 'moment':
                init = self.__get_initialization(x, len(components))
            elif init == 'joana':
                init = self.__get_initialization_joana(x)
        else:
            init = self.init_fct()
       
        print(('init:',
                ("{:.5f}".format(round(init[0][0][0], 5)),
                "{:.5f}".format(round(init[0][0][1], 5)),
                "{:.5f}".format(round(init[1][0], 5))),
                ("{:.5f}".format(round(init[0][1][0], 5)),
                "{:.5f}".format(round(init[0][1][1], 5)),
                "{:.5f}".format(round(init[1][1], 5))),
                ("{:.5f}".format(round(init[0][2][0], 5)),
                "{:.5f}".format(round(init[0][2][1], 5)),
                "{:.5f}".format(round(init[1][2], 5)))))

        if second_comp_uniform:
            (ab, pi, usedsteps) = self.__estimate_mixture_uniform(x, init, steps=steps, tolerance=tolerance, verbose=verbose)
        else:
            (ab, pi, usedsteps) = self.__estimate_mixture(x, init, steps=steps, tolerance=tolerance, verbose=verbose)
        return (ab, pi, usedsteps)

    def run(self, verbose=False, init='moment', second_comp_uniform=False):
        (ab, pi, us) = self.__estimate(self.data, self.components, steps=self.steps, tolerance=self.tolerance, verbose=verbose, init=init, second_comp_uniform=second_comp_uniform)
        self.w = self.__get_weights(self.data, ab, pi)
        self.assignments = np.argmax(self.w, axis=1)
        self.ab = ab
        self.pi = pi
        # self.__get_assignments()        

        if verbose:
            print("ab", end="\t")
            print(*[f"{j},{k}" for j, k in ab], sep="\t")

            print("pi", end="\t")
            print(*pi, sep=",")

        # write fit to file
        list_fit = list()
        for i in range(len(self.components)):
            list_fit.append(ab[i][0])
            list_fit.append(ab[i][1])
            list_fit.append(pi[i])

        if not self.resultpath is None:
            pd.DataFrame(list_fit).to_csv(self.resultpath, header=False, index=False)

    def __ecdf(self, sample):

        # convert sample to a numpy array, if it isn't already
        sample = np.atleast_1d(sample)

        # find the unique values and their corresponding counts
        quantiles, counts = np.unique(sample, return_counts=True)

        # take the cumulative sum of the counts and divide by the sample size to
        # get the cumulative probabilities between 0 and 1
        cumprob = np.cumsum(counts).astype(np.double) / sample.size

        return quantiles, cumprob

    # run KS test based on sampled data from the fit
    def goodness_of_fit(self, filename_result, plot_histograms=False):
        # sample data on fit
        samples_fit = list()
        print(range(len(self.w)))
        for i in range(len(self.w)):
            ind_arg_max = np.argmax(self.w[i])
            sample_component = np.random.beta(self.ab[ind_arg_max][0], self.ab[ind_arg_max][1], 1)
            samples_fit.append(sample_component)
        samples_fit = np.concatenate(samples_fit)


        # perform ks test between sampled fit and data
        def __beta_mixture_cdf(x, ab, pi):
            return(pi[0]*beta.cdf(x, ab[0][0], ab[0][1]) +
                   pi[1]*beta.cdf(x, ab[1][0], ab[1][1]) +
                   pi[2]*beta.cdf(x, ab[2][0], ab[2][1]))

        ks_stat, ks_pval = kstest(samples_fit, cdf=__beta_mixture_cdf, args=(self.ab, self.pi), mode='exact')
        self.ks_stat = ks_stat
        self.ks_pval = ks_pval

        # Cramer-von Mises test
        cvm_stat, cvm_pval = testCramerVonMises2Sample(self.data, samples_fit)

        # Anderson Darling test
        ad_stat, ad_critical_values, ad_signiflevel = anderson_ksamp([self.data, samples_fit])

        write_readme_file({'ks_stat' : [ks_stat], 'ks_pval' : [ks_pval],
                           'cvm_stat' : [cvm_stat], 'cvm_pval' : [cvm_pval],
                           'ad_stat' : [ad_stat], 'ad_critical_values' : ad_critical_values, 'ad_signiflevel' : ad_signiflevel},
                          filename_result)

        if plot_histograms:
            plt.figure()
            #plt.plot(
            #    [5, 4, 3], 
            #    [100, 200, 300] 
            #)
            plt.hist([self.data, samples_fit], label=['Data', 'Fit'], bins=20)
            plt.legend(loc='best', framealpha=0.8, fancybox=True, fontsize=18)
            plt.xlabel('q-values', fontsize=18)
            plt.ylabel('Frequency', fontsize=18)
            
            nameFile=filename_result.replace('.txt', '_hist.jpg')
            #plt.gcf().set_size_inches(10, 5)
            plt.tight_layout()
            plt.savefig(fname=nameFile)
            #plt.show()
            plt.close()

            # compute cdf functions of data and samples from the fit
            qe_data, pe_data = self.__ecdf(self.data)
            qe_sample_fit, pe_sample_fit = self.__ecdf(samples_fit)

            plt.figure()
            plt.plot(qe_data, pe_data, '-k', lw=2, label='Data')
            plt.plot(qe_sample_fit, pe_sample_fit, '--r', lw=2, label='Fit')
            plt.legend(loc='upper left', framealpha=0.8, fancybox=True, fontsize=18)
            plt.xlabel('q-values', fontsize=18)
            plt.ylabel('Frequency', fontsize=18)
            plt.tight_layout()
            plt.savefig(filename_result.replace('.txt', '_hist_cdf.jpg'))
            plt.close()

    

    # def __get_assignments(self):
    #     bpdf = beta.pdf
    #     n, c = len(self.data), len(self.ab)
    #     y = np.zeros((n, c), dtype=float)
    #     for (j, (a, b)) in zip(count(), self.ab):
    #         y[:, j] = bpdf(self.data, a, b)
    #     self.assignments = np.argmax(y, axis=1)

    def plot_components_density(self, filename_output):
        components = list()
        for i in range(len(self.assignments)):
            if self.assignments[i] == 0:
                components.append("Active")
            elif self.assignments[i] == 1:
                components.append("Inactive 1")
            elif self.assignments[i] == 2:
                components.append("Inactive 2")

        color_mapping = {"Active" : "green", "Inactive 1" : "coral", "Inactive 2" : "red"}

        kde_data = pd.DataFrame({"q-values" : self.data, "Components" : components})
        kde_data = kde_data.sort_values(by="Components")
        components_unique = kde_data["Components"].unique()
        palette = [color_mapping[x] for x in components_unique]
        plt.figure()
        sns.displot(data=kde_data, x='q-values', hue='Components', kind='kde', fill=True, palette=palette, height=5, aspect=1.5)
        plt.xlabel('q-values', fontsize=18)
        plt.ylabel("Density", fontsize=18)
        plt.xlim([0., 1.])
        plt.savefig(filename_output)
        plt.close()


    def plot_mixture_pdfs(self, filename_output):
        x = np.arange(0., 1., 0.01)
        beta_active = beta(self.ab[0][0], self.ab[0][1])
        beta_inactive_1 = beta(self.ab[1][0], self.ab[1][1])
        beta_inactive_2 = beta(self.ab[2][0], self.ab[2][1])

        pdf_active = beta_active.pdf(x)
        pdf_inactive_1 = beta_inactive_1.pdf(x)
        pdf_inactive_2 = beta_inactive_2.pdf(x)

        plt.figure()
        plt.plot(x, pdf_active, label='Active', color='green')
        plt.plot(x, pdf_inactive_1, label='Inactive 1', color='coral')
        plt.plot(x, pdf_inactive_2, label='Inactive 2', color='red')
        plt.xlabel('x', fontsize=18)
        plt.ylabel("PDF", fontsize=18)
        plt.legend(loc='best', fontsize=18)
        plt.xlim([0., 1.])
        plt.savefig(filename_output, bbox_inches='tight')
        plt.close()


