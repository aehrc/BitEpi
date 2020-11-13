import sys
import time
import numpy as np
import pandas as pd
# not compatible with pandas 1.1.0
# pip3 install --user -U pandas==1.0.5
# see https://github.com/limix/pandas-plink/issues/18
from pandas_plink import read_plink

import scipy
import scipy.stats

# Set the random seed
seedMax = 2**32 - 1
seed = int(time.time()*10000000) % seedMax
np.random.seed(seed)

# https://www.hackdeploy.com/fitting-probability-distributions-with-python/


class Distribution(object):

    def __init__(self, dist_names_list=[]):
        self.dist_names = ['gamma', 'norm', 'lognorm', 'expon']
        self.dist_results = []
        self.params = {}

        self.DistributionName = ""
        self.PValue = 0
        self.Param = None

        self.isFitted = False

    def Fit(self, y):
        self.dist_results = []
        self.params = {}
        for dist_name in self.dist_names:
            dist = getattr(scipy.stats, dist_name)
            param = dist.fit(y)

            self.params[dist_name] = param
            # Applying the Kolmogorov-Smirnov test
            D, p = scipy.stats.kstest(y, dist_name, args=param)
            self.dist_results.append((dist_name, p))

        # select the best fitted distribution
        sel_dist, p = (max(self.dist_results, key=lambda item: item[1]))
        # store the name of the best fit and its p value
        self.DistributionName = sel_dist
        self.PValue = p

        self.isFitted = True
        return self.DistributionName, self.PValue

    def Pvalue(self, obs):
        if self.isFitted:
            dist_name = self.DistributionName
            param = self.params[dist_name]
            # initiate the scipy distribution
            dist = getattr(scipy.stats, dist_name)
            return (1 - dist.cdf(obs, *param[:-2], loc=param[-2], scale=param[-1]))
        else:
            raise ValueError('Must first run the Fit method.')

    # fit the best distribution and compute pvalue for the observation
    def FitPval(self, y, obs, statName):
        self.dist_results = []
        self.params = {}
        for dist_name in self.dist_names:
            dist = getattr(scipy.stats, dist_name)
            param = dist.fit(y)

            self.params[dist_name] = param
            # Applying the Kolmogorov-Smirnov test
            D, p = scipy.stats.kstest(y, dist_name, args=param)
            self.dist_results.append((dist_name, p))

        # select the best fitted distribution
        sel_dist, p = (max(self.dist_results, key=lambda item: item[1]))
        # store the name of the best fit and its p value
        self.DistributionName = sel_dist
        self.PValue = p

        self.isFitted = True

        # calculate pvalue for the observation
        dist_name = self.DistributionName
        param = self.params[dist_name]
        # initiate the scipy distribution
        dist = getattr(scipy.stats, dist_name)
        PValueObs = (
            1 - dist.cdf(obs, *param[:-2], loc=param[-2], scale=param[-1]))
        return {statName+'-dist': self.DistributionName, statName+'-distPval': self.PValue, statName+'-pval': PValueObs}

# Return: A genotype datafream generated from plink bfile (BDF: Bfile DataFrame)


def BDF(prefix):
    (bim, fam, bed) = read_plink(prefix, verbose=True)
    bdf = pd.DataFrame(bed.compute().astype('int8')).join(bim[['snp']])    .set_index(
        'snp').append(fam.trait.astype('int8')).transpose().astype('category')
    bdf['cnt'] = 1
    return bdf

# Return the Contingency Table (CT) with sum (s) and row weight (w) given a number of SNPs and Bfile Data Frame


def CT(bdf, SNPs, pheno='trait'):
    ct = bdf.groupby([pheno]+SNPs).count()[['cnt']]
    ctrl = ct.loc[ct.index.get_level_values(pheno) == 1].droplevel(
        level=0).rename(columns={"cnt": "ctrl"})
    case = ct.loc[ct.index.get_level_values(pheno) == 2].droplevel(
        level=0).rename(columns={"cnt": "case"})
    ctx = ctrl.join(case)
    ctx['s'] = ctx.sum(axis=1)
    total = ctx['s'].sum()
    ctx['w'] = ctx['s']/total
    return ctx.fillna(0)

# Return Weighted Average Purity (WAP) given a contingency table


def WAP(ct):
    return ((((ct.ctrl/ct.s)**2) + ((ct.case/ct.s)**2)) * ct.w).sum()

# Return Mximum Lower Order WAP (MLOWAP). Lower Order means to exclude 1 SNP from combination


def MLOWAP(ct, SNPs):
    if len(SNPs)<2:
        return WAP(pd.DataFrame(ct.sum()).transpose())
    lowaps = list()
    for i in range(0, len(SNPs)):
        lo = SNPs.copy()
        del lo[i]
        ctx = ct.groupby(lo).sum()
        lowaps.append(WAP(ctx))
    return max(lowaps)

# return Alpha and Beta (BitEpi) given contingency table and list of SNPs


def AB(ct, SNPs):
    wap = WAP(ct)
    mlowap = MLOWAP(ct, SNPs)
    return {'beta': wap, 'alpha': (wap - mlowap)}


# Compute Alpha and Beta statistics as well as corresponding pvalue for a given interactive SNPs


def PvalueCnt(bdf, SNPs, numRepeat):
    ct = CT(bdf, SNPs)
    ab = AB(ct, SNPs)
    # rab for AB with random phenotype
    rab = list()
    for i in range(numRepeat):

        bdf['RandomPheno'] = np.random.permutation(bdf['trait'])
        ct = CT(bdf, SNPs, pheno='RandomPheno')
        rab.append(AB(ct, SNPs))

    ac = len([r for r in rab if r['alpha'] > ab['alpha']])
    bc = len([r for r in rab if r['beta'] > ab['beta']])

    pv = {'beta': bc/numRepeat, 'alpha': ac/numRepeat}
    return (ab, pv)


def PvalueDist(bdf, SNPs, numRepeat):
    ct = CT(bdf, SNPs)
    ab = AB(ct, SNPs)
    # rab for AB with random phenotype
    rab = list()
    for i in range(numRepeat):
        if(i % 100 == 99):
            print(">>> ", i, "permutation done")

        bdf['RandomPheno'] = np.random.permutation(bdf['trait'])
        ct = CT(bdf, SNPs, pheno='RandomPheno')
        rab.append(AB(ct, SNPs))

    dst = Distribution()

    a = dst.FitPval(list(map(lambda r: r['alpha'], rab)), ab['alpha'], 'alpha')
    b = dst.FitPval(list(map(lambda r: r['beta'], rab)), ab['beta'], 'beta')

    for d in (a, b):
        ab.update(d)

    return ab


def Pvalue(bdf, SNPs, numRepeat):
    return PvalueDist(bdf, SNPs, numRepeat)
    # return PvalueCnt(bdf, SNPs, numRepeat)

# Compute pvalue for list of interaction (only consider top numInteraction in the file)


def EpiPvalue(bfilePrefix, epiFile, numInteraction, numRepeat):

    # Read Plink bfile into pandas dataframe
    bdf = BDF(bfilePrefix)

    # Read BitEpi output into pandas dataframe
    epiInt = pd.read_csv(epiFile)

    # compute pvalue for the top interactions
    result = list()
    for i in range(numInteraction):
        print("Interaction:", i)
        row = epiInt.iloc[i]
        SNPs = list(row[1:].values)
        rd = {'firstCol': row[0], 'SNPs': "#".join(SNPs)}
        pvalue = Pvalue(bdf, SNPs, numRepeat)
        rd.update(pvalue)
        result.append(rd)

    return result

# Compute pvalue of random combination of SNPs


def RndPvalue(bfilePrefix, numSNPs, numInteraction, numRepeat):

    # Read Plink bfile into pandas dataframe
    bdf = BDF(bfilePrefix)
    varId = bdf.columns.values[0:-2]

    # compute pvalue for random interactions
    result = list()
    for i in range(numInteraction):
        print("Interaction:", i)
        np.random.shuffle(varId)
        SNPs = list(varId[0:numSNPs])
        rd = {'firstCol': 0, 'SNPs': "#".join(SNPs)}
        pvalue = Pvalue(bdf, SNPs, numRepeat)
        rd.update(pvalue)
        result.append(rd)

    return result


if __name__ == '__main__':

    argc = len(sys.argv)
    if(argc < 7):
        print("Enter a command and path to the inptu file")
        exit()

    command = sys.argv[1]
    bfilePrefix = sys.argv[2]
    numInteraction = int(sys.argv[3])
    numRepeat = int(sys.argv[4])
    outputFile = sys.argv[5]
    if(command == "rnd"):
        order = int(sys.argv[6])
        pvals = RndPvalue(bfilePrefix, order, numInteraction, numRepeat)
    if(command == "epi"):
        epiFile = sys.argv[6]
        pvals = EpiPvalue(bfilePrefix, epiFile, numInteraction, numRepeat)

    pd.DataFrame(pvals).to_csv(outputFile, sep='\t', index=None)

    # pd.DataFrame(list(map(lambda x: {'firstCol': x['firstCol'], 'SNPs': x['SNPs'],
    #                                  'beta': x['stat']['beta'], 'alpha': x['stat']['alpha'],
    #                                  'pval-beta': x['pvalue']['beta'], 'pval-alpha': x['pvalue']['alpha'], },
    #                       pvals))).to_csv(outputFile, sep='\t', index=None)
