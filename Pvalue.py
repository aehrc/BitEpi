import sys
import time
import numpy as np
import pandas as pd
# not compatible with pandas 1.1.0
# pip3 install --user -U pandas==1.0.5
# see https://github.com/limix/pandas-plink/issues/18
from pandas_plink import read_plink

# Set the random seed
seedMax = 2**32 - 1
seed = int(time.time()*10000000) % seedMax
np.random.seed(seed)

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


# Rather applying random phenotype with uniform distribution "np.random.choice([1,2], bdf.shape[0])"
# we can use a normal distribution for probability of being case
# performance wise it make it faster
# But "what standard deviation" reflect original unifrom distribution is a question I cannot answer it for now.
rndCT = False  # if true it will use this function to randomize phenotype
sd = 0.03


def RP(ct):
    # pc: probability of being case
    pc = np.random.normal(0.5, sd, ct.shape[0])
    pc[pc > 1] = 1
    pc[pc < 0] = 0
    ct['pc'] = pc
    ct.case = (ct.s * ct.pc).astype('int32')
    ct.ctrl = (ct.s - ct.case).astype('int32')
    ct.drop(['pc'], axis=1, inplace=True)

# Compute Alpha and Beta statistics as well as corresponding pvalue for a given interactive SNPs


def Pvalue(bdf, SNPs, numRepeat=100):
    ct = CT(bdf, SNPs)
    ab = AB(ct, SNPs)
    # rab for AB with random phenotype
    rab = list()
    for i in range(numRepeat):

        if rndCT:
            RP(ct)
        else:
            bdf['RandomPheno'] = np.random.choice([1, 2], bdf.shape[0])
            ct = CT(bdf, SNPs, pheno='RandomPheno')
        rab.append(AB(ct, SNPs))

    ac = len([r for r in rab if r['alpha'] > ab['alpha']])
    bc = len([r for r in rab if r['beta'] > ab['beta']])

    pv = {'beta': bc/numRepeat, 'alpha': ac/numRepeat}
    return (ab, pv)

# Compute pvalue for list of interaction (only consider top numInteraction in the file)


def EpiPvalue(bfilePrefix, epiFile, numInteraction, numRepeat):

    # Read Plink bfile into pandas dataframe
    bdf = BDF(bfilePrefix)

    # Read BitEpi output into pandas dataframe
    epiInt = pd.read_csv(epiFile)

    # compute pvalue for the top interactions
    result = list()
    for i in range(numInteraction):
        row = epiInt.iloc[i]
        firstCol = row[0]
        SNPs = list(row[1:].values)
        (stat, pvalue) = Pvalue(bdf, SNPs, numRepeat)
        result.append({'firstCol': firstCol, 'SNPs': SNPs,
                       'stat': stat, 'pvalue': pvalue})

    return result

# Compute pvalue of random combination of SNPs


def RndPvalue(bfilePrefix, numSNPs, numInteraction, numRepeat):

    # Read Plink bfile into pandas dataframe
    bdf = BDF(bfilePrefix)
    varId = bdf.columns.values[0:-2]

    # compute pvalue for the top interactions
    result = list()
    for i in range(numInteraction):
        firstCol = 0
        np.random.shuffle(varId)
        SNPs = list(varId[0:numSNPs])
        (stat, pvalue) = Pvalue(bdf, SNPs, numRepeat)
        result.append({'firstCol': firstCol, 'SNPs': SNPs,
                       'stat': stat, 'pvalue': pvalue})

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

    pd.DataFrame(list(map(lambda x: {'firstCol': x['firstCol'], 'SNPs': x['SNPs'],
                                     'beta': x['stat']['beta'], 'alpha': x['stat']['alpha'],
                                     'pval-beta': x['pvalue']['beta'], 'pval-alpha': x['pvalue']['alpha'], },
                          pvals))).to_csv(outputFile, sep='\t', index=None)
