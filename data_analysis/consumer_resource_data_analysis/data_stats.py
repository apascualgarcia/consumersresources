from scipy.stats import mannwhitneyu

def pvalue(sampleA, sampleB):
    (res, pvalue) = mannwhitneyu(sampleA, sampleB, alternative='greater')
    return pvalue
