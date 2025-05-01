import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.stats import norm
import logging
logger = logging.getLogger("app_logger")

def adjust_score(loadings, population_model, population_var_model):

    X = loadings[['PC1', 'PC2', 'PC3', 'PC4']]
    Xconst = sm.add_constant(X, has_constant='add')

    predldpred = population_model.predict(Xconst)
    residuals = loadings['ldpred'] - predldpred
    predicted_var = population_var_model.predict(Xconst)
    loadings['adjusted_score'] = residuals / np.sqrt(predicted_var)
    loadings['percentile'] = norm.cdf(loadings['adjusted_score'])
    return loadings


def calibrate(prscore, vcf_file, fileroot, chr):
    if chr == '0':
        map1 = list(pd.read_table(f"{fileroot}1000G_map.txt")['ID'])
        center = pd.read_table(f"{fileroot}1000G_center.txt")['out.center'].values
        scale = pd.read_table(f"{fileroot}1000G_scale.txt")['out.scale'].values
        V =  pd.read_table(f"{fileroot}1000G_PC1.txt").values
        D = pd.read_table(f"{fileroot}1000G_lambda.txt")['out.d'].values
    else:
        map1 = list(pd.read_table(f"{fileroot}1000G_map_chr{chr}.txt")['ID'])
        center = pd.read_table(f"{fileroot}1000G_center_chr{chr}.txt")['x'].values
        scale = pd.read_table(f"{fileroot}1000G_scale_chr{chr}.txt")['x'].values
        V =  pd.read_table(f"{fileroot}1000G_PC1_chr{chr}.txt").values
    
    pca = pd.read_table(f"{fileroot}1000G_PCA.txt")
    y = pca['ldpred']
    X = sm.add_constant(pca[['PC1', 'PC2', 'PC3', 'PC4']])
    #population_model = sm.GLM(y, X, family=sm.families.Gaussian()).fit()
    #Residuals
    #pca['residual_score'] = population_model.resid_response
    #squared residuals
    #pca['residual_score2'] = pca['residual_score'] ** 2
    #population_resid_mean = pca['residual_score'].mean()
    #population_resid_sd = pca['residual_score'].std()
    #population_var_model = sm.GLM(pca['residual_score2'], X, family=sm.families.Gaussian()).fit()

    #vcf_file = vcf_file[np.array([id in map1 for id in ids])]
    #Make sure order of dosages, center and scale factors are the same!
    vcf_file = vcf_file.set_index('combined_id').reindex(map1).dropna()
    dosages = vcf_file.apply(
    lambda row: dict(zip(
        row.iloc[8].split(":"),
        row.iloc[9].split(":")
    ))["DS"], axis=1) 
    dosages = np.array(dosages.values, dtype=float)
    if not (len(dosages) == len(center) == len(scale)):
        raise ValueError("Dosage, center, and scale vectors must be of the same length")

    r2 = vcf_file.apply(
        lambda row: float(dict(
            item.split("=") for item in row.iloc[7].split(";") if "=" in item).get("R2", None)), 
        axis=1)
    r2mean = np.mean(r2)
    r2median = np.median(r2)
    #Strand of 1kg SNPs are flipped
    dosages = 2 - dosages
    tdose = (dosages - center ) / scale
    loadings = np.dot(tdose, V)
    #Return loadings instead of adjusted score

    return {
        "chr": chr,
        "prs": prscore,
        "loadings": loadings,
        "r2mean": r2mean,
        "r2median":r2median,
        #"population_model": population_model,
        #"population_var_model": population_var_model
    }

def calc(vcf_file_path, fileroot, chr):
    
    if chr == '0':
        snp_weight = pd.read_table(fileroot + "trans_prs_Nov_19.txt")
    else:
        snp_weight = pd.read_table(fileroot + chr + ".trans_prs_snps.txt")

    wl = set(snp_weight['newid'])
    cn = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    vcf_file_o = pd.read_csv(vcf_file_path,
                        sep='\t',
                        comment='#',
                        engine='c',
                        names=cn,
                        header=None,
                        memory_map=True,
                        compression="gzip")
    #We drop duplicates. Only a few and mostly indels/CNAs. 
    vcf_file_o.drop_duplicates(subset=["ID", "REF", "ALT"], inplace=True)
    newids = vcf_file_o["ID"].astype(str) + ':' + vcf_file_o["REF"].astype(str) + ':'+ vcf_file_o["ALT"].astype(str)

    #vcf_file = vcf_file_o[np.array([newid in wl for newid in newids])]
    vcf_file_o['combined_id'] = newids
    #Ensure that the order of dosages and weights is the same!
    vcf_file = vcf_file_o.set_index('combined_id').reindex(snp_weight['newid']).dropna()
    ds_vals = vcf_file.apply(
    lambda row: dict(zip(
        row.iloc[8].split(":"),
        row.iloc[9].split(":")
    ))["DS"], axis=1)
    ds_vals = np.array(ds_vals.values, dtype=float)
    weight = snp_weight.loc[snp_weight['newid'].isin(newids), 'beta_grid4']
    sum=np.dot(ds_vals, weight.values)
    #Calibration
    caliobj = calibrate(sum, vcf_file_o, fileroot, chr)
    return caliobj
