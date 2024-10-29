import argparse
import sys

import pandas as pd
import pyranges as pr
from pyranges import PyRanges
import numpy as np
from scipy import optimize

class MethylSample:
    """Class to load, subset, and filter methylation data for a single sample"""
    def __init__(self, cpg_bedgraph_file: str,  panel_file: str, min_cov: int = 5, impute: bool = True):
        # the bedgraph file from methyl dackel (6 column format)
        print(cpg_bedgraph_file, panel_file)
        self.cpg_bedgraph_file = cpg_bedgraph_file
        # the panel file (bed format)
        self.panel_file = panel_file
        # minimum coverage for a CpG site to be included
        self.min_cov = min_cov
        # impute missing values at the panel region level
        self.impute = impute
        # load the methyldackyl file
        self.cpg_pr = self.load_cpg_results(self.cpg_bedgraph_file)
        self.filter_cpg_results()
        self.panel_pr = self.load_panel_file(panel_file)
        self.grouped_cpg_df = self.group_cpg_results()

    @staticmethod
    def load_cpg_results(fpath: str) -> PyRanges:
        """Loads methyl dackel output bedgraph"""
        # skip trackline (if exists)
        md_df = pd.read_table(fpath, comment='t', header=None)
        # note that here Percent is an int [0,100] with methyldackel - turn into a proportion [0,1]
        md_df.columns = ["Chromosome", "Start", "End", "Percent", "Meth_bases", "Unmeth_bases"]
        md_df['total_bases'] = md_df['Meth_bases'] + md_df['Unmeth_bases']
        md_df['Percent'] = md_df['Percent'] / 100.
        return PyRanges(md_df)

    @staticmethod
    def load_panel_file(fpath: str) -> PyRanges:
        """Loads a panel file into a pyranges"""
        panel_pr = pr.read_bed(fpath)
        # avoid typing issues with the below assign
        # add a 'key' for each panel region based on location for downstream grouping
        panel_pr = panel_pr.apply(lambda df:  (df.assign(key = lambda df: 
                                                    df['Chromosome'].astype(str) + ':' 
                                                    + df['Start'].apply(str) + '-' + df['End'].apply(str))))
        return panel_pr

    def filter_cpg_results(self):
        """Filters the cpg dataframe"""
        # filter the bp level cpg data to have sufficient coverage
        #self.cpg_pr = self.cpg_pr[self.cpg_pr.total_bases >= self.min_cov]
        pass

    def group_cpg_results(self) -> pd.DataFrame:
        """Groups the CpG results into windowed regions, index is window key"""
        # add the panel region to the cpg counts dataframe, will group by the panel region key
        join_df = self.cpg_pr.join(self.panel_pr, how='right').df
        grouped_df = join_df.groupby(by='key').agg({'Percent': 'mean', 'Meth_bases': 'sum', 'Unmeth_bases': 'sum'})
        # split the key back into genomic coordinates
        grouped_df[['Chromosome', 'Start', 'End']] = grouped_df.apply(
            lambda row: row.name.replace('-', ':').split(':'), result_type='expand', axis=1)
        grouped_df = grouped_df.astype({'Chromosome': 'category', 'Start': 'int', 'End': 'int'})
        grouped_df = grouped_df[["Chromosome", "Start", "End", "Percent", "Meth_bases", "Unmeth_bases"]]
        # the above join will set -1 values as default, convert to nan
        grouped_df[['Percent', 'Meth_bases', 'Unmeth_bases']] =  grouped_df[['Percent', 'Meth_bases', 'Unmeth_bases']].applymap(lambda x: np.nan if x < 0 else x)
        if self.impute:
            print(f'about to impute {grouped_df["Percent"].apply(pd.isnull).sum()} values')
            # impute missing values as the mean 
            grouped_df[['Percent', 'Meth_bases', 'Unmeth_bases']] = \
                grouped_df[['Percent', 'Meth_bases', 'Unmeth_bases']].fillna(grouped_df[['Percent', 'Meth_bases', 'Unmeth_bases']].mean())
        # drop any rows with na values
        grouped_df = grouped_df.dropna(axis=0)
        # note that percent here is the mean of the percentage values across all CpGs in the region,
        # and will not just be the percent methylated of these numbers
        grouped_df = grouped_df.astype({'Percent': 'float', 'Meth_bases': 'int', 'Unmeth_bases': 'int'})
        return grouped_df



class BulkMethylDecon:
    """Deconvolute methylation data into component tissue types"""
    def __init__(self, cpg_bedgraph_file: str, panel_file: str, reference_tissue_file: str):
        self.cpg_bedgraph_file = cpg_bedgraph_file
        self.panel_file = panel_file
        self.reference_tissue_file = reference_tissue_file
        # load and filter the CpG data into windows
        self.windowed_methyl_df = self.load_sample_methylation(self.cpg_bedgraph_file, self.panel_file)
        self.reference_tissue_df = self.load_ref_tissue(self.reference_tissue_file)
        self.tissue_decon_df = self.nnls_decon()

    @staticmethod
    def load_ref_tissue(ref_tissue_file) -> pd.DataFrame:
        """Loads the reference tissue matrix file"""
        return pd.read_table(ref_tissue_file, index_col=0)

    @staticmethod
    def load_sample_methylation(cpg_bedgraph_file, panel_file) -> pd.DataFrame:
        """Load the methylation data for the given sample"""
        # dont impute
        return MethylSample(cpg_bedgraph_file, panel_file, impute=False).grouped_cpg_df

    def common_region_filter(self) -> pd.Series:
        """Returns an index of common regions from the sample and the reference dataset"""
        shared_regions = self.windowed_methyl_df.index.intersection(self.reference_tissue_df.index)
        missing_region_count = len(self.reference_tissue_df) - len(shared_regions)
        print(f'Missing data from {missing_region_count} reference regions')
        return shared_regions

    def nnls_decon(self) -> pd.DataFrame:
        """Use non-negative least squares to deconvolute methylation data by tissue"""
        common_regions = self.common_region_filter()
        sample_query = self.windowed_methyl_df.loc[common_regions, 'Percent'].values
        ref_data = self.reference_tissue_df.loc[common_regions].values
        try:
            estimates, residuals =  optimize.nnls(ref_data, sample_query)
        except RuntimeError:
            print('NNLS failed, runtime error')
            # set to 1 for all tissues (normalized downstream)
            estimates = np.zeros(len(self.reference_tissue_df.columns)) + 1
        # normalize the coefficients to sum to 1
        estimates = estimates / estimates.sum()
        # normalize the coefficients to sum to 1
        estimates = estimates / estimates.sum()
        decon_df = pd.DataFrame([estimates], columns=self.reference_tissue_df.columns)
        return decon_df


def main():
    parser = argparse.ArgumentParser(description='Deconvolute MethylDackel CpG data into tissue types')
    parser.add_argument('-b', '--bedgraph', type=str, required=True, help='Path to the CpG bedgraph file')
    parser.add_argument('-p', '--panel', type=str, required=True, help='Path to the panel file')
    parser.add_argument('-r', '--reference', type=str, required=True, help='Path to the reference matrix')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Path to the output deconvolution file')

    args = parser.parse_args()

    bms = BulkMethylDecon(args.bedgraph, args.panel, args.reference)
    bms.tissue_decon_df.to_csv(args.outfile, sep='\t')


if __name__ == '__main__':
    main()