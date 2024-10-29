# fellow-genomics

## Methylation Deconvolution

Install requires python >= 3.9

within the tissue_decon folder, install the package within a fresh virtual environment:

```
pip install .
```

After the install finishes, you can run the methylation tissue deconvolution script using the command line like:

```
tissue-methyl-decon \
--bedgraph test/test.methyl_CpG.bedGraph \
--panel test/methyl_signature_panel.bed \
--reference test/methyl_sig_10tissue_090924.tsv \
--outfile test/test.methyl_decon.tsv 
```

Test data is available in the test folder within tissue_decon