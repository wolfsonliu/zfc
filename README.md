![ZFC](doc/zfc_card.png)

# ZFC #

ZFC is a software to calculate fold change z score of screening data.

ZFC can used for CRISPR library screening with or without [ibar][1], with
or without replicates. Please cite us [(Xu et al., 2021)][2]:

* Xu, P., Liu, Z., Liu, Y., Ma, H., Xu, Y., Bao, Y., Zhu, S., Cao, Z., Wu, Z., Zhou, Z., et al. (**2021**). Genome-wide interrogation of gene functions through base editor screens empowered by barcoded sgRNAs. *Nat Biotechnol*.


## Dependency ##

ZFC is designed with python3 and requires packages that are available
in Linux, Mac OS, and Windows.

* Python3.x
* numpy >= 1.10
* scipy >= 1.0
* pandas >= 0.16
* matplotlib >= 2.0
* sklearn >= 0.20

## Installation ##

Clone this repository, then install the software.

```{shell}
$ git clone https://github.com/wolfsonliu/zfc.git
$ cd zfc
$ python3 setup.py install
```

Or from pypi:

```{shell}
$ pip install --user zfc
```

It's advised to use
[virtualenv](https://virtualenv.pypa.io/en/stable/) or other software
to create virual environment for the ZFC software.

## Usage ##

The help of ZFC software:

```{shell}
usage: zfc [-h] [-i INPUT] [-o OUTPREFIX] [--punish-rate PUNISH_RATE]
           [--n-sd N_SD] [--null-iteration NULL_ITERATION] [--plot]
           [--version]

Calculate fold change of screening data (zscore log2 fold change).

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Raw count table with header should include: <gene>,
                        <guide>, <barcode>, <ctrl>, <exp>.<ctrl> is the raw
                        counts of control group, <exp> is the raw counts of
                        treatment group. For screening without barcode, the
                        barcode column can be the same with guide.
  -o OUTPREFIX, --outprefix OUTPREFIX
                        Output file prefix, can be the file directory path
                        with file name prefix. The directory in the outprefix
                        should be built before analysis.
  --normalization  {total,median,quantile,median_ratio,none}
                        Normalization of raw count data, default is total.
                        Support method: total (Total sum normalization);
                        median (Median normalization); quantile (Upper
                        quantile normalization (0.75)); median_ratio (Median
                        ratio normalization); none (No normalization). 
  --top-n-sgrna TOP_N_SGRNA
                        Only consider top n barcodes for each sgRNA. Default
                        to use all the data.
  --top-n-gene TOP_N_GENE
                        Only consider top n barcodes for each gene. Default to
                        use all the data.
  --null-iteration NULL_ITERATION
                        The iteration to generate null distribution in
                        calculating the p value of genes. The larger the
                        better, but slower in calculation, default to be 100.
  --plot                Output figures.
  --version             show program's version number and exit
```

## Example ##

* For screening data with replicates but without ibar, you should
  first format the replicates as barcode data and then use zfc to
  calculate.

```{shell}
gene	guide	barcode	ctrl	exp
A1BG	AAGAGCGCCTCGGTCCCAGC	R_A	213	0
A1BG	AAGAGCGCCTCGGTCCCAGC	R_B	213	1.03341
A1BG	AAGAGCGCCTCGGTCCCAGC	R_C	213	49.2844
A1BG	CAAGAGAAAGACCACGAGCA	R_A	647	679.474
A1BG	CAAGAGAAAGACCACGAGCA	R_B	647	295.554
A1BG	CAAGAGAAAGACCACGAGCA	R_C	647	472.941
A1BG	CACCTTCGAGCTGCTGCGCG	R_A	469	190.335
A1BG	CACCTTCGAGCTGCTGCGCG	R_B	469	62.0044
A1BG	CACCTTCGAGCTGCTGCGCG	R_C	469	280.542
A1BG	CACTGGCGCCATCGAGAGCC	R_A	678	188.288
A1BG	CACTGGCGCCATCGAGAGCC	R_B	678	165.345
A1BG	CACTGGCGCCATCGAGAGCC	R_C	678	202.824
A1BG	GCTCGGGCTTGTCCACAGGA	R_A	559	333.597
A1BG	GCTCGGGCTTGTCCACAGGA	R_B	559	103.341
A1BG	GCTCGGGCTTGTCCACAGGA	R_C	559	409.44
A1BG	TGGACTTCCAGCTACGGCGC	R_A	363	176.008
A1BG	TGGACTTCCAGCTACGGCGC	R_B	363	307.955
A1BG	TGGACTTCCAGCTACGGCGC	R_C	363	254.952
```

* For screening data with ibar but without replicates, you can use the
  software directely.

```{shell}
gene	guide	barcode	ctrl	exp
A1BG	ACTTCCAGCTGTTCAAGAAT	CTCGCT	597.0	659.0
A1BG	ACTTCCAGCTGTTCAAGAAT	GATGGT	1038.0	1233.0
A1BG	ACTTCCAGCTGTTCAAGAAT	GCACTG	884.0	855.0
A1BG	ACTTCCAGCTGTTCAAGAAT	TCCACT	807.0	870.0
A1BG	CGAGAGCCAGGGAGCAGGCA	CTCGCT	777.0	948.0
A1BG	CGAGAGCCAGGGAGCAGGCA	GATGGT	1448.0	1385.0
A1BG	CGAGAGCCAGGGAGCAGGCA	GCACTG	1225.0	1205.0
A1BG	CGAGAGCCAGGGAGCAGGCA	TCCACT	1030.0	1196.0
A1BG	GACTTCCAGCTGTTCAAGAA	CTCGCT	448.0	252.0
A1BG	GACTTCCAGCTGTTCAAGAA	GATGGT	685.0	700.0
A1BG	GACTTCCAGCTGTTCAAGAA	GCACTG	487.0	513.0
A1BG	GACTTCCAGCTGTTCAAGAA	TCCACT	383.0	409.0
```

* For screening data with ibar and replicates, you can analysis
  replicates separately, or make the replicates as ibar for analysis.

```{shell}
A1BG	ACTTCCAGCTGTTCAAGAAT	CTCGCT-R_1	402.548	399.666
A1BG	ACTTCCAGCTGTTCAAGAAT	CTCGCT-R_2	399.486	435.624
A1BG	ACTTCCAGCTGTTCAAGAAT	GATGGT-R_1	699.908	738.835
A1BG	ACTTCCAGCTGTTCAAGAAT	GATGGT-R_2	675.699	703.222
A1BG	ACTTCCAGCTGTTCAAGAAT	GCACTG-R_1	596.068	785.816
A1BG	ACTTCCAGCTGTTCAAGAAT	GCACTG-R_2	675.039	655.51
A1BG	ACTTCCAGCTGTTCAAGAAT	TCCACT-R_1	544.148	565.71
A1BG	ACTTCCAGCTGTTCAAGAAT	TCCACT-R_2	588.023	558.014
A1BG	CGAGAGCCAGGGAGCAGGCA	CTCGCT-R_1	523.92	643.584
A1BG	CGAGAGCCAGGGAGCAGGCA	CTCGCT-R_2	605.821	514.451
A1BG	CGAGAGCCAGGGAGCAGGCA	GATGGT-R_1	976.365	1022.66
A1BG	CGAGAGCCAGGGAGCAGGCA	GATGGT-R_2	922.905	1006.78
A1BG	CGAGAGCCAGGGAGCAGGCA	GCACTG-R_1	826	778.093
A1BG	CGAGAGCCAGGGAGCAGGCA	GCACTG-R_2	883.352	845.664
A1BG	CGAGAGCCAGGGAGCAGGCA	TCCACT-R_1	694.514	738.835
A1BG	CGAGAGCCAGGGAGCAGGCA	TCCACT-R_2	777.218	842.898
```

* For screening data without ibar and replicates, you can assign guide
  column as barcode column for analysis, which means each guide RNA
  have only one ibar.

```{shell}
gene	guide	barcode	ctrl	exp
A1BG	CACCTTCGAGCTGCTGCGCG	CACCTTCGAGCTGCTGCGCG	469	125
A1BG	AAGAGCGCCTCGGTCCCAGC	AAGAGCGCCTCGGTCCCAGC	213	0
A1BG	TGGACTTCCAGCTACGGCGC	TGGACTTCCAGCTACGGCGC	363	119
A1BG	CACTGGCGCCATCGAGAGCC	CACTGGCGCCATCGAGAGCC	678	122
A1BG	GCTCGGGCTTGTCCACAGGA	GCTCGGGCTTGTCCACAGGA	559	212
A1BG	CAAGAGAAAGACCACGAGCA	CAAGAGAAAGACCACGAGCA	647	464
A1CF	CGTGGCTATTTGGCATACAC	CGTGGCTATTTGGCATACAC	898	322
A1CF	GGTATACTCTCCTTGCAGCA	GGTATACTCTCCTTGCAGCA	199	94
A1CF	GACATGGTATTGCAGTAGAC	GACATGGTATTGCAGTAGAC	271	118
A1CF	GAGTCATCGAGCAGCTGCCA	GAGTCATCGAGCAGCTGCCA	158	33
A1CF	GGTGCAGCATCCCAACCAGG	GGTGCAGCATCCCAACCAGG	195	25
A1CF	CCAAGCTATATCCTGTGCGC	CCAAGCTATATCCTGTGCGC	353	367
```


## Algorithm ##

The ZFC analysis algorithm adopts the z-score of log2 fold change as
the judgement of the sgRNA and gene changes between reference group
(without treatment) and experiment group (with treatment). ZFC
supports screening with [iBAR][1] employed, as well as conventional
screening with replicates. The sgRNA with replicates and sgRNA-iBAR is
treated with similar procedure.

### Step 1: Normalization of raw counts ###

We use total counts for the normalization of raw counts, to rectify
the batch sequencing deptch. Because some sgRNAs in the reference have
very low raw counts, which can affect the fold change calculation of
the following analysis. We define sgRNAs counts less than 0.05
quantile both in reference group and experiment group as the small
count sgRNAs. The mean counts of all the small count sgNRAs were added
to the normalized counts. The normalized counts is calculated as
following expression:

```{latex}
$$Cn_{i} = \frac{Cr_{i}}{S_{ref}} \times 0.5 \times (S_{ref} + S_{exp})$$
$$C_{i} = Cn_{i} + Cm_{small}$$
```

where, `$Cn_{i}$` is the normalized count of ith sgRNA-iBAR in
reference group, `$Cr_{i}$` is the raw count of ith sgRNA-iBAR,
`$Cm_{small}$` is the mean counts of all the small count sgRNAs,
`$S_{ref}$` is the sum of raw counts of all the sgRNA-iBAR in
reference group, `$S_{exp}$` is the sum of raw counts of all the
sgRNA-iBAR in experiment group, `$C_{i}$` is the final normalized
counts for ith sgRNA-iBAR after small count adjustment. The normalized
counts for sgRNAs in experiment group are calculated similarly.


### Step 2: Calculate fold change ###

The raw fold change of each sgRNA-iBAR is calculated from the
normalized counts of reference and experiment groups.

```{latex}
$$fc_{i} = \frac{Cref_{i}}{Cexp_{i}}$$
$$lfc_{i} = \log_{2}fc_{i}$$
```

where, `$fc_{i}$` is the fold change (FC) of ith sgRNA-iBAR and `$lfc_{i}$`
is the log2 fold change (LFC) of ith sgRNA-iBAR, `$Cref_{i}$` and
`$Cexp_{i}$` are the normalized counts of reference and experiment
groups respectively.

### Step 3: Calculate fold change std ###

In order to calculate z-score of LFC (ZLFC), the standard deviation of
LFC should be calculated. The LFC of sgRNA-iBAR is related to the
normalized counts of reference group. So the standard deviations of
LFC are different for sgRNA-iBARs with different normalized counts of
reference group. All the sgRNA-iBARs are divided into several sets
according to the normalized counts of reference group. And the
standard deviations of log fold change are calculated among the
divided sets. The LFC standard diviation and the normalized counts of
the reference is linearly related. So, linear model is calculated for
the LFC sd and reference counts. And the linear model is used to
calculate the LFC standard diviation for all the sgRNA-iBAR.

### Step 4: Raw z score of log fold change ###

The raw z score of log fold change is calculated.

```{latex}
$$raw ZLFC = \frac{LFC}/{LFC std}$$
```

### Step 5: Calculate sgRNA mean z score of fold change ###

The sgRNA level ZLFCs are calculated as the mean of all the ZLFCs of
the relevant sgRNA$^{iBAR}$s.

```{latex}
$$ZLFC_{sgRNA} = \frac{\sum{ZLFC_{sgRNA-iBAR}}}{n}$$
```

where, the sgRNA has n sgRNA-iBAR.

Empirical P value is also calculated for the sgRNA ZLFCs. The p value
is adjusted considering control of False Discovery Rate.


### Step 6: Calculate gene mean zscore of fold change ###

The gene level ZLFCs are calculated as the mean of all the ZLFCs of
the relevant sgRNAs.

```{latex}
$$ZLFC_{gene} = \frac{\sum{ZLFC_{sgRNA}}}{m}$$
```

where, the gene has m sgRNAs.

Empirical P value is also calculated for the gene ZLFCs. The p value
is adjusted considering control of False Discovery Rate.


### Step 7: Robust rank aggregation analysis ###

[Robust rank aggregation][3] is utilized to calculate the rank
significance of the gene with the sgRNA-iBARs in the whole
library. Aside from robust rank aggregation, mean rank aggregation is
also calculated. The robust rank score is adjusted
considering control of Fault Discovery Rate.

***

[1]: <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1628-0> "Zhu, S. et al. Guide RNAs with embedded barcodes boost CRISPR-pooled screens. Genome Biology 20, (2019)."
[2]: <http://www.nature.com/articles/s41587-021-00944-1> "Xu, P., Liu, Z., Liu, Y., Ma, H., Xu, Y., Bao, Y., Zhu, S., Cao, Z., Wu, Z., Zhou, Z., et al. (2021). Genome-wide interrogation of gene functions through base editor screens empowered by barcoded sgRNAs. Nat Biotechnol."
[3]: <http://bioinformatics.oxfordjournals.org/content/28/4/573.abstract> "Kolde, R., Laur, S., Adler, P. & Vilo, J. Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics 28, 573–580 (2012)."
