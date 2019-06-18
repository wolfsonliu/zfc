# ZFC #

ZFC is a software to calculate fold change zscore of screening data.

ZFC can used for CRISPR library screening with or without [ibar][1], with
or without replicates.

## Dependency ##

ZFC is designed using python3 and require packages that is available
in Linux, Mac OS, and Windows.

* Python3.x
* numpy >= 1.10
* scipy >= 1.0
* pandas >= 0.16
* matplotlib >= 2.0.0

## Installation ##

Clone this repository, and install the software.

```{shell}
$ git clone https://github.com/wolfsonliu/zfc.git
$ cd zfc
$ python3 setup.py install
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
  --punish-rate PUNISH_RATE
                        The punish rate used in the punish of inconsistent
                        barcode reads. The value should better be in [0, 1],
                        the larger the harsher of the punishment, default to
                        be 0.5.
  --n-sd N_SD           The lfc within n-sd sds should be consider similar to
                        0. The value should better be in (0, 2], The smaller
                        the strict of the estimation, default to be 1.2.
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


### Reference ###

[1]: Zhu, S. et al. [Guide RNAs with embedded barcodes boost CRISPR-pooled screens](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1628-0). *Genome Biology* **20**, (2019).
