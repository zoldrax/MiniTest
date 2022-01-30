# MiniTest 179

A next-generation sequencing assay designed to cover 179 cancer-associated genes predicting sensitivity or resistance to
various types of anti-cancer therapy. The genes were selected from Cancer Gene Census based on their involvement in drug
response and include biomarkers from NCCN and ESMO guidelines. The panel is eligible for detecting all classes of
mutations (indels, substitutions, amplifications, and rearrangements) as well as TMB, LOH, and MSI in a single DNA
workflow. MiniTest 179 can be used for blood and FFPE DNA to assess germline and somatic variants.

Citation: [Iyevleva AG, Aleksakhina SN, Sokolenko AP, Baskina SV, Venina AR, Anisimova EI, Bizin IV, Ivantsov AO, Belysheva YV, Chernyakova AP, Togo AV, Imyanitov EN. Somatic loss of the remaining allele occurs approximately in half of CHEK2-driven breast cancers and is accompanied by a border-line increase of chromosomal instability. Breast Cancer Res Treat. 2022 Jan 12. doi: 10.1007/s10549-022-06517-3. Epub ahead of print. PMID: 35020107.][1]

## Installation

Install [Snakemake][2]:

```
pip install snakemake
```

Install [GATK][3]:

```
wget https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip
unzip gatk-4.2.2.0.zip
# export PATH="/path/to/gatk-4.2.2.0:$PATH"
```

Install [SnpEff][4]:

```
wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_3t_core.zip/download
unzip snpEff_v4_3t_core.zip
echo -e '#!/bin/sh\njava -jar -Xmx8G /path/to/snpEff/snpEff.jar $@' > ~/bin/snpEff
chmod a+x ~/bin/snpEff
echo -e '#!/bin/sh\njava -jar -Xmx8G /path/to/snpEff/SnpSift.jar $@' > ~/bin/SnpSift
chmod a+x ~/bin/SnpSift
```

Install [bedtools][5]:

```
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
mv bedtools.static.binary ~/bin/bedtools
```

Install [samtools and bcftools][6]:

```
wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2
tar jvfx samtools-1.14.tar.bz2
tar jvfx bcftools-1.14.tar.bz2
cd samtools-1.14
./configure
make
cp samtools ~/bin/
cd ../bcftools-1.14
./configure
make
cp bcftools ~/bin/
cd ..
```

Get [CADD][7], [GnomAD][8], and [ClinVar][9] databases and set path to them in the Snakefile.

## Usage

Run snakemake in the directory with ```Snakefile```,  ```Tumor.bam```, ```Blood.bam```, and Tumor.bed:

```
snakemake -j40 Tumor_vsn_Blood_TMB.txt
```

The ```Tumor.bed``` is a symbolic link to ```UTTV3_capture_targets_GRCh38.bed```

```
ln -s UTTV3_capture_targets_GRCh38.bed Tumor.bed
```

## Funding

The study has been supported by the [Russian Science Foundation (21-75-30015)][10].

[1]: https://doi.org/10.1007/s10549-022-06517-3

[2]: https://snakemake.github.io

[3]: https://gatk.broadinstitute.org

[4]: https://pcingola.github.io/SnpEff

[5]: https://bedtools.readthedocs.io

[6]: http://www.htslib.org/download

[7]: https://cadd.gs.washington.edu/download

[8]: https://gnomad.broadinstitute.org/downloads

[9]: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/

[10]: https://rscf.ru/en/project/21-75-30015/