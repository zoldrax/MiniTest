resultdir = "results"
homedir = ""
dbsdir = homedir+"/DBs"
REF = dbsdir+"/ucsc.GRCh38.fasta"
REFdict = dbsdir+"/ucsc.GRCh38.dict"
CADDdir = dbsdir+"/CADD_1.4_.GRCh38"

panel = "Minitest_pon.vcf.gz"
gnomad = dbsdir+"/gnomad.vcf.bgz"
clinvar = dbsdir+"/clinvar.vcf.gz"

thr=40

rule Mutect2:
    input:
        norm=resultdir+"/{norm_sample}_recal.bam",
        tum=resultdir+"/{tum_sample}_recal.bam",
        norm_idx=resultdir+"/{norm_sample}_recal.bam.bai",
        tum_idx=resultdir+"/{tum_sample}_recal.bam.bai",
        bed=resultdir+"/{tum_sample}.bed",
    output:
        vcf=(resultdir+"/{tum_sample}_vsn_{norm_sample}_rawsomatic.vcf.gz"),
        tbi=(resultdir+"/{tum_sample}_vsn_{norm_sample}_rawsomatic.vcf.gz.tbi"), 
        stats=(resultdir+"/{tum_sample}_vsn_{norm_sample}_rawsomatic.vcf.gz.stats")
    threads: thr
    shell:
        """
          gatk Mutect2 \
                --native-pair-hmm-threads {threads} \
                --dont-use-soft-clipped-bases true \
                --max-reads-per-alignment-start 10000 \
                -R {REF} \
                -L {input.bed} -ip 20 \
                -I {input.tum} \
                -I {input.norm} \
                -normal {wildcards.norm_sample} \
                -O {output.vcf}
         """


rule Mutect2_panel:
    input:
        tum=resultdir+"/{tum_sample}_recal.bam",
        tum_idx=resultdir+"/{tum_sample}_recal.bam.bai",
        bed=resultdir+"/{tum_sample}.bed",
    output:
        vcf=(resultdir+"/{tum_sample}_vsn_panel_rawsomatic.vcf.gz"),
        tbi=(resultdir+"/{tum_sample}_vsn_panel_rawsomatic.vcf.gz.tbi"),
        stats=(resultdir+"/{tum_sample}_vsn_panel_rawsomatic.vcf.gz.stats")
    threads: thr
    shell:
        """
          gatk Mutect2 \
                --native-pair-hmm-threads {threads} \
                --dont-use-soft-clipped-bases true \
                --max-reads-per-alignment-start 10000 \
                -R {REF} \
                -I {input.tum} \
                -L {input.bed} -ip 20 \
                --panel-of-normals {panel} \
                -O {output.vcf}
         """

rule Mutect2_nonorm:
    input:
        tum=resultdir+"/{tum_sample}_recal.bam",
        tum_idx=resultdir+"/{tum_sample}_recal.bam.bai",
        bed=resultdir+"/{tum_sample}.bed",
    output:
        vcf=(resultdir+"/{tum_sample}_vsn_nonorm_rawsomatic.vcf.gz"),
        tbi=(resultdir+"/{tum_sample}_vsn_nonorm_rawsomatic.vcf.gz.tbi"),
        stats=(resultdir+"/{tum_sample}_vsn_nonorm_rawsomatic.vcf.gz.stats")
    threads: thr
    shell:
        """
          gatk Mutect2 \
                --native-pair-hmm-threads {threads} \
                --dont-use-soft-clipped-bases true \
                --max-reads-per-alignment-start 10000 \
                -R {REF} \
                -I {input.tum} \
                -L {input.bed} -ip 20 \
                -O {output.vcf}
         """



rule FilterMutectCalls:
    input:
        vcf=resultdir+"/{tum_sample}_vsn_{norm_sample}_rawsomatic.vcf.gz",
        tbi=resultdir+"/{tum_sample}_vsn_{norm_sample}_rawsomatic.vcf.gz.tbi",
 #       stats=temp(resultdir+"/{tum_sample}_vsn_{norm_sample}_rawsomatic.vcf.gz.stats")
    output:
        vcf=resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic.vcf.gz",
        tbi=resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic.vcf.gz.tbi"
    threads: thr
    shell:
        """
           gatk FilterMutectCalls \
                 -R {REF} \
                 --min-reads-per-strand 2 \
                 -unique 2 \
                 -V {input.vcf} \
                 -O {output.vcf}
        """

rule SnpEff:
    input:
        resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic.vcf.gz"
    output:
        temp(resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic_ann.vcf.gz")
    threads: thr
    shell:
        """
          snpEff -noStats -v -cancer hg38 {input} | \
             SnpSift annotate {gnomad} | \
             SnpSift annotate {clinvar} | \
             bcftools view -Oz -o {output}
        """


rule CADDadd_and_filter:
    input: 
        resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic_ann.vcf.gz"
    output:
        resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic_annCADD.vcf.gz"
    threads: thr
    shell:
        """
           bcftools annotate --rename-chrs {dbsdir}/chrremove_hg19.tab {input} | \
              bcftools annotate -c CHROM,POS,REF,ALT,,PHRED \
                  --threads {threads} \
                  -h {CADDdir}/CADD.hdr -a {CADDdir}/whole_genome_SNVs.tsv.gz | \
              bcftools annotate --rename-chrs {dbsdir}/chradd_hg19.tab -Oz -o {output}
        """



rule OutTable:
    input:
        vcf=resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic_annCADD.vcf.gz"
    output:
        resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic.tsv"
    threads: thr
    shell:
        """
          bcftools view -i '(%FILTER!="haplotype")&&(%FILTER!="map_qual")&&(%FILTER!="normal_artifact")&&(%FILTER!="weak_evidence")&&(%FILTER!="strand_bias")&&(%FILTER!="strict_strand")&&(%FILTER!="slippage")&&(%FILTER!="germline")' -s {wildcards.tum_sample} {input.vcf} | bcftools view -i '(F1R2[0:1]>0)&&(F2R1[0:1]>0)' | \
          SnpSift extractFields -s ";" - CHROM POS REF ALT FILTER ANN[0].IMPACT ANN[0].EFFECT PHRED ANN[0].GENE ANN[0].HGVS_P ID AF_popmax CLNSIG CLNDN GEN[0].F1R2[1] GEN[0].F2R1[1] GEN[*].AD ANN[*].GENE ANN[*].HGVS_P ANN[*].EFFECT > {output}
        """


rule FullCoveredRegions:
    input:
        norm=resultdir+"/{norm_sample}_recal.bam",
        tum=resultdir+"/{tum_sample}_recal.bam",
        norm_idx=resultdir+"/{norm_sample}_recal.bam.bai",
        tum_idx=resultdir+"/{tum_sample}_recal.bam.bai",
        bed=resultdir+"/{tum_sample}.bed"
    output:
        resultdir+"/{tum_sample}_vsn_{norm_sample}_covered.bed"
    threads: 2
    shell:
        """
           samtools depth {input.tum} {input.norm} -b {input.bed} -q 10 -Q 10 |\
           awk 'BEGIN {{FS=OFS="\t"}} ($3>9)&&($4>9) {{print $1,$2-1,$2}}' |\
           bedtools merge > {output}
        """




rule CalculationTMB:
    input:
        vcf=resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic.vcf.gz",
        tbi=resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic.vcf.gz.tbi",
        bed=resultdir+"/{tum_sample}_vsn_{norm_sample}_covered.bed",
        tsv=resultdir+"/{tum_sample}_vsn_{norm_sample}_somatic.tsv"
    output:
        resultdir+"/{tum_sample}_vsn_{norm_sample}_TMB.txt"
    threads: 2
    shell:
        """
           paste \
                   <(bcftools view -f "PASS" -s {wildcards.tum_sample} -R {input.bed} {input.vcf} | bcftools view -i '(F1R2[0:1]>2)&&(F2R1[0:1]>2)' -H | wc -l) \
                   <(bedtools merge -i {input.bed} | \
                              awk -F '\t' 'BEGIN {{s=0}} {{s+=$3-$2}} END {{print s}}') | \
               awk -F '\t' '{{print "TMB = "$1*1000000/$2" mut/MB = ("$1"/"$2/1000000")"}}' > {output}
        """

rule index_bam:
    input:
        "{name1}.bam"
    output:
        "{name1}.bam.bai"
    threads: thr
    resources:
        mem_gb=8
    shell:
        "samtools index -@ {threads} {input} {output}"

