import math 
import os

configfile: "config.yaml"

READS       =   "{sample}/diploid/reads.fq.gz"
VCF         =   config["vcf"].replace(".vcf.gz", ".filtered.vcf.gz")
REF_PANEL   =   config["vcf"].replace(".vcf.gz", ".filtered.leave_out.vcf.gz")
DIR         =   "{sample}/diploid/{coverage}x"
snakemake.utils.makedirs(expand(DIR, sample=config["samples"], coverage=config["coverage"]))
GENOME1      =   os.path.join("{sample}/diploid", "genome1.fa")
GENOME2      =   os.path.join("{sample}/diploid", "genome2.fa")
BAM         =   os.path.join(DIR, "hg19_alns.bam") 
COUNTS    =   os.path.join(DIR, "genotype_by_count.vcf.gz")
FILTERED_COUNTS    =   os.path.join(DIR, "genotype_by_count.filtered.vcf.gz")
BEAGLE_GT   =   os.path.join(DIR, "genotype_by_impute.vcf.gz")
PG          =   os.path.join(DIR, "pg.{i}.fa")
LIFT        =   os.path.join(DIR, "pg.{i}.lft")
REF_LENGTHS="data/ref_lengths.txt"

rule all:
    input:
        expand(LIFT, sample=config["samples"], coverage=config["coverage"], i=[1,2]),
        expand(PG, sample=config["samples"], coverage=config["coverage"], i=[1,2])

rule filter_vcf:
    input:
        vcf=config["vcf"]
    output:
        vcf=VCF,
        idx=VCF + ".csi"
    shell:
        "bcftools view -O z -V mnps,other {input.vcf} > {output.vcf};\n"
        "bcftools index {output.vcf}"

rule make_ref_panel:
    input:
        vcf=VCF,
        idx=VCF + ".csi"
    output:
        vcf=REF_PANEL
    params:
        samples="|".join(config["samples"])
    shell:
        "bcftools view -Oz -S <(bcftools query -l {input.vcf} | grep -Pv {params.samples}) {input.vcf} > {output.vcf};\n"
        "bcftools index {output.vcf}\n"

rule index_hg19:
    input:
        config["index"]
    output:
        config["index"] + ".1.bt2",
        config["index"] + ".2.bt2",
        config["index"] + ".3.bt2",
        config["index"] + ".4.bt2",
        config["index"] + ".rev.1.bt2",
        config["index"] + ".rev.2.bt2"
    threads: 16
    shell:
        "bowtie2-build --threads {threads} {input} {input}"

# CHANGES: no simulates BOTH haplotypes!
rule simulate_haplotype:
    input:
        ref=config["index"],
        vcf=VCF,
        vcf_idx=VCF + ".csi"
    output:
        fa1=GENOME1,
        fa2=GENOME2
    params:
        sample = "{sample}"
    shell:
        "bcftools consensus -f {input.ref} -H 1 -s {params.sample} {input.vcf} > {output.fa1};\n"
        "bcftools consensus -f {input.ref} -H 2 -s {params.sample} {input.vcf} > {output.fa2}"

rule simulate_reads:
    input:
        fa1=GENOME1,
        fa2=GENOME2
    output:
        fa=temp("{sample}/diploid/full_genome.fa"),
        fq=READS,
        sam=READS.replace(".fq.gz", ".sam")
    params:
        nreads=config["nreads"],
        rlen=config["read_length"]
    threads: 16
    shell:
        "cat {input.fa1} {input.fa2} | awk 'BEGIN {{i=1}} /^>/ {{print \">genome\"i++}} /^[^>]/ {{print}}' > {output.fa};\n"
        "mason_simulator --num-threads {threads} "
                         "-ir {output.fa} "
                         "-n {params.nreads} "
                         "--illumina-read-length {params.rlen} "
                         "-o {output.fq} "
                         "-oa {output.sam}"

rule align_subset:
    input:
        fq=READS,
        idx1=config["index"] + ".1.bt2"
    params:
        fracs=lambda wildcards: float(wildcards.coverage) / config["total_coverage"], # multiply by 2 for diploid?
        index=config["index"]
    output:
        bam=BAM
    threads: 16
    shell:
        "bowtie2 --threads {threads} -x {params.index} -U <(seqtk sample {input.fq} {params.fracs}) | samtools sort -@ {threads} > {output.bam}"

rule count:
    input:
        bam=BAM,
        vcf=VCF
    output:
        vcf=COUNTS
    params:
        thres=lambda wildcards: math.ceil(float(wildcards.coverage) / 2),
        sample="{sample}"
    shell:
        "varcount/varcount -glikelihood -s {params.sample} {input.vcf} {input.bam} | "
        "bcftools sort -O z > {output.vcf};\n"
        "bcftools index {output.vcf}"


rule filter_counts:
    input:
        vcf=COUNTS
    output:
        vcf=FILTERED_COUNTS,
        tmp=temp(FILTERED_COUNTS + "_tmp.vcf")
    shell:
        "bcftools view -i 'GT~\"A\"' {input.vcf} > {output.tmp};\n"
        "bcftools view -H -i 'GT~\"RR\"' {input.vcf} | shuf -n 50000 >> {output.tmp};\n"
        "bcftools sort -Oz {output.tmp} > {output.vcf};\n"
        "bcftools index {output.vcf};\n"

rule beagle_impute:
    input:
        vcf=FILTERED_COUNTS,
        panel=REF_PANEL,
        gmap=config['beagle_map']
    output:
        vcf=BEAGLE_GT,
        idx=BEAGLE_GT + ".csi",
    params:
        prefix=BEAGLE_GT.replace(".vcf.gz", ""),
        old_vcf=temp(BEAGLE_GT.replace(".vcf.gz", "old.vcf.gz")),
        unzipped_vcf=temp(BEAGLE_GT.replace(".vcf.gz", ".vcf")),
    threads:
        16
    shell:
        "java -jar beagle.12Jul19.0df.jar \
        gt={input.vcf} \
        ref={input.panel} \
        map={input.gmap} \
        out={params.prefix};\n"
        "bcftools index {output.vcf};\n"

rule extract_lengths:
    input:
        vcf=REF_PANEL
    output:
        REF_LENGTHS
    shell:
        "bcftools view -h {input} | grep '##contig' | sed -E 's/##contig=<ID=(.*),assembly=.*,length=([0-9]+)>/\\1\t\\2/' > {output}"
        
rule index_personal:
    input:
        ref=config["index"],
        vcf=BEAGLE_GT,
        idx=BEAGLE_GT + ".csi"
    output:
        fa=PG,
        idx1=PG+".1.bt2",
        idx2=PG+".2.bt2",
        idx3=PG+".3.bt2",
        idx4=PG+".4.bt2",
        idx5=PG+".rev.1.bt2",
        idx6=PG+".rev.2.bt2"
    threads: 16
    params:
        sample="{sample}",
        i="{i}"
    shell:
        "bcftools consensus -f {input.ref} -H {params.i} -s {params.sample} {input.vcf} > {output.fa};\n"
        "bowtie2-build --threads {threads} {output.fa} {output.fa}"

rule serialize_liftover:
    input:
       vcf=BEAGLE_GT,
       lengths=REF_LENGTHS
    output:
       LIFT
    params:
        sample="{sample}",
        i=lambda w: int(w.i) - 1,
        prefix=LIFT.replace(".lft", ""),
    shell:
        "liftover/liftover serialize -v {input.vcf} --haplotype {params.i} -s {params.sample} -p {params.prefix} -k {input.lengths}"
