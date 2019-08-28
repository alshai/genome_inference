import math 
import os

configfile: "config.yaml"

READS       =   "{sample}/reads.fq.gz"
VCF         =   config["vcf"].replace(".vcf.gz", ".filtered.vcf.gz")
REF_PANEL   =   config["vcf"].replace(".vcf.gz", ".filtered.leave_out.vcf.gz")
DIR         =   "{sample}/{coverage}x"
snakemake.utils.makedirs(expand(DIR, 
                                sample=config["samples"], 
                                coverage=config["coverage"]))
GENOME      =   os.path.join("{sample}", "genome.fa")
BAM         =   os.path.join(DIR, "hg19_alns.bam") 
COUNT_GT    =   os.path.join(DIR, "genotype_by_count.vcf.gz")
IMPUTE_GT   =   os.path.join(DIR, "genotype_by_impute.vcf.gz")
PG          =   os.path.join(DIR, "pg.fa")
LIFT        =   os.path.join(DIR, "pg.lft")

print(expand(IMPUTE_GT, sample=config["samples"], coverage=config["coverage"]))
rule all:
    input:
        expand(IMPUTE_GT, sample=config["samples"], coverage=config["coverage"]),
        expand(PG, sample=config["samples"], coverage=config["coverage"]),
        expand(LIFT, sample=config["samples"], coverage=config["coverage"])


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

rule simulate_haplotype:
    input:
        ref=config["index"],
        vcf=VCF,
        vcf_idx=VCF + ".csi"
    output:
        fa=GENOME
    params:
        sample = "{sample}"
    shell:
        "bcftools consensus -f {input.ref} -H 1 -s {params.sample} {input.vcf} > {output.fa}"

rule simulate_reads:
    input:
        genome=GENOME
    output:
        fq=READS,
        sam=READS.replace(".fq.gz", ".sam")
    params:
        nreads=config["nreads"],
        rlen=config["read_length"]
    threads: 16
    shell:
        "mason_simulator --num-threads {threads} "
                         "-ir {input.genome} "
                         "-n {params.nreads} "
                         "--illumina-read-length {params.rlen} "
                         "-o {output.fq} "
                         "-oa {output.sam}"

rule align_subset:
    input:
        fq=READS,
        idx1=config["index"] + ".1.bt2"
    params:
        fracs=lambda wildcards: float(wildcards.coverage) / config["total_coverage"],
        index=config["index"]
    output:
        bam=BAM
    threads: 16
    shell:
        "bowtie2 --threads {threads} -x {params.index} -U <(seqtk sample {input.fq} {params.fracs}) | samtools sort -@ {threads} > {output.bam}"

rule genotype:
    input:
        bam=BAM,
        vcf=VCF
    output:
        vcf=COUNT_GT
    params:
        thres=lambda wildcards: math.ceil(float(wildcards.coverage) / 2),
        sample="{sample}"
    shell:
        "bin/varcount -s {params.sample} -dgc {params.thres} {input.vcf} {input.bam} | "
        "bcftools sort -O z > {output.vcf};\n"
        "bcftools index {output.vcf}"

# rule phase:
#     input:
#         vcf=COUNT_GT,
#         panel=IN_VCF,
#         gmap=config['shapeit_map']
#     output:
#         vcf=PHASE_GT,
#     params:
#         sample="{sample}"
#     threads:
#         16
#     benchmark:
#         "benchmarks/{sample}/{coverage}/impute.benchmark"
#     shell:
#         "bin/shapeit4 --input {input.vcf} \
#         --map {input.gmap} \
#         --region 21 \
#         --output {output.vcf} \
#         --thread {threads} \
#         --reference {input.panel};\n"
# 
#         "bcftools index {output.vcf};\n"

rule impute:
    input:
        vcf=COUNT_GT,
        panel=REF_PANEL,
        gmap=config['beagle_map']
    output:
        vcf=IMPUTE_GT,
        idx=IMPUTE_GT + ".csi",
    params:
        prefix=IMPUTE_GT.replace(".vcf.gz", ""),
        old_vcf=temp(IMPUTE_GT.replace(".vcf.gz", "old.vcf.gz")),
        unzipped_vcf=temp(IMPUTE_GT.replace(".vcf.gz", ".vcf")),
    threads:
        16
    shell:
        "java -jar beagle.12Jul19.0df.jar \
        gt={input.vcf} \
        ref={input.panel} \
        map={input.gmap} \
        out={params.prefix};\n"
        # fix the header
        "mv {output.vcf} {params.old_vcf};\n"
        "bcftools view -h {params.old_vcf} | sed '/##contig/d'  | head -n -1 >> {params.unzipped_vcf};\n"
        "bcftools view -h {input.vcf} | grep -P '##contig|#CHROM' >> {params.unzipped_vcf};\n"
        "bcftools view -H {params.old_vcf} >> {params.unzipped_vcf};\n"
        "bgzip {params.unzipped_vcf};\n"
        "bcftools index {output.vcf};\n"
        "rm {params.old_vcf};\n"

rule index_personal:
    input:
        ref=config["index"],
        vcf=IMPUTE_GT,
        idx=IMPUTE_GT + ".csi"
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
        sample="{sample}"
    shell:
        "bcftools consensus -f {input.ref} -H 1 -s {params.sample} {input.vcf} > {output.fa};\n"
        "bowtie2-build --threads {threads} {output.fa} {output.fa}"

rule serialize_liftover:
    input:
       vcf=IMPUTE_GT
    output:
       LIFT
    params:
        sample="{sample}",
        prefix=LIFT.replace(".lft", "")
    shell:
        "bin/liftover serialize -v {input.vcf} -s {params.sample} -p {params.prefix}"
