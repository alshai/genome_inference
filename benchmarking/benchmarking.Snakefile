import math
import os
import random
from command_dict import cmds
from scripts.summarize import *
from scripts.dist_windows import dist_windows

configfile: "config.yaml"

# TODO: figure out how to allow anything except for hg19 and personal
wildcard_constraints:
    sample="NA.+|HG.+",
    other_sample="NA.+|HG.+",

IMP_PREFIX    =  "{sample}/{coverage}x/{genotyper}/{filter}/{fname}"
PERS_MERGED_SCORE   =  "{sample}/personal/merged.score"
HG19_SCORE          =  "{sample}/GRCh37/alns.score"
MA_SCORE            = "{sample}/ma/lifted_md.score"
TEST_SAMPLES="data/test_samples.txt"
REF_SAMPLES="data/ref_samples.txt"

if os.path.isfile(TEST_SAMPLES) and os.path.isfile(REF_SAMPLES):
    SAMPLES = open(TEST_SAMPLES).read().strip().split("\n")
else:
    SAMPLES = []

rule all:
    input:
        # TEST_SAMPLES, REF_SAMPLES,
        # expand("{sample}/{x}/summary.txt",
        #        sample=SAMPLES,
        #        x=SAMPLES + ["personal", "GRCh37", "ma"]),
        # expand("{sample}/{x}/dists.txt",
        #        sample=SAMPLES,
        #        x=SAMPLES + ["personal", "GRCh37", "ma"]),
        expand(IMP_PREFIX,
               sample=["NA20534"],
               coverage=config["coverage"],
               genotyper=["likelihood_naive"],
               filter=["aa+ar"],
               # filter=["ar", "aa+ar+some_rr", "aa+ar"],
               fname=["merged.score", "summary.txt", "dists.txt"])




###### PREPARE DATA FOR LATER STAGES ############

VCF         =   config["vcf"].replace(".vcf.gz", ".filtered.vcf.gz")
ONEKG_GT    =   "{sample}/genotype.vcf.gz"
REF_PANEL   =   config["vcf"].replace(".vcf.gz", ".filtered.leave_out.vcf.gz")
TRAIN_READS =   "{sample}/train_reads.fq.gz"
TRAIN_TRUTH =   TRAIN_READS.replace(".fq.gz", ".lifted.bam")
TEST_READS  =   "{sample}/test_reads.fq.gz"
TEST_TRUTH  =   TEST_READS.replace(".fq.gz", ".lifted.bam")
REF_LENGTHS =   "data/ref_lengths.txt"
GENOME      =   "{sample}/genome{i}.fa"

rule filter_vcf:
    input:
        vcf=config["vcf"]
    output:
        vcf=VCF,
        idx=VCF + ".csi"
    shell:
        "bcftools view -O z -V mnps,other {input.vcf} > {output.vcf};\n"
        "bcftools index {output.vcf}"

rule choose_samples:
    input:
        vcf=VCF,
    output:
        t=TEST_SAMPLES,
        r=REF_SAMPLES
    params:
        ntest=config["n_test_samples"]
    run:
        samples = set()
        for line in shell("bcftools query -l {input.vcf}", iterable=True):
            samples.add(line.strip())
        test_samples = set(random.sample(samples, params.ntest))
        ref_samples = samples - test_samples
        with open(output.t, "w") as fp:
            fp.write("\n".join(test_samples))
        with open(output.r, "w") as fp:
            fp.write("\n".join(ref_samples))

rule make_ref_panel:
    input:
        vcf=VCF,
        idx=VCF + ".csi",
        samples=REF_SAMPLES
    output:
        vcf=REF_PANEL
    shell:
        "bcftools view -Oz -S {input.samples} {input.vcf} > {output.vcf};\n"
        "bcftools index {output.vcf}\n"


rule index_hg19:
    input:
        config["hg19_index"]
    output:
        config["hg19_index"] + ".1.bt2",
        config["hg19_index"] + ".2.bt2",
        config["hg19_index"] + ".3.bt2",
        config["hg19_index"] + ".4.bt2",
        config["hg19_index"] + ".rev.1.bt2",
        config["hg19_index"] + ".rev.2.bt2"
    threads: 16
    shell:
        "bowtie2-build --threads {threads} {input} {input}"

rule save_genotype:
    input:
        vcf=VCF
    params:
        sample="{sample}"
    output:
        vcf=ONEKG_GT,
        vcf_idx=ONEKG_GT + '.csi'
    shell:
        "bcftools view -Ou -s {params.sample} {input} | "
        "bcftools view -i 'GT~\"A\"' -Ou | "
        "bcftools norm -Oz -d indels > {output.vcf};\n"
        "bcftools index --force {output.vcf}"

rule simulate_haplotype:
    input:
        ref=config["hg19_index"],
        vcf=ONEKG_GT,
        vcf_idx=ONEKG_GT + ".csi"
    output:
        fa=GENOME
    params:
        sample="{sample}",
        i="{i}"
    shell:
        "bcftools consensus -f {input.ref} -H {params.i} -s {params.sample} {input.vcf} > {output.fa} 2>/dev/null;"

rule simulate_reads:
    input:
        "{sample}/genome1.fa",
        "{sample}/genome2.fa"
    output:
        fa=temp("{sample}/full_genome.fa"),
        fq1=TRAIN_READS,
        sam1=TRAIN_READS.replace(".fq.gz", ".sam"),
        fq2=TEST_READS,
        sam2=TEST_READS.replace(".fq.gz", ".sam")
    params:
        ntrain=config["n_train_reads"],
        ntest=config["n_test_reads"],
        rlen=config["read_length"]
    threads: 16
    shell:
        "cat {input} | awk 'BEGIN {{i=1}} /^>/ {{print \">genome\"i++}} /^[^>]/ {{print}}' > {output.fa};\n"
        "mason_simulator --num-threads {threads} "
                         "-ir {output.fa} "
                         "-n {params.ntrain} "
                         "--illumina-read-length {params.rlen} "
                         "-o {output.fq1} "
                         "-oa {output.sam1};\n"
        "mason_simulator --num-threads {threads} "
                         "-ir {output.fa} "
                         "-n {params.ntest} "
                         "--illumina-read-length {params.rlen} "
                         "-o {output.fq2} "
                         "-oa {output.sam2}"

rule split_test_read_alignments:
    input:
        sam=TEST_READS.replace(".fq.gz", ".sam")
    output:
        sam1=TEST_READS.replace(".fq.gz", ".1.sam"),
        sam2=TEST_READS.replace(".fq.gz", ".2.sam")
    script:
        "scripts/split_sam_by_ref.py"

rule lift_test_read_alignments:
    input:
        vcf=ONEKG_GT,
        sam=TEST_READS.replace(".fq.gz", ".{i}.sam"),
    output:
        sam=TEST_READS.replace(".fq.gz", ".{i}.lifted.sam")
    params:
        sample = "{sample}",
        h = lambda w: int(w.i) - 1,
        i = "{i}"
    shell:
        "liftover/liftover lift -v {input.vcf} -a {input.sam} -s {params.sample} --haplotype {params.h} -n <(echo '21	genome{params.i}') > {output}"

rule test_concat:
    input:
        one=TEST_READS.replace(".fq.gz", ".1.lifted.sam"),
        two=TEST_READS.replace(".fq.gz", ".2.lifted.sam")
    output:
        TEST_TRUTH
    shell:
        "samtools merge -u {output} {input.one} {input.two}"

rule extract_lengths:
    input:
        vcf=VCF
    output:
        REF_LENGTHS
    shell:
        "bcftools view -h {input} | grep '##contig' | sed -E 's/##contig=<ID=(.*),assembly=.*,length=([0-9]+)>/\\1\t\\2/' > {output}"

##### IMPUTATION START #######
genotyping_strategy = {
    'likelihood_naive': 'likelihood'
    # TODO: add majority-count and alt_sensitive here
}

# params:
#     cmd = lambda w: filtering_strategy[w.filter]
# run:
#     cmd = filtering_strategy[w.filter]
#     shell(cmd)
filtering_strategy = {
    'aa'                : "bcftools view -i 'GT~\"AA\"' {input.vcf} > {output.tmp};\n",
    'ar'                : "bcftools view -i 'GT~\"AR\"' {input.vcf} > {output.tmp};\n",
    'aa+ar'             : "bcftools view -i 'GT~\"A\"' {input.vcf} > {output.tmp};\n",
    'aa+ar+some_rr'     : "bcftools view -i 'GT~\"A\"' {input.vcf} > {output.tmp};\nbcftools view -H -i 'GT~\"RR\"' {input.vcf} | shuf -n 50000 >> {output.tmp};\n"
}

IMP_DIR               =   "{sample}/{coverage}x"
IMP_SUBSET_SAM        =   "{sample}/{coverage}x/subset_v_GRCh37.sam"
IMP_COUNTS            =   "{sample}/{coverage}x/counts.{genotyper}.vcf.gz"
IMP_FILTERED_COUNTS   =   "{sample}/{coverage}x/{genotyper}/{filter}/counts.vcf.gz"
IMP_BEAGLE_GT         =   "{sample}/{coverage}x/{genotyper}/{filter}/imputed.vcf.gz"
IMP_FILTERED_GT         =   "{sample}/{coverage}x/{genotyper}/{filter}/imputed.filtered.vcf.gz"
IMP_PG                =   "{sample}/{coverage}x/{genotyper}/{filter}/pg.{i}.fa"
IMP_LIFT              =   "{sample}/{coverage}x/{genotyper}/{filter}/pg.{i}.lft"

rule imp_align_subset:
    input:
        reads=TRAIN_READS,
        fa=config["hg19_index"],
        idx1=config["hg19_index"] + ".1.bt2"
    params:
        frac=lambda w: float(w.coverage) / config["total_coverage"], # multiply by 2 for diploid?
    output:
        sam=IMP_SUBSET_SAM
    threads: 16
    shell:
        "bowtie2 -p {threads} -x {input.fa} -U <(seqtk sample {input.reads} {params.frac}) > {output.sam}"

# genotyper is included here!!
rule imp_count:
    input:
        bam=IMP_SUBSET_SAM,
        vcf=VCF
    output:
        vcf=IMP_COUNTS
    params:
        sample="{sample}",
        gt=lambda w: genotyping_strategy[w.genotyper]
    shell:
        "varcount/varcount -g{params.gt} -s {params.sample} {input.vcf} {input.bam} | "
        "bcftools sort -O z > {output.vcf};\n"
        "bcftools index {output.vcf}"

rule imp_filter_counts:
    input:
        vcf=IMP_COUNTS
    output:
        vcf=IMP_FILTERED_COUNTS,
        tmp=temp(IMP_FILTERED_COUNTS + "_tmp.vcf")
    params:
        filter_method=lambda w: filtering_strategy[w.filter]
    run:
        cmd = params.filter_method
        shell(cmd)
        shell("bcftools sort -Oz {output.tmp} > {output.vcf};\n")
        shell("bcftools index {output.vcf};\n")

rule imp_beagle_impute:
    input:
        vcf=IMP_FILTERED_COUNTS,
        panel=REF_PANEL,
        gmap=config['beagle_map']
    output:
        vcf=IMP_BEAGLE_GT,
        idx=IMP_BEAGLE_GT + ".csi",
    params:
        prefix=IMP_BEAGLE_GT.replace(".vcf.gz", ""),
    threads:
        16
    shell:
        "java -jar beagle.12Jul19.0df.jar \
        nthreads={threads} \
        gt={input.vcf} \
        ref={input.panel} \
        map={input.gmap} \
        out={params.prefix};\n"
        "bcftools index {output.vcf};\n"

# outputs a filterd, sorted vcf file
rule imp_beagle_filter:
    input:
        vcf=IMP_BEAGLE_GT,
        vcf_idx=IMP_BEAGLE_GT + ".csi"
    output:
        vcf=IMP_FILTERED_GT,
        vcf_idx=IMP_FILTERED_GT + ".csi"
    threads: 16
    shell:
        "bcftools norm -Ou -d indels {input.vcf} | "
        "bcftools sort -Oz > {output.vcf};\n"
        "bcftools index --force {output.vcf};\n"

rule imp_index:
    input:
        ref=config["hg19_index"],
        vcf=IMP_FILTERED_GT,
        idx=IMP_FILTERED_GT + ".csi"
    output:
        fa=IMP_PG,
        idx1=IMP_PG+".1.bt2",
        idx2=IMP_PG+".2.bt2",
        idx3=IMP_PG+".3.bt2",
        idx4=IMP_PG+".4.bt2",
        idx5=IMP_PG+".rev.1.bt2",
        idx6=IMP_PG+".rev.2.bt2"
    threads: 16
    params:
        sample="{sample}",
        i="{i}"
    shell:
        "bcftools consensus -f {input.ref} -H {params.i} -s {params.sample} {input.vcf} > {output.fa} 2> /dev/null;\n"
        "bowtie2-build --threads {threads} {output.fa} {output.fa}"

# rule imp_liftover_serialize:
#     input:
#        vcf=IMP_BEAGLE_GT,
#        lengths=REF_LENGTHS
#     output:
#        IMP_LIFT
#     params:
#         sample="{sample}",
#         i=lambda w: int(w.i) - 1,
#         prefix=IMP_LIFT.replace(".lft", ""),
#     shell:
#         "liftover/liftover serialize -v {input.vcf} --haplotype {params.i} -s {params.sample} -p {params.prefix} -k {input.lengths}"

##### IMPUTATION ALIGNMENTS #####
IMP_UNLIFTED_ALNS     =   os.path.join(IMP_DIR, "{genotyper}/{filter}/unlifted.{i}.sam")
IMP_LIFTED_ALNS       =   os.path.join(IMP_DIR, "{genotyper}/{filter}/lifted.{i}.sam")
IMP_LIFTED_SCORE       =   os.path.join(IMP_DIR, "{genotyper}/{filter}/lifted.{i}.score")
IMP_MERGED_ALNS       =   os.path.join(IMP_DIR, "{genotyper}/{filter}/merged.sam")
IMP_MERGED_SCORE      =   os.path.join(IMP_DIR, "{genotyper}/{filter}/merged.score")

rule imp_align_all:
    input:
        fa=IMP_PG,
        idx1=IMP_PG+".1.bt2",
        idx2=IMP_PG+".2.bt2",
        idx3=IMP_PG+".3.bt2",
        idx4=IMP_PG+".4.bt2",
        idx5=IMP_PG+".rev.1.bt2",
        idx6=IMP_PG+".rev.2.bt2",
        reads=TEST_READS
    output:
        sam=IMP_UNLIFTED_ALNS
    threads: 16
    shell:
        "bowtie2 -p {threads} -x {input.fa} -U {input.reads} > {output.sam}"

rule imp_lift_alns:
    input:
        sam=IMP_UNLIFTED_ALNS,
        vcf=IMP_FILTERED_GT,
        lengths=REF_LENGTHS
    output:
        sam=IMP_LIFTED_ALNS
    params:
        i=lambda w: int(w.i) - 1,
        sample="{sample}",
        prefix=IMP_LIFTED_ALNS.replace(".sam","")
    shell:
        cmds["lift"]

rule imp_score_lifted:
    input:
        truth=TEST_TRUTH,
        sam=IMP_LIFTED_ALNS
    output:
        sam=IMP_LIFTED_SCORE
    shell:
        cmds["score"]

rule imp_merge_alns:
    input:
        a=IMP_LIFTED_ALNS.replace("{i}", "1"),
        b=IMP_LIFTED_ALNS.replace("{i}", "2"),
        fa=config["hg19_index"]
    output:
        IMP_MERGED_ALNS
    shell:
        "hts_utils/merge_sams {input.a} {input.b} | "
        "samtools calmd - {input.fa} > {output}"

rule imp_score_merged:
    input:
        truth=TEST_TRUTH,
        sam=IMP_MERGED_ALNS
    output:
        IMP_MERGED_SCORE
    shell:
        cmds["score"]

#### REFERENCE ALIGNMENTS ######
HG19_VCF               = "data/hg19.chr21.vcf.gz"
HG19_ALNS              =   "{sample}/GRCh37/alns.sam"
HG19_SCORE      =   "{sample}/GRCh37/alns.score"

rule HG19_vcf:
    input:
        VCF
    output:
        vcf=HG19_VCF,
        temp_vcf=temp(HG19_VCF.replace(".gz","")),
        idx=HG19_VCF.replace(".gz", ".gz.csi")
    run:
        ofp = open(output.temp_vcf, "w")
        for line in shell("bcftools view -Ou --force-samples -s' ' data/chr21.filtered.vcf.gz | bcftools annotate -h <(echo '##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')" , iterable=True):
            if "#CHROM" in line:
                line = line.rstrip() + "\tFORMAT\tGRCh37"
            elif line[0] != "#":
                line = line.rstrip() + "\tGT\t0|0"
            ofp.write(line + '\n')
        ofp.close()
        shell("bcftools view -Oz {output.temp_vcf} > {output.vcf}")
        shell("bcftools index {output.vcf}")

rule hg19_align_all:
    input:
        fa=config["hg19_index"],
        idx1=config["hg19_index"]+".1.bt2",
        idx2=config["hg19_index"]+".2.bt2",
        idx3=config["hg19_index"]+".3.bt2",
        idx4=config["hg19_index"]+".4.bt2",
        idx5=config["hg19_index"]+".rev.1.bt2",
        idx6=config["hg19_index"]+".rev.2.bt2",
        reads=TEST_READS
    output:
        sam=HG19_ALNS
    threads: 16
    shell:
        cmds["align"]

rule hg19_score_alns:
    input:
        truth=TEST_TRUTH,
        sam=HG19_ALNS
    output:
        HG19_SCORE
    shell:
        cmds["score"]

#### MAJOR ALLELE ALIGNMENTS ######
MA_VCF = "data/major_allele.chr21.vcf.gz"
MA_GENOME = "data/major_allele.chr21.fa"
MA_UNLIFTED_ALNS    =   "{sample}/ma/unlifted.sam"
MA_LIFTED_ALNS      =   "{sample}/ma/lifted.sam"
MA_ALNS             =   "{sample}/ma/lifted_md.sam"
MA_SCORE            =   "{sample}/ma/lifted_md.score"

rule ma_vcf:
    input:
        VCF
    output:
        vcf=MA_VCF,
        temp_vcf=temp(MA_VCF.replace(".gz","")),
        idx=MA_VCF.replace(".gz", ".gz.csi")
    run:
        ofp = open(output.temp_vcf, "w")
        for line in shell("bcftools view -Ou --force-samples -s' ' -i 'INFO/AF>0.5' data/chr21.filtered.vcf.gz | bcftools annotate -h <(echo '##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')" , iterable=True):
            if "#CHROM" in line:
                line = line.rstrip() + "\tFORMAT\tMA"
            elif line[0] != "#":
                line = line.rstrip() + "\tGT\t1|1"
            ofp.write(line + '\n')
        ofp.close()
        shell("bcftools view -Oz {output.temp_vcf} > {output.vcf}")
        shell("bcftools index {output.vcf}")

rule ma_fasta_index:
    input:
        gt=MA_VCF,
        idx=MA_VCF.replace(".vcf.gz", ".vcf.gz.csi"),
        fa=config["hg19_index"]
    output:
        fa=MA_GENOME,
        idx1=MA_GENOME+".1.bt2",
        idx2=MA_GENOME+".2.bt2",
        idx3=MA_GENOME+".3.bt2",
        idx4=MA_GENOME+".4.bt2",
        idx5=MA_GENOME+".rev.1.bt2",
        idx6=MA_GENOME+".rev.2.bt2"
    params:
    run:
        shell("bcftools consensus -f {input.fa} -H 1 {input.gt} > {output.fa}")
        shell("bowtie2-build {output.fa} {output.fa}")

rule ma_align_all:
    input:
        fa=MA_GENOME,
        idx1=MA_GENOME+".1.bt2",
        idx2=MA_GENOME+".2.bt2",
        idx3=MA_GENOME+".3.bt2",
        idx4=MA_GENOME+".4.bt2",
        idx5=MA_GENOME+".rev.1.bt2",
        idx6=MA_GENOME+".rev.2.bt2",
        reads=TEST_READS
    output:
        sam=MA_UNLIFTED_ALNS
    threads: 16
    shell:
        cmds["align"]

rule ma_lift_alns:
    input:
        vcf=MA_VCF,
        sam=MA_UNLIFTED_ALNS,
        lengths=REF_LENGTHS
    output:
        sam=MA_LIFTED_ALNS
    params:
        sample="MA",
        prefix=MA_LIFTED_ALNS.replace(".sam", ""),
        i=0
    shell:
        cmds["lift"]

rule ma_md_lifted:
    input:
        sam=MA_LIFTED_ALNS,
        fa=config["hg19_index"]
    output:
        MA_ALNS
    shell:
        "samtools calmd {input.sam} {input.fa} > {output}"

rule ma_score_alns:
    input:
        truth=TEST_TRUTH,
        sam=MA_ALNS
    output:
        MA_SCORE
    shell:
        cmds["score"]


###### PERSONAL GENOME #######
PERS_UNLIFTED_ALNS     =   "{sample}/personal/unlifted.{i}.sam"
PERS_LIFTED_ALNS       =   "{sample}/personal/lifted.{i}.sam"
PERS_LIFTED_SCORES       =   "{sample}/personal/lifted.{i}.score"
PERS_MERGED_ALNS       =   "{sample}/personal/merged.sam"
PERS_MERGED_SCORE      =   "{sample}/personal/merged.score"

rule personal_index:
    input:
        fa=GENOME
    output:
        idx1=GENOME+".1.bt2",
        idx2=GENOME+".2.bt2",
        idx3=GENOME+".3.bt2",
        idx4=GENOME+".4.bt2",
        idx5=GENOME+".rev.1.bt2",
        idx6=GENOME+".rev.2.bt2"
    threads: 16
    shell:
        "bowtie2-build --threads {threads} {input} {input}"

rule personal_align_all:
    input:
        fa=GENOME,
        reads=TEST_READS,
        idx1=GENOME+".1.bt2",
    output:
        sam=PERS_UNLIFTED_ALNS
    threads: 16
    shell:
        cmds["align"]

rule personal_lift_alns:
    input:
        vcf=ONEKG_GT,
        sam=PERS_UNLIFTED_ALNS,
        lengths=REF_LENGTHS
    output:
        sam=PERS_LIFTED_ALNS
    params:
        sample="{sample}",
        i=lambda w: int(w.i) - 1,
        prefix=PERS_LIFTED_ALNS.replace(".sam", "")
    shell:
        cmds["lift"]

rule personal_score_lifted:
    input:
        truth=TEST_TRUTH,
        sam=PERS_LIFTED_ALNS
    threads: 16
    output:
        PERS_LIFTED_SCORES
    shell:
        cmds["score"]


rule personal_merge_alns:
    input:
        a=PERS_LIFTED_ALNS.replace("{i}", "1"),
        b=PERS_LIFTED_ALNS.replace("{i}","2"),
        fa=config["hg19_index"]
    output:
        PERS_MERGED_ALNS
    shell:
        "hts_utils/merge_sams {input.a} {input.b} | "
        "samtools calmd - {input.fa} > {output}"


rule personal_score_merged:
    input:
        truth=TEST_TRUTH,
        sam=PERS_MERGED_ALNS
    threads: 16
    output:
        PERS_MERGED_SCORE
    shell:
        cmds["score"]


### OTHER GENOME ALIGNMENTS ####
OTHER_GENOME            =   "{other_sample}/genome{i}.fa"
OTHER_VCF               =   "{other_sample}/genotype.vcf.gz"
OTHER_UNLIFTED_ALNS     =   "{sample}/{other_sample}/unlifted.{i}.sam"
OTHER_LIFTED_ALNS       =   "{sample}/{other_sample}/lifted.{i}.sam"
OTHER_MERGED_ALNS       =   "{sample}/{other_sample}/merged.sam"
OTHER_MERGED_SCORE      =   "{sample}/{other_sample}/merged.score"
OTHER_SCORE             =   "{sample}/{other_sample}/alns.score"

rule other_align_all:
    input:
        fa=OTHER_GENOME,
        reads=TEST_READS,
        idx1=OTHER_GENOME+".1.bt2",
    output:
        sam=OTHER_UNLIFTED_ALNS
    threads: 16
    shell:
        cmds["align"]

rule other_lift_alns:
    input:
        vcf=OTHER_VCF,
        sam=OTHER_UNLIFTED_ALNS,
        lengths=REF_LENGTHS
    output:
        sam=OTHER_LIFTED_ALNS
    params:
        sample="{other_sample}",
        i=lambda w: int(w.i) - 1,
        prefix=OTHER_LIFTED_ALNS.replace(".sam", "")
    shell:
        cmds["lift"]

rule other_merge_alns:
    input:
        a=OTHER_LIFTED_ALNS.replace("{i}", "1"),
        b=OTHER_LIFTED_ALNS.replace("{i}","2"),
        fa=config["hg19_index"]
    output:
        OTHER_MERGED_ALNS
    shell:
        "hts_utils/merge_sams {input.a} {input.b} | "
        "samtools calmd - {input.fa} > {output}"

rule other_score_alns:
    input:
        truth=TEST_TRUTH,
        sam=OTHER_MERGED_ALNS
    threads: 16
    output:
        OTHER_MERGED_SCORE
    shell:
        cmds["score"]


### SUMMARIZE EXPERIMENTS ###
rule imp_summarize_exp:
    input:
        vcf1=ONEKG_GT,
        vcf2=IMP_FILTERED_GT,
        scores=IMP_MERGED_SCORE,
        alns=IMP_MERGED_ALNS,
    output:
        out="{sample}/{coverage}x/{genotyper}/{filter}/summary.txt",
        out2="{sample}/{coverage}x/{genotyper}/{filter}/dists.txt",
    params:
        s1="{sample}",
        s2="{sample}",
        experiment="IMP",
        ploidy="diploid",
        coverage="{coverage}",
        imp_input="{filter}",
    run:
        summarize(input, output, params, shell)
        dist_windows(input, output, params, shell)

rule personal_summarize_exp:
    input:
        vcf1=ONEKG_GT,
        vcf2=ONEKG_GT,
        scores=PERS_MERGED_SCORE,
        alns=PERS_MERGED_ALNS,
    output:
        out="{sample}/personal/summary.txt",
        out2="{sample}/personal/dists.txt",
    params:
        s1="{sample}",
        s2="{sample}",
        experiment="PERSONAL",
        ploidy="diploid",
        coverage=0,
        imp_input="PERSONAL",
    run:
        summarize(input, output, params, shell)
        dist_windows(input, output, params, shell)

rule hg19_summarize_exp:
    input:
        vcf1=ONEKG_GT,
        vcf2=HG19_VCF,
        alns=HG19_ALNS,
        scores=HG19_SCORE,
    output:
        out="{sample}/GRCh37/summary.txt",
        out2="{sample}/GRCh37/dists.txt",
    params:
        s1="{sample}",
        s2="GRCh37",
        experiment="GRCh37",
        ploidy="haploid",
        coverage="0",
        imp_input="GRCh37",
    run:
        summarize(input, output, params, shell)
        dist_windows(input, output, params, shell)

rule other_summarize_exp:
    input:
        vcf1=ONEKG_GT,
        vcf2=OTHER_VCF,
        scores=OTHER_MERGED_SCORE,
        alns=OTHER_MERGED_ALNS,
    output:
        out="{sample}/{other_sample}/summary.txt",
        out2="{sample}/{other_sample}/dists.txt"
    params:
        s1="{sample}",
        s2="{other_sample}",
        experiment="OTHER",
        ploidy="diploid",
        coverage="0",
        imp_input="{other_sample}",
    run:
        summarize(input, output, params, shell)
        dist_windows(input, output, params, shell)

rule ma_summarize_exp:
    input:
        vcf1=ONEKG_GT,
        vcf2=MA_VCF,
        scores=MA_SCORE,
        alns=MA_ALNS,
    output:
        out="{sample}/ma/summary.txt",
        out2="{sample}/ma/dists.txt",
    params:
        s1="{sample}",
        s2="MA",
        experiment="MA",
        ploidy="haploid",
        coverage="0",
        imp_input="MA",
    run:
        summarize(input, output, params, shell)
        dist_windows(input, output, params, shell)

#### ANALYZE RESCUED READS #####
IMP_RESCUED_READS      =   os.path.join(IMP_DIR, "{genotyper}/{filter}/rescued_reads.bed")
IMP_RESCUED_REPEATS      =   os.path.join(IMP_DIR, "{genotyper}/{filter}/rescued_repeats.bed")
rule rescued_reads:
    input:
        imp_score=IMP_MERGED_SCORE,
        ref_score=HG19_SCORE
    output:
        IMP_RESCUED_READS
    script:
        "scripts/get_rescued_reads.py"

rule rescued_repeats:
    input:
        a=IMP_RESCUED_READS,
        b=config["repeat_bed"]
    output:
        IMP_RESCUED_REPEATS
    shell:
        "bedtools intersect -wb -a {input.a} -b {input.b} > {output}"
