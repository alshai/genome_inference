import math 
import os
import random
from command_dict import cmds

configfile: "config.yaml"

# TODO: figure out how to allow anything except for hg19 and personal
wildcard_constraints:
    sample="NA.+|HG.+",
    other_sample="NA.+|HG.+",

IMP_MERGED_SCORE    =  "{sample}/{coverage}x/{genotyper}/{filter}/merged.score"
PERS_MERGED_SCORE   =  "{sample}/personal/merged.score"
HG19_SCORE          =  "{sample}/GRCh37/alns.score"
TEST_SAMPLES="data/test_samples.txt"
REF_SAMPLES="data/ref_samples.txt"

if os.path.isfile(TEST_SAMPLES) and os.path.isfile(REF_SAMPLES):
    SAMPLES = open(TEST_SAMPLES).read().strip().split("\n")
else:
    SAMPLES = []

print(TEST_SAMPLES)
rule all:
    input:
        TEST_SAMPLES, REF_SAMPLES,
        expand("{sample}/GRCh37/summary.txt", sample=SAMPLES),
        expand("{sample}/personal/summary.txt", sample=SAMPLES),
        expand("{sample}/{coverage}x/{genotyper}/{filter}/summary.txt",
               genotyper=["likelihood_naive"],
               filter=['aa+ar+some_rr', 'ar', 'aa+ar'],
               sample=SAMPLES, 
               coverage=config["coverage"]
              ),
        expand("{sample}/{other_sample}/summary.txt", sample=SAMPLES, other_sample=SAMPLES),
        expand("{sample}/ma/summary.txt", sample=SAMPLES, other_sample=SAMPLES),
        # expand("{sample}/{coverage}x/{genotyper}/{filter}/rescued_repeats.bed",
        #        genotyper=["likelihood_naive"],
        #        filter=['ar', 'aa+ar', 'aa+ar+some_rr'],
        #        sample=["HG02374"], 
        #        coverage=[10]
        #        ),
 
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
        config["hg19_index"] + ".amb",
        config["hg19_index"] + ".ann",
        config["hg19_index"] + ".bwt",
        config["hg19_index"] + ".pac",
        config["hg19_index"] + ".sa",
    threads: 16
    shell:
        # "bowtie2-build --threads {threads} {input} {input}"
        "bwa index {input} {input}"

rule save_genotype:
    input:
        vcf=VCF
    params:
        sample="{sample}"
    output:
        ONEKG_GT
    shell:
        "bcftools view -Ou -s {params.sample} {input} | bcftools view -i 'GT~\"A\"' -Oz > {output}"

rule simulate_haplotype:
    input:
        ref=config["hg19_index"],
        vcf=VCF,
        vcf_idx=VCF + ".csi"
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
        vcf=VCF,
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
IMP_PG                =   "{sample}/{coverage}x/{genotyper}/{filter}/pg.{i}.fa"
IMP_LIFT              =   "{sample}/{coverage}x/{genotyper}/{filter}/pg.{i}.lft"

rule imp_align_subset:
    input:
        reads=TRAIN_READS,
        fa=config["hg19_index"],
        idx1=config["hg19_index"] + ".amb"
    params:
        frac=lambda w: float(w.coverage) / config["total_coverage"], # multiply by 2 for diploid?
    output:
        sam=IMP_SUBSET_SAM
    threads: 16
    shell:
        # "bowtie2 -p {threads} -x {input.fa} -U <(seqtk sample {input.reads} {params.frac}) > {output.sam}"
        "bwa mem -t {threads} {input.fa} <(seqtk sample {input.reads} {params.frac}) > {output.sam}"

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
        gt={input.vcf} \
        ref={input.panel} \
        map={input.gmap} \
        out={params.prefix};\n"
        "bcftools index {output.vcf};\n"

rule imp_index:
    input:
        ref=config["hg19_index"],
        vcf=IMP_BEAGLE_GT,
        idx=IMP_BEAGLE_GT + ".csi"
    output:
        fa=IMP_PG,
        idx1=IMP_PG+".amb",
        idx2=IMP_PG+".ann",
        idx3=IMP_PG+".bwt",
        idx4=IMP_PG+".pac",
        idx5=IMP_PG+".sa",
    threads: 16
    params:
        sample="{sample}",
        i="{i}"
    shell:
        "bcftools consensus -f {input.ref} -H {params.i} -s {params.sample} {input.vcf} > {output.fa} 2> /dev/null;\n"
        # "bowtie2-build --threads {threads} {output.fa} {output.fa}"
        "bwa index {output.fa} {output.fa}"

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
IMP_MERGED_ALNS       =   os.path.join(IMP_DIR, "{genotyper}/{filter}/merged.sam")
IMP_MERGED_SCORE      =   os.path.join(IMP_DIR, "{genotyper}/{filter}/merged.score")

rule imp_align_all:
    input:
        fa=IMP_PG,
        idx1=IMP_PG+".amb",
        reads=TEST_READS
    output:
        sam=IMP_UNLIFTED_ALNS
    threads: 16
    shell:
        # "bowtie2 -p {threads} -x {input.fa} -U {input.reads} > {output.sam}"
        "bwa mem -t {threads} {input.fa} {input.reads} > {output.sam}"

rule imp_lift_alns:
    input:
        sam=IMP_UNLIFTED_ALNS,
        vcf=IMP_BEAGLE_GT,
        lengths=REF_LENGTHS
    output:
        sam=IMP_LIFTED_ALNS
    params:
        i=lambda w: int(w.i) - 1,
        sample="{sample}",
        prefix=IMP_LIFTED_ALNS.replace(".sam","")
    shell:
        cmds["lift"]

rule imp_merge_alns:
    input:
        IMP_LIFTED_ALNS.replace("{i}", "1"),
        IMP_LIFTED_ALNS.replace("{i}", "2")
    output:
        IMP_MERGED_ALNS
    shell:
        "hts_utils/merge_sams {input} > {output}"
        

rule imp_score_alns:
    input:
        truth=TEST_TRUTH,
        sam=IMP_MERGED_ALNS
    output:
        IMP_MERGED_SCORE
    shell:
        cmds["score"]

#### REFERENCE ALIGNMENTS ######
HG19_ALNS              =   "{sample}/GRCh37/alns.sam"
HG19_SCORE      =   "{sample}/GRCh37/alns.score"

rule hg19_align_all:
    input:
        fa=config["hg19_index"],
        idx1=config["hg19_index"]+".amb",
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
MA_UNLIFTED_ALNS              =   "{sample}/ma/unlifted.sam"
MA_LIFTED_ALNS              =   "{sample}/ma/lifted.sam"
MA_SCORE      =   "{sample}/ma/lifted.score"

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
        idx1=MA_GENOME+".amb",
        idx2=MA_GENOME+".ann",
        idx3=MA_GENOME+".bwt",
        idx4=MA_GENOME+".pac",
        idx5=MA_GENOME+".sa",
    params:
    run:
        shell("bcftools consensus -f {input.fa} -H 1 {input.gt} > {output.fa}")
        # shell("bowtie2-build {output.fa} {output.fa}")
        shell("bwa index {output.fa} {output.fa}")

rule ma_align_all:
    input:
        fa=MA_GENOME,
        idx1=MA_GENOME+".amb",
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

rule ma_score_alns:
    input:
        truth=TEST_TRUTH,
        sam=MA_LIFTED_ALNS
    output:
        MA_SCORE
    shell: 
        cmds["score"]


###### PERSONAL GENOME #######
PERS_UNLIFTED_ALNS     =   "{sample}/personal/unlifted.{i}.sam"
PERS_LIFTED_ALNS       =   "{sample}/personal/lifted.{i}.sam"
PERS_MERGED_ALNS       =   "{sample}/personal/merged.sam"
PERS_MERGED_SCORE      =   "{sample}/personal/merged.score"

rule personal_index:
    input:
        fa=GENOME
    output:
        idx1=GENOME+".amb",
        idx2=GENOME+".ann",
        idx3=GENOME+".bwt",
        idx4=GENOME+".pac",
        idx5=GENOME+".sa",
    threads: 16
    shell:
        # "bowtie2-build --threads {threads} {input} {input}"
        "bwa index {input} {input}"
        
rule personal_align_all:
    input:
        fa=GENOME,
        reads=TEST_READS,
        idx1=GENOME+".amb",
    output:
        sam=PERS_UNLIFTED_ALNS
    threads: 16
    shell:
        cmds["align"]

rule personal_lift_alns:
    input:
        vcf=VCF,
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

rule personal_merge_alns:
    input:
        PERS_LIFTED_ALNS.replace("{i}", "1"),PERS_LIFTED_ALNS.replace("{i}","2")
    output:
        PERS_MERGED_ALNS
    shell:
        "hts_utils/merge_sams {input} > {output}"

rule personal_score_alns:
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
OTHER_UNLIFTED_ALNS     =   "{sample}/{other_sample}/unlifted.{i}.sam"
OTHER_LIFTED_ALNS       =   "{sample}/{other_sample}/lifted.{i}.sam"
OTHER_MERGED_ALNS       =   "{sample}/{other_sample}/merged.sam"
OTHER_MERGED_SCORE      =   "{sample}/{other_sample}/merged.score"
OTHER_SCORE             =   "{sample}/{other_sample}/alns.score"

rule other_align_all:
    input:
        fa=OTHER_GENOME,
        reads=TEST_READS,
        idx1=OTHER_GENOME+".amb",
    output:
        sam=OTHER_UNLIFTED_ALNS
    threads: 16
    shell:
        cmds["align"]
    
rule other_lift_alns:
    input:
        vcf=VCF,
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
        OTHER_LIFTED_ALNS.replace("{i}", "1"),
        OTHER_LIFTED_ALNS.replace("{i}","2")
    output:
        OTHER_MERGED_ALNS
    shell:
        "hts_utils/merge_sams {input} > {output}"

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
        ref_vcf=VCF,
        imp_vcf=IMP_BEAGLE_GT, 
        scores=IMP_MERGED_SCORE,
    output:
        "{sample}/{coverage}x/{genotyper}/{filter}/summary.txt",
    params:
        sample="{sample}",
        coverage="{coverage}",
        genotyper="{genotyper}",
        filter="{filter}",
        ploidy="diploid"
    run:
        # get hamm of imputed
        hamm1 = 0
        hamm2 = 0
        for line in shell("varcount/vcf_score --score-only -r {params.sample} -p {params.sample} {input.ref_vcf} {input.imp_vcf}", iterable=True):
            hamm1 = int(line)
            break
        for line in shell("varcount/vcf_score --flip-gt --score-only -r {params.sample} -p {params.sample} {input.ref_vcf} {input.imp_vcf}", iterable=True):
            hamm2 = int(line)
            break
        # get score of alignments
        c_alns = 0
        t_alns = 0
        for line in open(input.scores):
            if line[-2] == "1":
                c_alns += 1
            t_alns += 1
        with open(output[0], "w") as fp:
            to_write = "{sample}\t{coverage}\t{genotyper}\t{filter}\t{hamm1}\t{hamm2}\t{correct:.5f}\t{ploidy}\n".format(
                sample=params.sample,
                coverage=params.coverage,
                genotyper=params.genotyper,
                filter=params.filter,
                hamm1=hamm1,
                hamm2=hamm2,
                correct=float(c_alns)/t_alns,
                ploidy=params.ploidy
            )
            fp.write(to_write)

rule personal_summarize_exp:
    input:
        ref_vcf=VCF,
        scores=PERS_MERGED_SCORE,
    output:
        "{sample}/personal/summary.txt",
    params:
        sample="{sample}",
        coverage=config["total_coverage"],
        genotyper="PERSONAL",
        filter="PERSONAL",
        ploidy="diploid"
    run:
        hamm1 = 0
        hamm2 = 0
        c_alns = 0
        t_alns = 0
        for line in open(input.scores):
            if line[-2] == "1":
                c_alns += 1
            t_alns += 1
        with open(output[0], "w") as fp:
            to_write = "{sample}\t{coverage}\t{genotyper}\t{filter}\t{hamm1}\t{hamm2}\t{correct:.5f}\t{ploidy}\n".format(
                sample=params.sample,
                coverage=params.coverage,
                genotyper=params.genotyper,
                filter=params.filter,
                hamm1=hamm1,
                hamm2=hamm2,
                correct=float(c_alns)/t_alns,
                ploidy=params.ploidy
            )
            fp.write(to_write)

rule hg19_summarize_exp:
    input:
        ref_vcf=VCF,
        scores=HG19_SCORE,
    output:
        "{sample}/GRCh37/summary.txt",
    params:
        sample="{sample}",
        coverage="0",
        genotyper="GRCh37",
        filter="GRCh37",
        ploidy="haploid"
    run:
        hamm1 = 0
        hamm2 = 0
        talts = 0
        c_alns = 0
        t_alns = 0
        for line in shell("bcftools view -s {params.sample} {input.ref_vcf} | wc -l", iterable=True):
            talts = int(line)
            break
        for line in shell("bcftools view -s {params.sample} {input.ref_vcf} | bcftools view -i 'GT~\"0|\"'", iterable=True):
            hamm1 += 1
        for line in shell("bcftools view -s {params.sample} {input.ref_vcf} | bcftools view -i 'GT~\"|0\"'", iterable=True):
            hamm2 += 1
        hamm1 = talts - hamm1
        hamm2 = talts - hamm2
        for line in open(input.scores):
            if line[-2] == "1":
                c_alns += 1
            t_alns += 1

        with open(output[0], "w") as fp:
            to_write = "{sample}\t{coverage}\t{genotyper}\t{filter}\t{hamm1}\t{hamm2}\t{correct:.5f}\t{ploidy}\n".format(
                sample=params.sample,
                coverage=params.coverage,
                genotyper=params.genotyper,
                filter=params.filter,
                hamm1=hamm1,
                hamm2=hamm2,
                correct=float(c_alns)/t_alns,
                ploidy=params.ploidy
            )
            fp.write(to_write)

rule other_summarize_exp:
    input:
        ref_vcf=VCF,
        scores=OTHER_MERGED_SCORE,
    output:
        "{sample}/{other_sample}/summary.txt",
    params:
        other_sample="{other_sample}",
        sample="{sample}",
        coverage="0",
        genotyper="{other_sample}",
        filter="{other_sample}",
        ploidy="diploid"
    run:
        # get hamm of imputed
        hamm1 = 0
        hamm2 = 0
        for line in shell("varcount/vcf_score --score-only -r {params.sample} -p {params.other_sample} {input.ref_vcf} {input.ref_vcf}", iterable=True):
            hamm1 = int(line)
            break
        for line in shell("varcount/vcf_score --flip-gt --score-only -r {params.sample} -p {params.other_sample} {input.ref_vcf} {input.ref_vcf}", iterable=True):
            hamm2 = int(line)
            break
        # get score of alignments
        c_alns = 0
        t_alns = 0
        for line in open(input.scores):
            if line[-2] == "1":
                c_alns += 1
            t_alns += 1
        with open(output[0], "w") as fp:
            to_write = "{sample}\t{coverage}\t{genotyper}\t{filter}\t{hamm1}\t{hamm2}\t{correct:.5f}\t{ploidy}\n".format(
                sample=params.sample,
                coverage=params.coverage,
                genotyper=params.genotyper,
                filter=params.filter,
                hamm1=hamm1,
                hamm2=hamm2,
                correct=float(c_alns)/t_alns,
                ploidy=params.ploidy
            )
            fp.write(to_write)

rule ma_summarize_exp:
    input:
        ref_vcf=VCF,
        ma_vcf=MA_VCF,
        scores=MA_SCORE,
    output:
        "{sample}/ma/summary.txt",
    params:
        sample="{sample}",
        coverage="0",
        genotyper="MA",
        filter="MA",
        ploidy="haploid"
    run:
        # get hamm of imputed
        hamm1 = 0
        for line in shell("varcount/vcf_score --score-only -r {params.sample} -p MA {input.ref_vcf} {input.ma_vcf}", iterable=True):
            hamm1 = int(line)
            break
        hamm2 = hamm1
        # get score of alignments
        c_alns = 0
        t_alns = 0
        for line in open(input.scores):
            if line[-2] == "1":
                c_alns += 1
            t_alns += 1
        with open(output[0], "w") as fp:
            to_write = "{sample}\t{coverage}\t{genotyper}\t{filter}\t{hamm1}\t{hamm2}\t{correct:.5f}\t{ploidy}\n".format(
                sample=params.sample,
                coverage=params.coverage,
                genotyper=params.genotyper,
                filter=params.filter,
                hamm1=hamm1,
                hamm2=hamm2,
                correct=float(c_alns)/t_alns,
                ploidy=params.ploidy
            )
            fp.write(to_write)

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
