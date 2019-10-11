import math 
import os

configfile: "config.yaml"


cmds = {
    'align':
        "bowtie2 -p {threads} -x {input.fa} -U {input.reads} > {output.sam}",
    'score': 
        "hts_utils/score_sam {input.truth} {input.sam} > {output}" ,
    'lift': 
        "liftover/liftover lift -v {input.vcf} -a {input.sam} -k {input.lengths} --haplotype {params.i} -s {params.sample} -p {params.prefix}"
}


IMP_LIFTED_SCORE    =  "{sample}/{coverage}x/{genotyper}/{filter}/merged.score"
PERS_MERGED_SCORE   =  "{sample}/personal/merged.score"
HG19_SCORE          =  "{sample}/hg19/alns.score"

rule all:
    input:
        expand("{sample}/hg19/summary.txt", sample=config["samples"]),
        expand("{sample}/personal/summary.txt", sample=config["samples"]),
        expand("{sample}/{coverage}x/{genotyper}/{filter}/summary.txt",
               genotyper=["likelihood_naive"],
               filter=["aa", 'ar', 'aa+ar', 'aa+ar+some_rr'],
               sample=config["samples"], 
               coverage=config["coverage"]
              ),
        expand(IMP_LIFTED_SCORE, 
               genotyper=["likelihood_naive"],
               filter=["aa", 'ar', 'aa+ar', 'aa+ar+some_rr'],
               sample=config["samples"], 
               coverage=config["coverage"]),
        expand(PERS_MERGED_SCORE, sample=config["samples"]),
        expand(HG19_SCORE, sample=config["samples"])
###### PREPARE DATA FOR LATER STAGES ############

VCF         =   config["vcf"].replace(".vcf.gz", ".filtered.vcf.gz")
REF_PANEL   =   config["vcf"].replace(".vcf.gz", ".filtered.leave_out.vcf.gz")
READS       =   "{sample}/reads.fq.gz"
TRUE_ALNS =   READS.replace(".fq.gz", ".lifted.bam")
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
        "bcftools consensus -f {input.ref} -H {params.i} -s {params.sample} {input.vcf} > {output.fa};"

rule simulate_reads:
    input:
        "{sample}/genome1.fa",
        "{sample}/genome2.fa"
    output:
        fa=temp("{sample}/full_genome.fa"),
        fq=READS,
        sam=READS.replace(".fq.gz", ".sam")
    params:
        nreads=config["nreads"],
        rlen=config["read_length"]
    threads: 16
    shell:
        "cat {input} | awk 'BEGIN {{i=1}} /^>/ {{print \">genome\"i++}} /^[^>]/ {{print}}' > {output.fa};\n"
        "mason_simulator --num-threads {threads} "
                         "-ir {output.fa} "
                         "-n {params.nreads} "
                         "--illumina-read-length {params.rlen} "
                         "-o {output.fq} "
                         "-oa {output.sam}"

rule split_true_read_alignments:
    input:
        sam=READS.replace(".fq.gz", ".sam")
    output:
        sam1=READS.replace(".fq.gz", ".1.sam"),
        sam2=READS.replace(".fq.gz", ".2.sam")
    run:
        fp1 = open(output[0], "w")
        fp2 = open(output[1], "w")
        for line in open(input[0], "r"):
            if line[0] == "@":
                fp1.write(line)
                fp2.write(line)
            else:
                fields = line.split()
                if fields[2] == "genome1":
                    fp1.write(line)
                elif fields[2] == "genome2":
                    fp2.write(line)
        fp1.close()
        fp2.close()

rule lift_true_read_alignments:
    input:
        vcf=VCF,
        sam=READS.replace(".fq.gz", ".{i}.sam"),
    output:
        sam=READS.replace(".fq.gz", ".{i}.lifted.sam")
    params:
        sample = "{sample}",
        h = lambda w: int(w.i) - 1,
        i = "{i}"
    shell:
        "liftover/liftover lift -v {input.vcf} -a {input.sam} -s {params.sample} --haplotype {params.h} -n <(echo '21	genome{params.i}') > {output}"
        
rule concat_lifted_truth:
    input:
        one=READS.replace(".fq.gz", ".1.lifted.sam"),
        two=READS.replace(".fq.gz", ".2.lifted.sam")
    output:
        TRUE_ALNS
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
IMP_SUBSET_SAM        =   "{sample}/{coverage}x/subset_v_hg19.sam"
IMP_COUNTS            =   "{sample}/{coverage}x/counts.{genotyper}.vcf.gz"
IMP_FILTERED_COUNTS   =   "{sample}/{coverage}x/{genotyper}/{filter}/counts.vcf.gz"
IMP_BEAGLE_GT         =   "{sample}/{coverage}x/{genotyper}/{filter}/imputed.vcf.gz"
IMP_PG                =   "{sample}/{coverage}x/{genotyper}/{filter}/pg.{i}.fa"
IMP_LIFT              =   "{sample}/{coverage}x/{genotyper}/{filter}/pg.{i}.lft"

rule imp_align_subset:
    input:
        reads=READS,
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
        "bcftools consensus -f {input.ref} -H {params.i} -s {params.sample} {input.vcf} > {output.fa};\n"
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
IMP_MERGED_ALNS       =   os.path.join(IMP_DIR, "{genotyper}/{filter}/merged.sam")
IMP_LIFTED_SCORE      =   os.path.join(IMP_DIR, "{genotyper}/{filter}/merged.score")

rule imp_align_all:
    input:
        fa=IMP_PG,
        idx1=IMP_PG+".1.bt2",
        idx2=IMP_PG+".2.bt2",
        idx3=IMP_PG+".3.bt2",
        idx4=IMP_PG+".4.bt2",
        idx5=IMP_PG+".rev.1.bt2",
        idx6=IMP_PG+".rev.2.bt2",
        reads=READS
    output:
        sam=IMP_UNLIFTED_ALNS
    threads: 16
    shell:
        "bowtie2 -p {threads} -x {input.fa} -U {input.reads} > {output.sam}"

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
        truth=TRUE_ALNS,
        sam=IMP_MERGED_ALNS
    output:
        IMP_LIFTED_SCORE
    shell:
        cmds["score"]

#### REFERENCE ALIGNMENTS ######
HG19_ALNS              =   "{sample}/hg19/alns.sam"
HG19_SCORE      =   "{sample}/hg19/alns.score"

rule hg19_align_all:
    input:
        fa=config["hg19_index"],
        idx1=config["hg19_index"]+".1.bt2",
        idx2=config["hg19_index"]+".2.bt2",
        idx3=config["hg19_index"]+".3.bt2",
        idx4=config["hg19_index"]+".4.bt2",
        idx5=config["hg19_index"]+".rev.1.bt2",
        idx6=config["hg19_index"]+".rev.2.bt2",
        reads=READS
    output:
        sam=HG19_ALNS
    threads: 16
    shell:
        cmds["align"]

rule hg19_score_alns:
    input:
        truth=TRUE_ALNS,
        sam=HG19_ALNS
    output:
        HG19_SCORE
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
        reads=READS,
        idx1=GENOME+".1.bt2",
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
        truth=TRUE_ALNS,
        sam=PERS_MERGED_ALNS
    threads: 16
    output:
        PERS_MERGED_SCORE
    shell:
        cmds["score"]


rule imp_summarize_exp:
    input:
        ref_vcf=VCF,
        imp_vcf=IMP_BEAGLE_GT, 
        scores=IMP_LIFTED_SCORE,
    output:
        "{sample}/{coverage}x/{genotyper}/{filter}/summary.txt",
    params:
        sample="{sample}",
        coverage="{coverage}",
        genotyper="{genotyper}",
        filter="{filter}"
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
            to_write = "{sample}\t{coverage}\t{genotyper}\t{filter}\t{hamm1}\t{hamm2}\t{correct:.5f}\n".format(
                sample=params.sample,
                coverage=params.coverage,
                genotyper=params.genotyper,
                filter=params.filter,
                hamm1=hamm1,
                hamm2=hamm2,
                correct=float(c_alns)/t_alns
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
        coverage="30x",
        genotyper="NA",
        filter="NA"
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
            to_write = "{sample}\t{coverage}\t{genotyper}\t{filter}\t{hamm1}\t{hamm2}\t{correct:.5f}\n".format(
                sample=params.sample,
                coverage=params.coverage,
                genotyper=params.genotyper,
                filter=params.filter,
                hamm1=hamm1,
                hamm2=hamm2,
                correct=float(c_alns)/t_alns
            )
            fp.write(to_write)

rule hg19_summarize_exp:
    input:
        ref_vcf=VCF,
        scores=HG19_SCORE,
    output:
        "{sample}/hg19/summary.txt",
    params:
        sample="{sample}",
        coverage="30x",
        genotyper="NA",
        filter="NA"
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
            to_write = "{sample}\t{coverage}\t{genotyper}\t{filter}\t{hamm1}\t{hamm2}\t{correct:.5f}\n".format(
                sample=params.sample,
                coverage=params.coverage,
                genotyper=params.genotyper,
                filter=params.filter,
                hamm1=hamm1,
                hamm2=hamm2,
                correct=float(c_alns)/t_alns
            )
            fp.write(to_write)

rule gather_summaries:
    input:
        "{sample}/{everything}/summary.txt"
    output:
        "{sample}/summaries.txt"
    shell:
        "cat {input} > {output}"
