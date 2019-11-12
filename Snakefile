configfile: "snakefiles/config.yaml"

if 'threads' not in config:
    MAX_THREADS = 32
else:
    MAX_THREADS = config['threads']

wildcard_constraints:
    i="[12]"

rule all:
    input:
        expand("{sample}/pg.{cov}x.{ext}",
               sample=config["samples"],
               cov=config["coverages"],
               ext=["sam"]
        )

rule imp_align_subset:
    input:
        fq1="{sample}/reads.1.fq",
        fq2="{sample}/reads.2.fq",
        fa=config["ref_fa"],
        idx1=config["ref_idx_pre"] + ".1.bt2"
    params:
        frac=lambda w: float(w.cov) / config["total_coverage"], # multiply by 2 for diploid?
    output:
        tempsam=temp("{sample}/pg.{cov}x.subsample.temp.bam"),
        sam="{sample}/pg.{cov}x.subsample.bam",
        sam_idx="{sample}/pg.{cov}x.subsample.bam.bai",
    threads: MAX_THREADS
    log: "{sample}/pg.{cov}x.align_subset.log"
    benchmark: "{sample}/pg.{cov}x.align_subset.benchmark"
    shell:
        "bowtie2 -p {threads} -x {input.fa} "
        "-U <(seqtk sample -s100 <(cat {input.fq1} {input.fq2}) {params.frac}) "
        "| samtools view -@ {threads} -b > {output.tempsam};\n"
        "samtools sort -@{threads} {output.tempsam} > {output.sam};\n"
        "samtools index {output.sam};"

# we only count for chr21!!
rule imp_count:
    input:
        sam="{sample}/pg.{cov}x.subsample.bam",
        sam_idx="{sample}/pg.{cov}x.subsample.bam.bai",
        vcf=config["ref_panel"]
    output:
        vcf="{sample}/pg.{cov}x.counts.vcf.gz"
    params:
        sample="{sample}",
        gt="likelihood_naive"
    log: "{sample}/pg.{cov}x.count.log"
    benchmark: "{sample}/pg.{cov}x.count.benchmark"
    shell:
        "varcount/varcount -g{params.gt} -s {params.sample} {input.vcf} <(samtools view -u -q 40 {input.sam} 21) | "
        "bcftools sort -O z > {output.vcf};\n"
        "bcftools index {output.vcf}"

rule imp_filter_counts:
    input:
        vcf="{sample}/pg.{cov}x.counts.vcf.gz"
    output:
        vcf="{sample}/pg.{cov}x.counts.filtered.vcf.gz",
    log: "{sample}/pg.{cov}x.count_filter.log"
    benchmark: "{sample}/pg.{cov}x.count_filter.benchmark"
    shell:
        "bcftools view -i 'GT~\"A\"' {input.vcf} | bcftools sort -Oz > {output.vcf};\n"
        "bcftools index {output.vcf};"

rule imp_beagle_impute:
    input:
        vcf="{sample}/pg.{cov}x.counts.filtered.vcf.gz",
        panel=config["ref_panel"],
        gmap=config['beagle_map'],
        beagle_jar=config["beagle_jar"]
    output:
        vcf="{sample}/pg.{cov}x.imputed.vcf.gz",
        idx="{sample}/pg.{cov}x.imputed.vcf.gz.csi",
    params:
        prefix="{sample}/pg.{cov}x.imputed",
    benchmark: "{sample}/pg.{cov}x.impute.benchmark"
    threads:
        32
    log: "{sample}/pg.{cov}x.impute.log"
    shell:
        "java -jar {input.beagle_jar} \
        nthreads={threads} \
        gt={input.vcf} \
        ref={input.panel} \
        map={input.gmap} \
        out={params.prefix};\n"
        "bcftools index {output.vcf};\n"

# outputs a filterd, sorted vcf file
rule imp_beagle_filter:
    input:
        vcf="{sample}/pg.{cov}x.imputed.vcf.gz",
        vcf_idx="{sample}/pg.{cov}x.imputed.vcf.gz.csi"
    output:
        vcf="{sample}/pg.{cov}x.imputed.filtered.vcf.gz",
        vcf_idx="{sample}/pg.{cov}x.imputed.filtered.vcf.gz.csi"
    log:
        "{sample}/pg.{cov}x.impute_filter.log"
    benchmark: "{sample}/pg.{cov}x.impute_filter.benchmark"
    shell:
        "bcftools norm -Ou -d indels {input.vcf} | "
        "bcftools view -Ou -i 'GT~\"A\"' | "
        "bcftools sort -Oz > {output.vcf};\n"
        "bcftools index --force {output.vcf};\n"

#  bowtie2 indexing and alignment
rule imp_index:
    input:
        ref=config["ref_fa"],
        vcf="{sample}/pg.{cov}x.imputed.filtered.vcf.gz",
        idx="{sample}/pg.{cov}x.imputed.filtered.vcf.gz.csi"
    output:
        fa="{sample}/pg.{cov}x.{i}.fa",
        idx1="{sample}/pg.{cov}x.{i}.fa.1.bt2",
        idx2="{sample}/pg.{cov}x.{i}.fa.2.bt2",
        idx3="{sample}/pg.{cov}x.{i}.fa.3.bt2",
        idx4="{sample}/pg.{cov}x.{i}.fa.4.bt2",
        idx5="{sample}/pg.{cov}x.{i}.fa.rev.1.bt2",
        idx6="{sample}/pg.{cov}x.{i}.fa.rev.2.bt2"
    threads: MAX_THREADS
    params:
        sample="{sample}",
        i="{i}"
    benchmark: "{sample}/pg.{cov}x.{i}.bowtie2_build.benchmark"
    shell:
        "bcftools consensus -f {input.ref} -H {params.i} -s {params.sample} {input.vcf} > {output.fa};\n"
        "bowtie2-build --threads {threads} {output.fa} {output.fa};\n"

rule contig_lengths:
    input:
        config["ref_panel"]
    output:
        "{sample}/contig_lengths"
    shell:
        "bcftools view -h {input} | grep '##contig' | sed -E 's/##contig=<ID=(.*),assembly=.*,length=([0-9]+)>/\\1\t\\2/' > {output}"

rule imp_align_all:
    input:
        fa="{sample}/pg.{cov}x.{i}.fa",
        idx1="{sample}/pg.{cov}x.{i}.fa.1.bt2",
        idx2="{sample}/pg.{cov}x.{i}.fa.2.bt2",
        idx3="{sample}/pg.{cov}x.{i}.fa.3.bt2",
        idx4="{sample}/pg.{cov}x.{i}.fa.4.bt2",
        idx5="{sample}/pg.{cov}x.{i}.fa.rev.1.bt2",
        idx6="{sample}/pg.{cov}x.{i}.fa.rev.2.bt2",
        fq1="{sample}/reads.1.fq",
        fq2="{sample}/reads.2.fq",
    output:
        tempsam=temp("{sample}/pg.{cov}x.{i}.unlifted.tmp.bam"),
        sam="{sample}/pg.{cov}x.{i}.unlifted.bam",
        sam_idx="{sample}/pg.{cov}x.{i}.unlifted.bam.bai",
    threads: MAX_THREADS
    benchmark: "{sample}/pg.{cov}x.{i}.bowtie2_align.benchmark"
    shell:
        "bowtie2 -p {threads} -x {input.fa} -U {input.fq1},{input.fq2} | samtools view -b -@ {threads} > {output.tempsam};\n"
        "samtools sort -@{threads} {output.tempsam} > {output.sam};\n"
        "samtools index -@ {threads} {output.sam};\n"

# we only lift alignments from chr21!!
rule imp_lift_all:
    input:
        sam="{sample}/pg.{cov}x.{i}.unlifted.bam",
        sam_idx="{sample}/pg.{cov}x.{i}.unlifted.bam.bai",
        vcf="{sample}/pg.{cov}x.imputed.filtered.vcf.gz",
        ref_panel=config["ref_panel"],
        contig_lengths = "{sample}/contig_lengths"
    output:
        sam="{sample}/pg.{cov}x.{i}.lifted.sam"
    params:
        i=lambda w: int(w.i) - 1,
        sample="{sample}",
    benchmark:
        "{sample}/pg.{cov}x.{i}.lift.benchmark"
    log:
        "{sample}/pg.{cov}x.{i}.lift.log"
    shell:
        "~/.local/bin/time -v liftover/liftover lift "
            "-v {input.vcf} "
            "-a <(samtools view -u {input.sam} 21) "
            "-k  {input.contig_lengths} "
            "--haplotype {params.i} "
            "-s {params.sample} "
            "> {output.sam} 2>{log}"


rule imp_merge_alns:
    input:
        a="{sample}/pg.{cov}x.1.lifted.sam",
        b="{sample}/pg.{cov}x.2.lifted.sam",
        fa=config["ref_fa"]
    output:
        "{sample}/pg.{cov}x.sam"
    log:
        "{sample}/pg.{cov}x.merge.log"
    benchmark:
        "{sample}/pg.{cov}x.merge.benchmark"
    shell:
        "hts_utils/merge_sams {input.a} {input.b} | samtools calmd - {input.fa} > {output}"

# vg indexing and alignment

rule imp_vg_construct:
    input:
        fa=config["ref_fa"],
        vcf="{sample}/pg.{cov}x.imputed.filtered.vcf.gz",
    output:
        vg="{sample}/pg.{cov}x.vg",
    threads:
        MAX_THREADS
    log:
        "{sample}/pg.{cov}x.vg_construct.log"
    benchmark:
        "{sample}/pg.{cov}x.vg_construct.benchmark"
    shell: 
        "~/.local/bin/time -v vg construct -t {threads} -r {input.fa} -v {input.vcf} -a > {output.vg} 2>{log}"

rule imp_vg_index:
    input:
        vcf="{sample}/pg.{cov}x.imputed.filtered.vcf.gz",
        vg="{sample}/pg.{cov}x.vg",
    output:
        xg="{sample}/pg.{cov}x.xg",
        gcsa="{sample}/pg.{cov}x.gcsa"
    log:
        "{sample}/pg.{cov}x.vg_index.log"
    benchmark:
        "{sample}/pg.{cov}x.vg_index.benchmark"
    params:
        tmpdir="{sample}"
    threads: MAX_THREADS
    shell: 
        "~/.local/bin/time -v vg index -b {params.tmpdir}/ -t {threads} -x {output.xg} -g {output.gcsa} {input.vg} 2>{log}"

rule imp_vg_map:
    input:
        vg="{sample}/pg.{cov}x.vg",
        xg="{sample}/pg.{cov}x.xg",
        gcsa="{sample}/pg.{cov}x.gcsa",
        fq1="{sample}/reads.1.fq",
        fq2="{sample}/reads.2.fq"
    params:
        prefix="{sample}/pg.{cov}x",
    output:
        "{sample}/pg.{cov}x.gam"
    threads:
        MAX_THREADS
    log:
        "{sample}/pg.{cov}x.vg_map.log"
    benchmark:
        "{sample}/pg.{cov}x.vg_map.benchmark"
    shell:
        "~/.local/bin/time -v vg map -t {threads} -d {params.prefix} -f {input.fq1} -f {input.fq2} > {output} 2>{log}"

rule imp_vg_surject:
    input:
        gam="{sample}/pg.{cov}x.gam",
        vg="{sample}/pg.{cov}x.vg",
        xg="{sample}/pg.{cov}x.xg",
        gcsa="{sample}/pg.{cov}x.gcsa",
    params:
        prefix="{sample}/pg.{cov}x",
    output:
        sam="{sample}/pg.{cov}x.vg.bam"
    threads:
        MAX_THREADS
    log:
        "{sample}/pg.{cov}x.vg_surject.log"
    benchmark:
        "{sample}/pg.{cov}x.vg_surject.benchmark"
    shell:
        "~/.local/bin/time -v vg surject -b -t {threads} -x {input.xg} {input.gam} > {output.sam} 2>{log}"
