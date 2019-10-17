cmds = {
    'align': 
        # "bowtie2 -p {threads} -x {input.fa} -U {input.reads} > {output.sam}",
        "bwa mem -t {threads} {input.fa} {input.reads} > {output.sam}",
    'score': 
    "hts_utils/score_sam -b 25 {input.truth} {input.sam} > {output}" , # TODO: DEAL WITH THIS WORKAROUND
    'lift': 
        "liftover/liftover lift -v {input.vcf} -a {input.sam} -k {input.lengths} --haplotype {params.i} -s {params.sample} -p {params.prefix}",
    'align_and_lift' :
        "bowtie2 \
                -p {threads} \
                -x {input.fa} \
                -U {input.reads}  | \
        liftover/liftover lift \
        -v {input.vcf} \
        -a {input.sam} \
        -k {input.lengths} \
        --haplotype {params.i} \
        -s {params.sample} \
        -p {params.prefix}"
}
