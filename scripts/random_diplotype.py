prob_alt = float(snakemake.params.prob_alt)
vfp = open(snakemake.output.vcf.replace(".gz", ""), "w")
for line in shell("bcftools view {input.vcf}", iterable=True):
    if "#CHROM" in line:
        vfp.write('\t'.join(line.strip().split()[:8]) + "\tFORMAT\trandom" + snakemake.params.j + "\n") 
    elif line[0] == "#":
        vfp.write(line + "\n")
    else:
        fields = line.strip().split('\t')
        a1 = '0'
        a2 = '0'
        if random.random() < prob_alt:
            a1 = '1'
        if random.random() < prob_alt:
            a2 = '1'
        if a1 == '1' or a2 == '1':
            vfp.write('\t'.join(fields[0:8] + ['GT\t' + a1 + '|' + a2 + '\n']))
vfp.close()
shell("bcftools view -Oz " + snakemake.output.vcf.replace(".gz", "") + " > {output.vcf}")
shell("bcftools index --force {output.vcf}")
shell("bcftools consensus -f {input.fa} -H 1 {output.vcf} > {output.fa1}")
shell("bcftools consensus -f {input.fa} -H 2 {output.vcf} > {output.fa2}")
shell("bowtie2-build {output.fa1} {output.fa1}")
shell("bowtie2-build {output.fa2} {output.fa2}")
