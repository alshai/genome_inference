def vcf_score(input,output,params,shell):
    dist1 = 0
    dist2 = 0
    for line in shell("varcount/vcf_score --score-only -r {params.s1} -p {params.s2} {input.vcf1} {input.vcf2}", iterable=True):
        dist1 = int(line)
        break
    for line in shell("varcount/vcf_score --flip-gt --score-only -r {params.s1} -p {params.s2} {input.vcf1} {input.vcf2}", iterable=True):
        dist2 = int(line)
        break
    return min(dist1, dist2)

def aln_score(input,output,params,shell):
    c_alns = 0
    t_alns = 0
    for line in open(input.scores):
        line = line.strip().split()
        if line[5] == "1":
            c_alns += 1
        t_alns += 1
    return float(c_alns)/t_alns

def ref_bias(input, output, params, shell):
    refc = 0
    altc = 0
    for line in shell("~/geninf/varcount/varcount <(bcftools norm -Ou -d all {input.vcf1} | bcftools view -Ou -V indels -g het -m2 -M2) {input.alns} | bcftools view -H - ", iterable=True):
        r,a = list(map(int, line.strip().split("\t")[9].split(",")))
        refc += r
        altc += a
    return refc, altc

def summarize(input, output, params, shell):
    rc, ac = ref_bias(input, output, params, shell)
    with open(output.out, "w") as fp:
        to_write = "{sample}\t{exp}\t{ploidy}\t{coverage}\t{imp_input}\t{dist}\t{correct:.5f}\t{het_ref_count}\t{het_alt_count}\t{ref_over_alt:.5f}\n".format(
            sample=params.s1,
            exp=params.experiment,
            ploidy=params.ploidy,
            coverage=params.coverage,
            imp_input=params.imp_input,
            dist=vcf_score(input,output,params,shell),
            correct=aln_score(input,output,params,shell),
            het_ref_count=rc,
            het_alt_count=ac,
            ref_over_alt=float(rc)/ac
        )
        fp.write(to_write)
