from scripts.summarize import vcf_score

def dist_windows(input, output, params, shell):
    with open(output.out2, "w") as fp:
        for w in [2**i for i in range(5,25)]:
            to_write = "{sample}\t{exp}\t{ploidy}\t{coverage}\t{imp_input}\t{w}\t{dist}\n".format(
                sample=params.s1,
                exp=params.experiment,
                ploidy=params.ploidy,
                coverage=params.coverage,
                imp_input=params.imp_input,
                w=w,
                dist=vcf_score(input,output,params,shell, w),
            )
            fp.write(to_write)

