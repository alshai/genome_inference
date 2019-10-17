import glob
import os.path

def get_reads_set(fname, is_correct):
    reads = set()
    for line in open(fname):
        fields = line.strip().split('\t')
        if int(fields[5]) == is_correct:
            reads.add(fields[0])
    return reads

def get_reads_dict(fname, is_correct):
    reads = dict()
    for line in open(fname):
        fields = line.strip().split('\t')
        if int(fields[5]) == is_correct:
            reads[fields[0]] = (fields[1], int(fields[2]))
    return reads


ref_fname = snakemake.input.ref_score
ref_incorrect = get_reads_set(ref_fname, 0)

imp_fname = snakemake.input.imp_score
imp_correct = get_reads_set(imp_fname, 1)

imp_dict = get_reads_dict(imp_fname, 1)

out_arr = []
ixn = imp_correct & ref_incorrect
for r in ixn:
    pos = imp_dict[r]
    out_arr.append((pos[0], pos[1], pos[1]+150, r))
out_arr.sort(key=lambda x:x[1])

with open(snakemake.output[0], "w") as fp:
    for c, s, e, r in out_arr:
        fp.write("chr{chr}\t{s}\t{e}\t{r}\n".format(chr=c, s=s, e=e, r=r))
