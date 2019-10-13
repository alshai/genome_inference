fp1 = open(snakemake.output[0], "w")
fp2 = open(snakemake.output[1], "w")
for line in open(snakemake.input[0], "r"):
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
