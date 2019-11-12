# Impute a personalized genome to improve read alignment

A Snakemake workflow (https://snakemake.readthedocs.io)

Use some of your sample's genomic data with a haplotype reference panel to
create a personalized genome reference.

## Usage:

``` 
snakemake -j [threads]
```

OR 

```
snakemake -j [threads] {sample}/pg.{coverage}x.sam
```

`{sample}` is the name of the sample from which your reads originate.

`{coverage}` is the amount of information from your sample that you wish to use
to impute the personalized reference.

## Prerequisites:

- Paired end short reads files should be located at `{sample}/reads.1.fq`
  and `{sample}/reads.2.fq`.

- `config.yaml` must be modified to include, among other things:

    + the reference genome location

    + haplotype reference panel VCF location (download the 1000 Genomes
      Reference Panel from
      [here](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)). The
      coordinates of this VCF must be with respect to the specified reference
      genome.

    + location of the Beagle JAR file for imputation (download
      [here](https://faculty.washington.edu/browning/beagle/beagle.08Nov19.3ec.jar))

## The workflow outputs:

- two fasta files representing each haploid of an approximate personalized
  diploid genome (`{sample}/pg.{coverage}x.1.fa`, `{sample}/pg.{coverage}x.2.fa`)

- bowtie2 indexes for each of the two fasta files

- Merged alignments against both haploid indexes. For each read, only the
  better alignment against the two haploids is reported. These alignments are
  lifted over to the specified reference genome (`{sample}/pg.{coverage}x.sam`)
