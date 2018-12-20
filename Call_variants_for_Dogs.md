# Call variants at all non-conserved and non-CDS regions for Basenji

Basenji bam file (mapped reads) location on Phoenix: `/data/acad/ysouilmi/dingo/01-data/Basenji/Basenji.bam`

CanFam (dog reference genome) location: `/data/acad/ysouilmi/dingo/CanisRef/Canis_familiaris.CanFam3.1.dna.toplevel.fasta`

## Restrict bed file of non-coding and non-conserved regions to only chromosome 1:

```
awk '$1 = 1 {OFS="\t" ; print $0}' Canis_familiaris.CanFam3.1.94.nonConserved_nonCoding.bed > chr1.bed
```

## Call variants at all positions in bed file

```
ref=/data/acad/ysouilmi/dingo/CanisRef/Canis_familiaris.CanFam3.1.dna.toplevel.fasta
bed=chr1.bed # see above
bam=/data/acad/ysouilmi/dingo/01-data/Basenji/Basenji.bam

bcftools mpileup --threads 4 \
        -f $ref \
        --ignore-RG \
        -q 30 -Q 30 \
        -R $bed \
        $bam
        | bcftools call --threads 4 \
        -m \
        --skip-variants indels \
        -O z \
        -o Basenji_chr1.vcf.gz
```

## Merge Basenji VCF with other dogs
```
bcftools merge \
    -O z \
    -o merged_chr1.vcf.gz \
    -R $bed \
    --threads 4 \
    -0 \
    Basenji_chr1.vcf.gz other_dogs_chr1.vcf.gz # remember to change file name where necessary

```