# *Phytophthora rubi*
Commands used in the analysis of *P. rubi* genomes
SCRP249 SCRP324 SCRP333
====================

Commands used during analysis of *Phytophthora rubi* genomes. Note - all this work was performed in the directory: /home/groups/harrisonlab/project_files/phytophthora_rubi

The following is a summary of the work presented in this Readme:
Data organisation:
  * Preparing data  
Draft Genome assembly
  * Data qc
  * Genome assembly
  * Repeatmasking
  * Gene prediction
  * Functional annotation
Genome analysis
  * Homology between predicted genes & published effectors


#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```bash
```


#Data qc

programs: fastqc fastq-mcf kmc

Data quality was visualised using fastqc:


```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for RawData in $(ls raw_dna/paired/P.rubi/$Strain/*/*.fastq.gz)
    do
        echo $RawData
        ProgDir=/home/adamst/git_repos/tools/seq_tools/dna_qc
        qsub $ProgDir/run_fastqc.sh $RawData
    done
done
```

Trimming was performed on data to trim adapters from sequences and remove poor quality data.
This was done with fastq-mcf


```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    echo $Strain
    Read_F=$(ls raw_dna/paired/P.rubi/$Strain/F/*.fastq.gz)
    Read_R=$(ls raw_dna/paired/P.rubi/$Strain/R/*.fastq.gz)
    IluminaAdapters=/home/adamst/git_repos/tools/seq_tools/ncbi_adapters.fa
    ProgDir=/home/adamst/git_repos/tools/seq_tools/rna_qc
    echo $Read_F
    echo $Read_R
    qsub $ProgDir/rna_qc_fastq-mcf.sh $Read_F $Read_R $IluminaAdapters DNA
done
```

Data quality was visualised once again following trimming:

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for RawData in $(ls qc_dna/paired/P.rubi/$Strain/*/*.fq.gz)
    do
        echo $RawData
        ProgDir=/home/adamst/git_repos/tools/seq_tools/dna_qc
        qsub $ProgDir/run_fastqc.sh $RawData
    done
done
```


kmer counting was performed using kmc.
This allowed estimation of sequencing depth and total genome size:

```bash  
for Strain in SCRP249 SCRP324 SCRP333
do
    echo $Strain
    Trim_F1=$(ls qc_dna/paired/P.rubi/$Strain/F/*.fq.gz)
    Trim_R1=$(ls qc_dna/paired/P.rubi/$Strain/R/*.fq.gz)
    ProgDir=/home/adamst/git_repos/tools/seq_tools/dna_qc
    qsub $ProgDir/kmc_kmer_counting.sh $Trim_F1 $Trim_R1
done
```

** Estimated Genome Size is:
SCRP249: 105406182
SCRP324: 104021199
SCRP333: 102239327 **

** Esimated Coverage is:
SCRP249: 26
SCRP324: 27
SCRP333: 27 **

Target coverage is 20.
The ones at value 5 are errors from filtering of error kmers, estimate from plots follow in ().

##Coverage assessed using count_nucl.pl

```bash
for DataDir in $(ls -d qc_dna/paired/P.rubi/*)
do
    F_Read=$(ls $DataDir/F/*.gz)
    R_Read=$(ls $DataDir/R/*.gz)
    Strain=$(echo $DataDir | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $DataDir | rev | cut -f2 -d '/' | rev)
    WorkDir=tmp_dir/$Strain
    mkdir -p $WorkDir
    cp -r $F_Read $WorkDir
    cp -r $R_Read $WorkDir
    cd $WorkDir
    Read1=*R1*
    Read2=*R2*
    gunzip $Read1
    gunzip $Read2
    Sub1=*R1*.fq
    Sub2=*R2*.fq
    echo "$Organism - $Strain"
    count_nucl.pl -i $Sub1 -i $Sub2 -g 96
    cd /home/groups/harrisonlab/project_files/phytophthora_rubi
done
```

```
Estimated coverage is:
SCRP249: 51.00
SCRP324: 50.36
SCRP333: 49.79
```

# Assembly
Assembly was performed using: Spades and Discovar

## Spades Assembly

For single runs

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    F_Read=$(ls qc_dna/paired/P.rubi/$Strain/F/*.fq.gz)
    R_Read=$(ls qc_dna/paired/P.rubi/$Strain/R/*.fq.gz)
    CovCutoff='10'
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/spades
    Species=$(echo $F_Read | rev | cut -f4 -d '/' | rev)
    OutDir=assembly/spades/$Species/$Strain
    echo $Species
    echo $Strain
    qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct $CovCutoff
done
```

### Quast

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    OutDir=$(ls -d assembly/spades/P.rubi/$Strain/filtered_contigs)
    AssFiltered=$OutDir/contigs_min_500bp.fasta
    AssRenamed=$OutDir/contigs_min_500bp_renamed.fasta
    echo $AssFiltered
    printf '.\t.\t.\t.\n' > editfile.tab
    $ProgDir/remove_contaminants.py --inp $AssFiltered --out $AssRenamed --coord_file editfile.tab
    rm editfile.tab
done
```

### QUAST used to summarise assembly statistics

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Strain in SCRP249 SCRP324 SCRP333
do
    for Assembly in $(ls assembly/spades/P.rubi/$Strain/filtered_contigs/*_500bp_renamed.fasta)
    do
        Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
        Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
        OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
        qsub $ProgDir/sub_quast.sh $Assembly $OutDir
    done
done
```

```
N50:
SCRP249: 16613
SCRP324: 19133
SCRP333: 16946

L50:
SCRP249: 1233
SCRP324: 1057
SCRP333: 1210

Number of contigs > 1kb:
SCRP249: 9215
SCRP324: 9110
SCRP333: 8986
```

####Run Deconseq to remove contaminents

Contigs were identified that had BLAST hits to non-phytophthora genomes

```bash
for Assembly in $(ls assembly/spades/P.rubi/*/filtered_contigs/contigs_min_500bp_renamed.fasta)
do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    # Exclude_db="bact,virus,hsref"
    Exclude_db="paenibacillus"
    Good_db="phytoph"
    AssemblyDir=$(dirname $Assembly)
    # OutDir=$AssemblyDir/../deconseq
    OutDir=$AssemblyDir/../deconseq_Paen
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
done
```

Results were summarised using the following commands

```bash
# for File in $(ls assembly/spades/P.*/*/deconseq/log.txt); do
for File in $(ls assembly/spades/P.*/*/deconseq_Paen/log.txt)
do
    Name=$(echo $File | rev | cut -f3 -d '/' | rev)
    Good=$(cat $File |cut -f2 | head -n1 | tail -n1)
    Both=$(cat $File |cut -f2 | head -n2 | tail -n1)
    Bad=$(cat $File |cut -f2 | head -n3 | tail -n1)
    printf "$Name\t$Good\t$Both\t$Bad\n"
done
```

```
SCRP249	14020	3	1
SCRP324	13943	4	36
SCRP333	13677	2	1
```

Contaminent organisms identified by NCBI BLAST

```
SCRP249
PhiX

SCRP324
Paenibacillus spp (inc. Paenibacillus xylanexedens - PAMC 22703 https://www.ncbi.nlm.nih.gov/nucleotide/1120692005?report=genbank&log$=nuclalign&blast_rank=1&RID=R3JYNVSY014)
PhiX

SCRP333
PhiX
```

Assembly stats were collected on filtered assemblies

```bash
for Assembly in $(ls assembly/spades/P.*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta)
do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    OutDir=$(dirname $Assembly)
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

Assembly size was summarised and compared to previous assembly results

```bash
for Assembly in $(ls assembly/spades/P.*/*/deconseq_Paen/report.tsv)
do  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/'| rev)
    Size=$(cat $Assembly | grep 'Total length' | head -n1 | cut -f2)
    OldAssembly=$(ls assembly/spades/P.*/$Strain/filtered_contigs*/report.tsv)
    OldSize=$(cat $OldAssembly | grep 'Total length' | head -n1 | cut -f2)
    printf "$Strain\t$Size\t$OldSize\n"
done
```

```
SCRP249	77788370	77793883
SCRP324	78256282	84618946
SCRP333	77817047	77822560
```

```
N50:
SCRP249: 16,620
SCRP324: 17,037
SCRP333: 16,946

L50:
SCRP249: 1,232
SCRP324: 1,209
SCRP333: 1,210

Number of contigs > 1kb:
SCRP249: 9,214
SCRP324: 9,081
SCRP333: 8,985
```

Extra corrections for submitting to genbank performed in Genbank_corrections.md

<!-- ##Discovar

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    Organism=P.rubi
    F_Read=qc_dna/paired/$Organism/$Strain/F/*.gz
    R_Read=qc_dna/paired/$Organism/$Strain/R/*.gz
    OutDir=assembly/discovar/$Organism/$Strain
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/discovar
    qsub $ProgDir/sub_discovar.sh $F_Read $R_Read $OutDir
done
```

filter out contigs < 500 bp

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for AssemblyDir in $(ls -d assembly/discovar/P.rubi/$Strain/assembly/a.final)
    do
        echo $Strain
        echo "Filtering contigs smaller than 500bp"
        mkdir -p $AssemblyDir/filtered_contigs
        FilterDir=/home/adamst/git_repos/tools/seq_tools/assemblers/abyss
        $FilterDir/filter_abyss_contigs.py $AssemblyDir/a.lines.fasta 500 > $AssemblyDir/filtered_contigs/contigs_min_500bp.fasta
    done
done
```

### Quast

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    OutDir=$(ls -d assembly/discovar/*/$Strain/assembly/a.final/filtered_contigs)
    AssFiltered=$OutDir/contigs_min_500bp.fasta
    AssRenamed=$OutDir/contigs_min_500bp_renamed.fasta
    echo $AssFiltered
    printf '.\t.\t.\t.\n' > editfile.tab
    $ProgDir/remove_contaminants.py --inp $AssFiltered --out $AssRenamed --coord_file editfile.tab
    rm editfile.tab
done
```

### QUAST used to summarise assembly statistics


```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Strain in SCRP249 SCRP324 SCRP333
do
    for Assembly in $(ls assembly/discovar/*/$Strain/assembly/a.final/filtered_contigs/contigs_min_500bp_renamed.fasta)
    do
        Organism=P.rubi
        OutDir=assembly/discovar/$Organism/$Strain/assembly/a.final/filtered_contigs/QUAST
        qsub $ProgDir/sub_quast.sh $Assembly $OutDir
    done
done
```

** N50:
SCRP249: 16638
SCRP324: 18682
SCRP333: 16946 **

** L50:
SCRP249: 1332
SCRP324: 1185
SCRP333: 1321 **

** Number of contigs > 1kb:
SCRP249: 10353
SCRP324: 10480
SCRP333: 10260 ** -->

#Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assemblies were used to perform repeatmasking

<!-- for discovar assembly:

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/repeat_masking
for Strain in SCRP249 SCRP324 SCRP333
do
    for BestAss in $(ls assembly/discovar/*/$Strain/assembly/a.final/filtered_contigs/contigs_min_500bp_renamed.fasta)
    do
        qsub $ProgDir/rep_modeling.sh $BestAss
        qsub $ProgDir/transposonPSI.sh $BestAss
    done
done
``` -->

for SPAdes assembly:

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/repeat_masking
for Strain in SCRP249 SCRP324 SCRP333
do
    for BestAss in $(ls assembly/spades/*/$Strain/deconseq_Paen/*_500bp_filtered_renamed.fasta)
    do
        qsub $ProgDir/rep_modeling.sh $BestAss
        qsub $ProgDir/transposonPSI.sh $BestAss
    done
done   
```

For manually cleaned assemblies

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/repeat_masking
for Strain in SCRP324
do
    for BestAss in $(ls assembly/spades/*/$Strain/manual_edits/contigs_min_500bp_renamed.fasta)
    do
        qsub $ProgDir/rep_modeling.sh $BestAss
        qsub $ProgDir/transposonPSI.sh $BestAss
    done
done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised using the following commands:

```bash
for RepDir in $(ls -d repeat_masked/P.*/*/deconseq_Paen_repmask)
do
    Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
    RepMaskGff=$(ls $RepDir/"$Strain"_contigs_hardmasked.gff)
    TransPSIGff=$(ls $RepDir/"$Strain"_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    printf "$Organism\t$Strain\n"
    printf "The number of bases masked by RepeatMasker:\t"
    sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The number of bases masked by TransposonPSI:\t"
    sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The total number of masked bases are:\t"
    cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done
```

For manualy cleaned assemblies

```bash
for RepDir in $(ls -d repeat_masked/P.*/*/manual_edits_repmask)
do
    Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
    RepMaskGff=$(ls $RepDir/"$Strain"_contigs_hardmasked.gff)
    TransPSIGff=$(ls $RepDir/"$Strain"_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    printf "$Organism\t$Strain\n"
    printf "The number of bases masked by RepeatMasker:\t"
    sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The number of bases masked by TransposonPSI:\t"
    sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The total number of masked bases are:\t"
    cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done
```

```
discovar:
P.rubi	SCRP249
The number of bases masked by RepeatMasker:	31,215,924
The number of bases masked by TransposonPSI:	7,577,505
The total number of masked bases are:	33,062,283
P.rubi	SCRP324
The number of bases masked by RepeatMasker:	31,949,318
The number of bases masked by TransposonPSI:	7,607,213
The total number of masked bases are:	33,859,484
P.rubi	SCRP333
The number of bases masked by RepeatMasker:	31,597,333
The number of bases masked by TransposonPSI:	7,644,529
The total number of masked bases are:	33,426,973

SPAdes:
P.rubi	SCRP249
The number of bases masked by RepeatMasker:	23,906,929
The number of bases masked by TransposonPSI:	5,953,026
The total number of masked bases are:	25,659,011
P.rubi	SCRP324
The number of bases masked by RepeatMasker:	23,969,622
The number of bases masked by TransposonPSI:	5,940,402
The total number of masked bases are:	25,657,326
P.rubi	SCRP333
The number of bases masked by RepeatMasker:	23,126,516
The number of bases masked by TransposonPSI:	5,961,557
The total number of masked bases are:	25,001,991
```

#Merging RepeatMasker and TransposonPSI outputs

```bash
for File in $(ls -d repeat_masked/P.*/*/deconseq_Paen_repmask/*_contigs_softmasked.fa)
do
    OutDir=$(dirname $File)
    TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
    bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
    echo "$OutFile"
    echo "Number of masked bases:"
    cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
```

For manually cleaned assemblies

```bash
for File in $(ls -d repeat_masked/P.*/*/manual_edits_repmask/*_contigs_softmasked.fa)
do
    OutDir=$(dirname $File)
    TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
    bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
    echo "$OutFile"
    echo "Number of masked bases:"
    cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
```

#Gene Prediction
Gene prediction followed three steps: Pre-gene prediction - Quality of genome assemblies were assessed using BUSOC to see how many eukaryotic single copy orthologs are present, and quantify levels of duplication in the assembly. Gene model training - Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline Gene prediction - Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

##Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    Organism=P.rubi
    if [ -f repeat_masked/$Organism/$Strain/manualedits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    elif [ -f repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    elif [ -f repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    fi
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/busco
    BuscoDB=Eukaryotic
    OutDir=assembly/spades/P.rubi/$Strain/Busco
    mkdir -p $OutDir
    qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```
Of 303 BUSCOs searched:

SCRP249:
Complete and single-copy BUSCOs: 273
Complete and duplicated BUSCOs: 8
Fragmented BUSCOs: 3
Missing BUSCOs: 19

SCRP324:
Complete and single-copy BUSCOs: 274
Complete and duplicated BUSCOs: 8
Fragmented BUSCOs: 3
Missing BUSCOs: 18

SCRP333:
Complete and single-copy BUSCOs: 275
Complete and duplicated BUSCOs: 7
Fragmented BUSCOs: 3
Missing BUSCOs: 18
```

#Gene prediction
Gene prediction was performed for the P. rubi genomes. Two gene prediction approaches were used:

Gene prediction using Braker1 and Prediction of all putative ORFs in the genome using the ORF finder (atg.pl) approach.

##Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to P. rubi genomes.

RNAseq data was obtained from the phytophthora sequencing consortium - Nick Grünwald's lab at OSU

RNA-Seq data aligned using STAR

```bash
for Assembly in $(ls repeat_masked/P.rubi/*/deconseq_Paen_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
do
    Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for FileF in $(ls qc_rna/qc_rna/raw_rna/consortium/P.rubi/F/*_trim.fq.gz)
    do
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw' | wc -l)
        while [ $Jobs -gt 1 ]
        do
            sleep 1m
            printf "."
            Jobs=$(qstat | grep 'sub_sta' | grep 'qw' | wc -l)
        done
        printf "\n"
        FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1/_2/g')
        echo $FileF
        echo $FileR
        Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
        OutDir=alignment/star/$Organism/$Strain/$Sample_Name
        ProgDir=/home/adamst/git_repos/tools/seq_tools/RNAseq
        qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
    done
done
```

Repeat for manually cleaned assemblies

```bash
for Assembly in $(ls repeat_masked/P.rubi/*/manual_edits_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
do
    Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for FileF in $(ls qc_rna/qc_rna/raw_rna/consortium/P.rubi/F/*_trim.fq.gz)
    do
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw' | wc -l)
        while [ $Jobs -gt 1 ]
        do
            sleep 1m
            printf "."
            Jobs=$(qstat | grep 'sub_sta' | grep 'qw' | wc -l)
        done
        printf "\n"
        FileR=$(echo $FileF | sed 's&/F/&/R/&g'| sed 's/_1/_2/g')
        echo $FileF
        echo $FileR
        Sample_Name=$(echo $FileF | rev | cut -d '/' -f1 | rev | sed 's/_1_trim.fq.gz//g')
        OutDir=alignment/star/$Organism/$Strain/$Sample_Name
        ProgDir=/home/adamst/git_repos/tools/seq_tools/RNAseq
        qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
    done
done
```

#Braker prediction

```bash
for Assembly in $(ls repeat_masked/P.rubi/*/deconseq_Paen_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'SCRP324')
do
    Jobs=$(qstat | grep 'braker' | grep -w 'r' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'braker' | grep -w 'r' | wc -l)
    done
    printf "\n"
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    mkdir -p alignment/star/$Organism/$Strain/concatenated
    samtools merge -f alignment/star/$Organism/$Strain/concatenated/concatenated.bam \
    alignment/star/$Organism/$Strain/4671V8/star_aligmentAligned.sortedByCoord.out.bam \
    alignment/star/$Organism/$Strain/Pr4671PB/star_aligmentAligned.sortedByCoord.out.bam
    OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    AcceptedHits=alignment/star/$Organism/$Strain/concatenated/concatenated.bam
    GeneModelName="$Organism"_"$Strain"_braker
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

Repeat for manually cleaned assemblies

```bash
for Assembly in $(ls repeat_masked/P.rubi/*/manual_edits_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
do
    Jobs=$(qstat | grep 'braker' | grep -w 'r' | wc -l)
    while [ $Jobs -gt 1 ]
    do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'braker' | grep -w 'r' | wc -l)
    done
    printf "\n"
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    mkdir -p alignment/star/$Organism/$Strain/concatenated
    samtools merge -f alignment/star/$Organism/$Strain/concatenated/concatenated.bam \
    alignment/star/$Organism/$Strain/4671V8/star_aligmentAligned.sortedByCoord.out.bam \
    alignment/star/$Organism/$Strain/Pr4671PB/star_aligmentAligned.sortedByCoord.out.bam
    OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    AcceptedHits=alignment/star/$Organism/$Strain/concatenated/concatenated.bam
    GeneModelName="$Organism"_"$Strain"_braker
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

#Supplementing Braker gene models with CodingQuarry genes

Additional genes were added to Braker gene predictions, using CodingQuarry in pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/deconseq_Paen_repmask/*_contigs_unmasked.fa | grep -v 'SCRP324')
do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/star/cufflinks/$Organism/$Strain/concatenated
    mkdir -p $OutDir
    AcceptedHits=alignment/star/$Organism/$Strain/concatenated/concatenated.bam
    ProgDir=/home/adamst/git_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Repeat for manually cleaned assemblies

```bash
for Assembly in $(ls repeat_masked/*/*/manual_edits_repmask/*_contigs_unmasked.fa)
do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/star/cufflinks/$Organism/$Strain/concatenated
    mkdir -p $OutDir
    AcceptedHits=alignment/star/$Organism/$Strain/concatenated/concatenated.bam
    ProgDir=/home/adamst/git_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Secondly, genes were predicted using CodingQuarry:

```bash
for Assembly in $(ls repeat_masked/*/*/deconseq_Paen_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'SCRP324')
do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquarry/$Organism/$Strain
    CufflinksGTF=gene_pred/star/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```

Repeated on manually cleaned assemblies

```bash
for Assembly in $(ls repeat_masked/*/*/manual_edits_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquarry/$Organism/$Strain
    CufflinksGTF=gene_pred/star/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```

Then, additional transcripts were added to Braker1 gene models, when CodingQuarry genes were predicted in regions of the genome, not containing Braker1 gene models:

```bash
for BrakerGff in $(ls gene_pred/braker/P.rubi/*_braker/*/augustus.gff3)
do
    Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
    Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    if [ -f repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    elif [ -f repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    else [ -f repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
        Assembly=$(ls repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    fi
    CodingQuarryGff=gene_pred/codingquarry/$Organism/$Strain/out/PredictedPass.gff3
    PGNGff=gene_pred/codingquarry/$Organism/$Strain/out/PGN_predictedPass.gff3
    AddDir=gene_pred/codingquarry/$Organism/$Strain/additional
    FinalDir=gene_pred/final/$Organism/$Strain/final
    AddGenesList=$AddDir/additional_genes.txt
    AddGenesGff=$AddDir/additional_genes.gff
    FinalGff=$AddDir/combined_genes.gff
    mkdir -p $AddDir
    mkdir -p $FinalDir

    bedtools intersect -v -a $CodingQuarryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
    bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuarryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
    $ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
    # -
    # This section is edited
    $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuarry_unspliced.gff3
    $ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuarry_unspliced.gff3 > $FinalDir/final_genes_CodingQuarry.gff3
    # -
    $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuarry.gff3 $FinalDir/final_genes_CodingQuarry
    cp $BrakerGff $FinalDir/final_genes_Braker.gff3
    $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
    cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuarry.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
    cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuarry.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
    cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuarry.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
    cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuarry.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta


    GffBraker=$FinalDir/final_genes_Braker.gff3
    GffQuarry=$FinalDir/final_genes_CodingQuarry.gff3
    GffAppended=$FinalDir/final_genes_appended.gff3
    cat $GffBraker $GffQuarry > $GffAppended
done
```

Then, duplicated transcripts were removed

```bash
for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3)
do
    Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    FinalDir=gene_pred/final/$Organism/$Strain/final
    GffFiltered=$FinalDir/filtered_duplicates.gff
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
    $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
    GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
    LogFile=$FinalDir/final_genes_appended_renamed.log
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
    $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
    rm $GffFiltered
    if [ -f repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    elif [ -f repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    elif [ -f repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
        echo $Assembly
    fi
    $ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed
    # The proteins fasta file contains * instead of Xs for stop codons, these should
    # be changed
    sed -i 's/\*/X/g' gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta
done
```

The final number of genes per isolate was observed using:

```bash
for DirPath in $(ls -d gene_pred/final/P.*/*/final)
do
    echo $DirPath
    echo Braker:
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l
    echo CodingQuarry:
    cat $DirPath/final_genes_CodingQuarry.pep.fasta | grep '>' | wc -l
    echo Total:
    cat $DirPath/final_genes_appended_renamed.pep.fasta | grep '>' | wc -l
    echo ""
done
```

```
discovar:
P.rubi - SCRP249
Braker:
33629
CodingQuarry:
2490
Total:
36119

P.rubi - SCRP324
Braker:
39285
CodingQuarry:
2909
Total:
42194

P.rubi - SCRP333
Braker:
34276
CodingQuarry:
3236
Total:
37512


SPAdes:
P.rubi - SCRP249
Braker:
29,583
CodingQuarry:
2,225
Total:
31,808

P.rubi - SCRP324
Braker:
20,707
CodingQuarry:
10,524
Total:
31,231

P.rubi - SCRP333
Braker:
29,696
CodingQuarry:
3,250
Total:
32,946
```

##Predicted gene set assessed using BUSCO to assess completeness

```bash
for Transcriptome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gene.fasta )
do
    Strain=$(echo $Transcriptome| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Transcriptome | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/busco
    BuscoDB=Eukaryotic
    OutDir=gene_pred/busco/$Organism/$Strain/
    qsub $ProgDir/sub_busco3.sh $Transcriptome $BuscoDB $OutDir
done
```

```
SCRP249
Complete and single copy genes: 269
Complete and duplicated genes: 9
Fragmented genes: 8
Missing genes: 17

SCRP324
Complete and single copy genes: 265
Complete and duplicated genes: 11
Fragmented genes: 10
Missing genes: 17

SCRP333
Complete and single copy genes: 270
Complete and duplicated genes: 10
Fragmented genes: 7
Missing genes: 16
```

Changes with respect to genome sequence

```
SCRP249
Complete and single copy genes: -4
Complete and duplicated genes: +1
Fragmented genes: +5
Missing genes: -2

SCRP324
Complete and single copy genes: -5
Complete and duplicated genes: +2
Fragmented genes: +4
Missing genes: -1

SCRP333
Complete and single copy genes: -5
Complete and duplicated genes: +3
Fragmented genes: +4
Missing genes: -2
```

#Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the path_pipe.sh pipeline. This pipeline also identifies open reading frames containing Signal peptide sequences and RxLRs. This pipeline was run with the following commands:


```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    Organism=P.rubi
    if [ -f repeat_masked/$Organism/$Strain/manual_edits_repmask/*_unmasked.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/manual_edits_repmask/*_unmasked.fa)
        echo $Assembly
    elif [ -f repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_unmasked.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_unmasked.fa)
        echo $Assembly
    elif [ -f repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_unmasked.fa ]
    then
        Assembly=$(ls repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_unmasked.fa)
        echo $Assembly
    fi
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    qsub $ProgDir/run_ORF_finder.sh $Assembly
done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were corrected using the following commands:

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
for OrfGff in $(ls gene_pred/ORF_finder/P.*/*/*_ORF.gff | grep -v 'atg')
do
    echo "$OrfGff"
    OrfGffMod=$(echo $OrfGff | sed 's/.gff/.gff3/g')
    $ProgDir/gff_corrector.pl $OrfGff > $OrfGffMod
done
```

The final number of genes per isolate were determined using

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for DirPath in $(ls -d gene_pred/ORF_finder/P.rubi/$Strain)
    do
        echo $DirPath
        cat $DirPath/"$Strain".aa_cat.fa | grep '>' | wc -l
        echo ""
    done
done
```

```
discovar:

SCRP249
716,788

SCRP324
795,934

SCRP333
720,807

SPAdes:

SCRP249
645,852

SCRP324
647,415

SCRP333
646,159
```

#Genomic analysis

##Effector genes

Putative effector genes were identified within Augustus gene models and ORF fragments using a number of approaches:

A) From Augustus gene models - Signal peptides
B) From Augustus gene models - Regex evidence of RxLRs
C) From Augustus gene models - Hmm evidence of RxLR effectors
D) From Augustus gene models - Hmm evidence of CRN effectors
E) From Augustus gene models - ApoplastP prediction of apoplastic effectors
F) From ORF fragments - Signal peptides
G) From ORF fragments - Regex evidence of RxLRs
H) From ORF fragments - Hmm evidence of RxLR effectors
I) From ORF fragments - Hmm evidence of CRN effectors
J) From ORF fragments - ApoplastP prediction of apoplastic effectors

###A) From Augustus gene models - Signal peptide

Required programs:

SigP
biopython

####A.1) Signal peptide prediction using SignalP 2.0, SignalP 3.0 & SignalP4.1

Proteins that were predicted to contain signal peptides were identified using the following commands:

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for Proteome in $(ls gene_pred/final/*/$Strain/final/final_genes_appended_renamed.pep.fasta)
    do
        SplitfileDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
        ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
        Organism=P.rubi
        SplitDir=gene_pred/final_split/$Organism/$Strain
        mkdir -p $SplitDir
        BaseName="$Organism""_$Strain"
        $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
        for File in $(ls $SplitDir/*_"$Strain"_*)
        do
            Jobs=$(qstat | grep 'pred_sigP' | wc -l)
            while [ $Jobs -gt 20 ]
            do
                sleep 1
                printf "."
                Jobs=$(qstat | grep 'pred_sigP' | wc -l)
            done  
            printf "\n"
            echo $File
            qsub $ProgDir/pred_sigP_modified.sh $File
            qsub $ProgDir/pred_sigP_modified.sh $File signalp-3.0
            qsub $ProgDir/pred_sigP_modified.sh $File signalp-4.1
        done
    done
done
```

The batch files of predicted secreted proteins needed to be combined into a single file for each strain. This was done with the following commands:

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for SplitDir in $(ls -d gene_pred/final_split/P.*/$Strain)
    do
        Organism=P.rubi
        echo "$Organism - $Strain"
        for SigpDir in $(ls -d gene_pred/final_sig* | cut -f2 -d'/')
        do
            InStringAA=''
            InStringNeg=''
            InStringTab=''
            InStringTxt=''
            for GRP in $(ls -l $SplitDir/*_"$Strain"_*.fa | rev | cut -d '_' -f1 | rev | sort -n)
            do  
                InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_$GRP""_sp.aa"
                InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_$GRP""_sp_neg.aa"
                InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_$GRP""_sp.tab"
                InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_$GRP""_sp.txt"
            done
            cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
            cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
            tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
            cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
        done
    done
done
```

####A.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for Proteome in $(ls gene_pred/final/*/$Strain/final/final_genes_appended_renamed.pep.fasta)
    do
        Organism=P.rubi
        echo "$Organism - $Strain"
        OutDir=analysis/phobius_CQ/$Organism/$Strain
        mkdir -p $OutDir
        phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
        ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
        $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
    done
done
```

Secreted proteins from different sources were combined into a single file:

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for Proteome in $(ls gene_pred/final/*/$Strain/final/final_genes_appended_renamed.pep.fasta)
    do
        Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
        echo "$Organism - $Strain"
        OutDir=gene_pred/combined_sigP_CQ/$Organism/$Strain
        mkdir -p $OutDir
        echo "The following number of sequences were predicted as secreted:"
        cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius_CQ/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
        cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | wc -l
        echo "This represented the following number of unique genes:"
        cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius_CQ/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
        ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
        $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
        cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
    done
done
```

```
discovar:

P.rubi - SCRP249
The following number of sequences were predicted as secreted:
10,696
This represented the following number of unique genes:
3,748
P.rubi - SCRP324
The following number of sequences were predicted as secreted:
12,476
This represented the following number of unique genes:
4,434
P.rubi - SCRP333
The following number of sequences were predicted as secreted:
10,397
This represented the following number of unique genes:
3,714

SPAdes:

P.rubi - SCRP249
The following number of sequences were predicted as secreted:
11,040
This represented the following number of unique genes:
3,697
P.rubi - SCRP324
The following number of sequences were predicted as secreted:
10,995
This represented the following number of unique genes:
3,706
P.rubi - SCRP333
The following number of sequences were predicted as secreted:
11,358
This represented the following number of unique genes:
3,832
```

####B) The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP_CQ/*/*/*_all_secreted.fa)
do
    Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev)
    Proteome=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta)
    Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
    OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain"
    mkdir -p $OutDir
    printf "\nstrain: $Strain\tspecies: $Organism\n"
    printf "the total number of SigP gene is:\t"
    cat $Secretome | grep '>' | wc -l
    printf "the number of unique SigP gene is:\t"
    cat $Secretome | grep '>' | cut -f1 | tr -d ' '| sort | uniq | wc -l
    printf "the number of SigP-RxLR genes are:\t"
    ProgDir=/home/adamst/git_repos/tools/pathogen/RxLR_effectors
    $ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_all_secreted_RxLR_regex.fa
    cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_regex.txt
    cat $OutDir/"$Strain"_RxLR_regex.txt | wc -l
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa
    printf "the number of SigP-RxLR-EER genes are:\t"
    cat $OutDir/"$Strain"_all_secreted_RxLR_regex.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $OutDir/"$Strain"_RxLR_EER_regex.txt
    cat $OutDir/"$Strain"_RxLR_EER_regex.txt | wc -l
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.fa
    printf "\n"
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_regex.txt
    sed -i -r 's/\.t.*//' $OutDir/"$Strain"_RxLR_EER_regex.txt
    cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_regex.txt> $OutDir/"$Strain"_RxLR_regex.gff3
    cat $Gff | grep -w -f $OutDir/"$Strain"_RxLR_EER_regex.txt > $OutDir/"$Strain"_RxLR_EER_regex.gff3
    echo "$Strain complete"
done
```

```
discovar:

strain: SCRP249	species: P.rubi
the total number of SigP gene is:	10,696
the number of unique SigP gene is:	3,748
the number of SigP-RxLR genes are:	379
the number of SigP-RxLR-EER genes are:	200

strain: SCRP324	species: P.rubi
the total number of SigP gene is:	12,476
the number of unique SigP gene is:	4,434
the number of SigP-RxLR genes are:	377
the number of SigP-RxLR-EER genes are:	192

strain: SCRP333	species: P.rubi
the total number of SigP gene is:	10,397
the number of unique SigP gene is:	3,714
the number of SigP-RxLR genes are:	354
the number of SigP-RxLR-EER genes are:	179

SPAdes:

strain: SCRP249	species: P.rubi
the total number of SigP gene is:	11,039
the number of unique SigP gene is:	3,696
the number of SigP-RxLR genes are:	363
the number of SigP-RxLR-EER genes are:	188

strain: SCRP324	species: P.rubi
the total number of SigP gene is:	10,995
the number of unique SigP gene is:	3,706
the number of SigP-RxLR genes are:	348
the number of SigP-RxLR-EER genes are:	174

strain: SCRP333	species: P.rubi
the total number of SigP gene is:	11,351
the number of unique SigP gene is:	3,832
the number of SigP-RxLR genes are:	359
the number of SigP-RxLR-EER genes are:	188
```

####C) From Secreted gene models - HMM evidence of RxLR effectors

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for Proteome in $(ls gene_pred/final/*/$Strain/final/final_genes_appended_renamed.pep.fasta)
    do
        ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer
        HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
        Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
        OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
        mkdir -p $OutDir
        HmmResults="$Strain"_RxLR_hmmer.txt
        hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
        echo "$Organism $Strain"
        cat $OutDir/$HmmResults | grep 'Initial search space'
        cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
        HmmFasta="$Strain"_RxLR_hmmer.fa
        $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
        Headers="$Strain"_RxLR_hmmer_headers.txt
        cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' | sort | uniq > $OutDir/$Headers
        Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
        cat $Gff | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
    done
    echo "$Strain complete"
done
```

```
discovar:

P.rubi SCRP249
Initial search space (Z):             36,943  [actual number of targets]
Domain search space  (domZ):             202  [number of targets reported over threshold]
P.rubi SCRP324
Initial search space (Z):             42,604  [actual number of targets]
Domain search space  (domZ):             207  [number of targets reported over threshold]
P.rubi SCRP333
Initial search space (Z):             36,843  [actual number of targets]
Domain search space  (domZ):             189  [number of targets reported over threshold]

SPAdes:

P.rubi SCRP249
Initial search space (Z):             31,808  [actual number of targets]
Domain search space  (domZ):             197  [number of targets reported over threshold]
P.rubi SCRP324
Initial search space (Z):             31,231  [actual number of targets]
Domain search space  (domZ):             194  [number of targets reported over threshold]
P.rubi SCRP333
Initial search space (Z):             32,946  [actual number of targets]
Domain search space  (domZ):             207  [number of targets reported over threshold]
```

####Combining RxLRs from Regex and hmm searches

```bash
#Without EER
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_regex.txt)
do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
    Proteome=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta)
    HmmRxLR=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer_headers.txt
    echo "$Organism - $Strain"
    echo "Number of RxLRs identified by Regex:"
    cat $RegexRxLR | sort | uniq | wc -l
    echo "Number of RxLRs identified by Hmm:"
    cat $HmmRxLR | sort | uniq | wc -l
    echo "Number of RxLRs in combined dataset:"
    cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
    # echo "Number of RxLRs in both datasets:"
    # cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
    echo ""
    # echo "Extracting RxLRs from datasets"
    OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
    mkdir -p $OutDir
    cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_RxLR_headers.txt
    Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
    cat $Gff | grep -w -f $OutDir/"$Strain"_total_RxLR_headers.txt > $OutDir/"$Strain"_total_RxLR.gff
    echo "Number of genes in the extracted gff file:"
    cat $OutDir/"$Strain"_total_RxLR.gff | grep -w 'gene' | wc -l
    echo "$Strain complete without EER"
done

#With EER
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_RxLR_EER_regex.txt)
do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
    Proteome=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta)
    HmmRxLR=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/*_RxLR_hmmer_headers.txt
    echo "$Organism - $Strain"
    echo "Number of RxLRs identified by Regex:"
    cat $RegexRxLR | sort | uniq | wc -l
    echo "Number of RxLRs identified by Hmm:"
    cat $HmmRxLR | sort | uniq | wc -l
    echo "Number of RxLRs in combined dataset:"
    cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
    # echo "Number of RxLRs in both datasets:"
    # cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
    echo ""
    # echo "Extracting RxLRs from datasets"
    OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
    mkdir -p $OutDir
    cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_RxLR_EER_headers.txt
    Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
    cat $Gff | grep -w -f $OutDir/"$Strain"_total_RxLR_EER_headers.txt > $OutDir/"$Strain"_total_RxLR_EER.gff
    echo "Number of genes in the extracted gff file:"
    cat $OutDir/"$Strain"_total_RxLR_EER.gff | grep -w 'gene' | wc -l
    echo "$Strain complete with EER"
done
```

```
With EER:
discovar:

P.rubi - SCRP249
Number of RxLRs identified by Regex:
200
Number of RxLRs identified by Hmm:
202
Number of RxLRs in combined dataset:
245
Number of genes in the extracted gff file:
245

P.rubi - SCRP324
Number of RxLRs identified by Regex:
192
Number of RxLRs identified by Hmm:
207
Number of RxLRs in combined dataset:
246
Number of genes in the extracted gff file:
246

P.rubi - SCRP333
Number of RxLRs identified by Regex:
179
Number of RxLRs identified by Hmm:
189
Number of RxLRs in combined dataset:
229
Number of genes in the extracted gff file:
229

SPAdes:

With EER:
P.rubi - SCRP249
Number of RxLRs identified by Regex:
188
Number of RxLRs identified by Hmm:
197
Number of RxLRs in combined dataset:
238
Number of genes in the extracted gff file:
238

P.rubi - SCRP324
Number of RxLRs identified by Regex:
174
Number of RxLRs identified by Hmm:
193
Number of RxLRs in combined dataset:
229
Number of genes in the extracted gff file:
229

P.rubi - SCRP333
Number of RxLRs identified by Regex:
188
Number of RxLRs identified by Hmm:
207
Number of RxLRs in combined dataset:
244
Number of genes in the extracted gff file:
244

Without EER:
P.rubi - SCRP249
Number of RxLRs identified by Regex:
362
Number of RxLRs identified by Hmm:
197
Number of RxLRs in combined dataset:
407
Number of genes in the extracted gff file:
407
P.rubi - SCRP324
Number of RxLRs identified by Regex:
346
Number of RxLRs identified by Hmm:
193
Number of RxLRs in combined dataset:
393
Number of genes in the extracted gff file:
393
P.rubi - SCRP333
Number of RxLRs identified by Regex:
359
Number of RxLRs identified by Hmm:
207
Number of RxLRs in combined dataset:
410
Number of genes in the extracted gff file:
410
```

####D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers in Augustus gene models. This was done with the following commands:

```bash
HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
for Proteome in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.pep.fasta)
do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
    mkdir -p $OutDir
    echo "$Organism - $Strain"
    # Run hmm searches LFLAK domains
    CrinklerProts_LFLAK=$OutDir/"$Strain"_pub_CRN_LFLAK_hmm.txt
    hmmsearch -T0 $LFLAK_hmm $Proteome > $CrinklerProts_LFLAK
    cat $CrinklerProts_LFLAK | grep 'Initial search space'
    cat $CrinklerProts_LFLAK | grep 'number of targets reported over threshold'
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer
    $ProgDir/hmmer2fasta.pl $CrinklerProts_LFLAK $Proteome > $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa
    # Run hmm searches DWL domains
    CrinklerProts_DWL=$OutDir/"$Strain"_pub_CRN_DWL_hmm.txt
    hmmsearch -T0 $DWL_hmm $Proteome > $CrinklerProts_DWL
    cat $CrinklerProts_DWL | grep 'Initial search space'
    cat $CrinklerProts_DWL | grep 'number of targets reported over threshold'
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer
    $ProgDir/hmmer2fasta.pl $CrinklerProts_DWL $Proteome > $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa
    # Identify the genes detected in both models
    cat $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt
    echo "Total number of CRNs from both models"
    cat $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt | wc -l
    echo "$Strain done"
done
```

```
discovar

P.rubi - SCRP249
Initial search space (Z):              36943  [actual number of targets]
Domain search space  (domZ):             161  [number of targets reported over threshold]
Initial search space (Z):              36943  [actual number of targets]
Domain search space  (domZ):             152  [number of targets reported over threshold]
Identify the genes detected in both models:     135

P.rubi - SCRP324
Initial search space (Z):              42604  [actual number of targets]
Domain search space  (domZ):             159  [number of targets reported over threshold]
Initial search space (Z):              42604  [actual number of targets]
Domain search space  (domZ):             151  [number of targets reported over threshold]
Identify the genes detected in both models:     135

P.rubi - SCRP333
Initial search space (Z):              36843  [actual number of targets]
Domain search space  (domZ):             145  [number of targets reported over threshold]
Initial search space (Z):              36843  [actual number of targets]
Domain search space  (domZ):             136  [number of targets reported over threshold]
Identify the genes detected in both models:     123

SPAdes:

P.rubi - SCRP249
Initial search space (Z):             31,808  [actual number of targets]
Domain search space  (domZ):             156  [number of targets reported over threshold]
Initial search space (Z):             31,808  [actual number of targets]
Domain search space  (domZ):             139  [number of targets reported over threshold]
Identify the genes detected in both models:     128

P.rubi - SCRP324
Initial search space (Z):             31,231  [actual number of targets]
Domain search space  (domZ):              90  [number of targets reported over threshold]
Initial search space (Z):             31,231  [actual number of targets]
Domain search space  (domZ):             100  [number of targets reported over threshold]
Identify the genes detected in both models:     69

P.rubi - SCRP333
Initial search space (Z):             32,946  [actual number of targets]
Domain search space  (domZ):             154  [number of targets reported over threshold]
Initial search space (Z):             32,946  [actual number of targets]
Domain search space  (domZ):             142  [number of targets reported over threshold]
Identify the genes detected in both models:     128
```

Extract gff annotations for Crinklers:

```bash
for CRNlist in $(ls analysis/CRN_effectors/hmmer_CRN/*/*/*_pub_CRN_LFLAK_DWL.txt)
do
    Strain=$(echo $CRNlist | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $CRNlist | rev | cut -f3 -d '/' | rev)
    OutName=$(echo $CRNlist | sed 's/.txt/.gff/g')
    echo "$Organism - $Strain"
    Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
    cat $CRNlist | sed -r 's/\.t.$//g' > tmp.txt
    cat $Gff | grep -w -f tmp.txt > $OutName
    rm tmp.txt
done
```

####E) From Augustus gene models - ApoplastP prediction of apoplastic effectors

```bash
for Secretome in $(ls gene_pred/combined_sigP_CQ/*/*/*_all_secreted.fa)
do
    Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    BaseName="$Organism"_"$Strain"_ApoplastP_raw
    OutDir=analysis/ApoplastP/$Organism/$Strain
    ProgDir=/home/adamst/git_repos/tools/seq_tools/apoplastic_effectors
    qsub $ProgDir/pred_apoplastP.sh $Secretome $BaseName $OutDir
done
```

TEMPORARY COMMANDS TO RUN ON THE HEAD NODE

```bash
for Secretome in $(ls gene_pred/combined_sigP_CQ/*/*/*_all_secreted.fa)
do
    Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    BaseName="$Organism"_"$Strain"_ApoplastP
    OutDir=analysis/ApoplastP/$Organism/$Strain
    mkdir -p $OutDir
    ApoplastP.py -o "$OutDir"/"$BaseName".txt -A "$OutDir"/"$BaseName"_raw.fa -i $Secretome
done
```

The number of proteins predicted as being apoplastic effectors were summarised using the following commands

```bash
for File in $(ls analysis/ApoplastP/*/*/*_ApoplastP_raw.fa)
do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo $File | sed 's/_ApoplastP_raw.fa/_ApoplastP_headers.txt/g')
    Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
    Fasta=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta)
    echo "Creating Headers file"
    cat $File | grep '>' | sed 's/>//g' | cut -f1 | tr -d '>' | tr -d ' ' | sort -g | uniq > $Headers
    echo "The number of genes predicted as Apoplastic effectors is:"
    cat $Headers | wc -l
    echo "Creating GFF3 file"
    OutName=$(echo $File | sed 's/_raw.fa/.gff/g')
    Test=analysis/ApoplastP/$Organism/$Strain/tmp.txt
    cat $Headers | sed 's/\..*$//g' > $Test
    cat $Gff | grep -w -f $Test > $OutName
    rm $Test
    echo "Number of genes extracted into GFF3 file is:"
    cat $OutName | grep -w 'gene' | wc -l
    echo "Creating parsed fasta file"
    Parsed_Fasta=$(echo $File | sed 's/_ApoplastP_raw.fa/_ApoplastP.fa/g')
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Fasta --headers $Headers > $Parsed_Fasta
    echo "Number of genes extracted to Fasta file is:"
    cat $Parsed_Fasta | grep '>' | wc -l
done
```

```
P.rubi - SCRP249
The number of genes predicted as Apoplastic effectors is:
1,017
Number of genes extracted into GFF file is:
1,017
Number of genes extracted to Fasta file is:
1,017

P.rubi - SCRP324
The number of genes predicted as Apoplastic effectors is:
1,069
Number of genes extracted into GFF file is:
1,065
Number of genes extracted to Fasta file is:
1,069

P.rubi - SCRP333
The number of genes predicted as Apoplastic effectors is:
1,090
Number of genes extracted into GFF file is:
1,090
Number of genes extracted to Fasta file is:
1,090
```

####F) From ORF gene models - Signal peptides

Required programs:

SigP
Phobius
biopython

#####F.1) Prediction using SignalP

Proteins that were predicted to contain signal peptides were identified using the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa)
do
    SplitfileDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    SplitDir=gene_pred/ORF_split/$Organism/$Strain
    mkdir -p $SplitDir
    BaseName="$Organism""_$Strain"_ORF_preds
    $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
    for File in $(ls $SplitDir/*_ORF_preds_*)
    do
        Jobs=$(qstat | grep 'pred_sigP' | wc -l)
        while [ $Jobs -gt 20 ]
        do
            sleep 1
            printf "."
            Jobs=$(qstat | grep 'pred_sigP' | wc -l)
        done        
        printf "\n"
        echo $File
        qsub $ProgDir/pred_sigP_modified.sh $File
        qsub $ProgDir/pred_sigP_modified.sh $File signalp-3.0
        qsub $ProgDir/pred_sigP_modified.sh $File signalp-4.1
    done
done
```

The batch files of predicted secreted proteins needed to be combined into a single file for each strain. This was done with the following commands:

```bash
for SplitDir in $(ls -d gene_pred/ORF_split/P.rubi/*)
do
    Strain=$(echo $SplitDir | cut -d '/' -f4)
    Organism=$(echo $SplitDir | cut -d '/' -f3)
    echo "$Organism - $Strain"
    for SigpDir in $(ls -d gene_pred/ORF_sig* | cut -f2 -d'/')
    do
        InStringAA=''
        InStringNeg=''
        InStringTab=''
        InStringTxt=''
        for GRP in $(ls -l $SplitDir/*_ORF_*.fa | rev | cut -d '_' -f1 | rev | sort -n)
        do  
            InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.aa"
            InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp_neg.aa"
            InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.tab"
            InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.txt"
        done
        cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
        cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
        tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
        cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
    done
done
```

#####F.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/phobius_ORF/$Organism/$Strain
    mkdir -p $OutDir
    phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
    $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
done
```

Because of the way ORF_finder predicts proteins, phobius predictions cannot be used downstream as there is no way to remove overlapping features.

Secreted proteins from different sources were combined into a single file:

```bash
for Proteome in $(ls gene_pred/ORF_finder/P.*/*/*.aa_cat.fa)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/combined_sigP_ORF/$Organism/$Strain
    mkdir -p $OutDir
    echo "The following number of sequences were predicted as secreted:"
    # cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius_ORF/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
    cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa > $OutDir/"$Strain"_all_secreted.fa
    cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | tr -d '>' | tr -d ' ' | sed "s/HMM_score\t/HMM_score=\t/g" > $OutDir/"$Strain"_all_secreted_headers.txt
    cat $OutDir/"$Strain"_all_secreted_headers.txt | wc -l
    echo "This represented the following number of unique genes:"
    # cat gene_pred/final_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    cat gene_pred/ORF_sig*/$Organism/$Strain/*_aug_sp.aa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
    cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
done
```

<!-- For spades assembly of SCRP333, the production of all_secreted_headers.txt required the removal of one tr step shown below:

```bash
for Proteome in $(ls gene_pred/ORF_finder/spades/P.*/SCRP333/*.aa_cat.fa)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Assembler=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Assembler - $Organism - $Strain"
    OutDir=gene_pred/combined_sigP_ORF/$Assembler/$Organism/$Strain
    mkdir -p $OutDir
    echo "The following number of sequences were predicted as secreted:"
    # cat gene_pred/"$Assembler"_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Assembler/$Organism/$Strain/"$Strain"_phobius_ORF.fa > $OutDir/"$Strain"_all_secreted.fa
    cat gene_pred/"$Assembler"_sig*/$Organism/$Strain/*_aug_sp.aa > $OutDir/"$Strain"_all_secreted.fa
    cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | tr -d '>' | sed "s/HMM_score\t/HMM_score=\t/g" > $OutDir/"$Strain"_all_secreted_headers.txt
    cat $OutDir/"$Strain"_all_secreted_headers.txt | wc -l
    echo "This represented the following number of unique genes:"
    # cat gene_pred/"$Assembler"_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Assembler/$Organism/$Strain/"$Strain"_phobius_ORF.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    cat gene_pred/"$Assembler"_sig*/$Organism/$Strain/*_aug_sp.aa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
    cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
done
``` -->

```
Discovar:
P.rubi - SCRP249
The following number of sequences were predicted as secreted:
81,145
This represented the following number of unique genes:
38,987
P.rubi - SCRP324
The following number of sequences were predicted as secreted:
89,423
This represented the following number of unique genes:
42,890
P.rubi - SCRP333
The following number of sequences were predicted as secreted:
81,434
This represented the following number of unique genes:
39,145

SPAdes:
P.rubi - SCRP249
The following number of sequences were predicted as secreted:
82,268
This represented the following number of unique genes:
40,893
P.rubi - SCRP324
The following number of sequences were predicted as secreted:
82,082
This represented the following number of unique genes:
40,843
P.rubi - SCRP333
The following number of sequences were predicted as secreted:
82,132
This represented the following number of unique genes:
40,825
```

####G) From ORF fragments - The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

Names of ORFs containing signal peptides were extracted from fasta files. This included information on the position and hmm score of RxLRs.

```bash
for FastaFile in $(ls gene_pred/combined_sigP_ORF/P.rubi/*/*_all_secreted.fa)
do
    Strain=$(echo $FastaFile | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $FastaFile | rev | cut -f3 -d '/' | rev)
    SigP_headers=gene_pred/ORF_sigP/$Organism/$Strain/"$Strain"_ORF_sp_names.txt
    cat $FastaFile | grep '>' | sed -r 's/>//g' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | sed 's/--//g' > $SigP_headers
    echo "$Strain done"
done
```

<!-- Due to errors in spades SCRP333 fasta file, script was modified to look the same as those produced above

```bash
for FastaFile in $(ls gene_pred/combined_sigP_ORF/spades/*/SCRP333/SCRP333_all_secreted.fa)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Assembler=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SigP_headers=gene_pred/combined_sigP_ORF/$Assembler/$Organism/$Strain/"$Strain"_all_secreted_headers.txt
    cat $FastaFile | grep '>' | sed -r 's/>//g' | sed "s/HMM_score\t/HMM_score=\t/g" | sed -r 's/\s+/\t/g' > $SigP_headers
done
``` -->

Due to the nature of predicting ORFs, some features overlapped with one another. A single ORF was selected from each set of overlapped ORFs. This was was selected on the basis of its SignalP Hmm score. Biopython was used to identify overlaps and identify the ORF with the best signalP score.

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
    for ORF_Gff in $(ls gene_pred/ORF_finder/P.rubi/$Strain/"$Strain"_ORF.gff3)
    do
        Organism=$(echo $ORF_Gff | rev |  cut -d '/' -f3 | rev)
        OutDir=$(ls -d gene_pred/combined_sigP_ORF/$Organism/$Strain)
        echo "$Organism - $Strain"
        SigP_fasta=$(ls $OutDir/"$Strain"_all_secreted.fa)
        SigP_headers=$(ls $OutDir/"$Strain"_all_secreted_headers.txt)
        ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)

        SigP_Gff=$OutDir/"$Strain"_all_secreted_unmerged.gff
        SigP_Merged_Gff=$OutDir/"$Strain"_all_secreted_merged.gff
        SigP_Merged_txt=$OutDir/"$Strain"_all_secreted_merged.txt
        SigP_Merged_AA=$OutDir/"$Strain"_all_secreted_merged.aa

        ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
        $ProgDir/extract_gff_for_sigP_hits.pl $SigP_headers $ORF_Gff SigP Name > $SigP_Gff
        echo "Extracting Gff done"
        ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
        $ProgDir/make_gff_database.py --inp $SigP_Gff --db sigP_ORF.db
        echo "Gff database made"
        ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
        $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF.db --id sigP_ORF --out sigP_ORF_merged.db --gff > $SigP_Merged_Gff
        echo "Merging complete"
        cat $SigP_Merged_Gff | grep 'transcript' | rev | cut -f1 -d'=' | rev > $SigP_Merged_txt
        $ProgDir/extract_from_fasta.py --fasta $SigP_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
        # $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $SigP_Merged_txt > $SigP_Merged_AA
    done
    echo "$Strain complete"
done
```

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa)
do
    ProgDir=/home/adamst/git_repos/tools/pathogen/RxLR_effectors
    Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev)
    OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain"
    mkdir -p $OutDir
    printf "\nstrain: $Strain\tspecies: $Organism\n"
    printf "the number of SigP gene is:\t"
    cat $Secretome | grep '>' | wc -l
    printf "the number of SigP-RxLR genes are:\t"
    $ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa
    cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt
    cat $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
    printf "the number of SigP-RxLR-EER genes are:\t"
    cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' '> $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt
    cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
    ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
    # $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt  $SigP_Gff   RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
    $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt  $SigP_Gff   RxLR_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
    RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_regex_merged.gff
    RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_regex_merged.txt
    RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_regex_merged.aa
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff --db sigP_ORF_RxLR.db
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR.db --id sigP_ORF_RxLR --out sigP_ORF_RxLR_merged.db --gff > $RxLR_Merged_Gff
    cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
    $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
    printf "Merged RxLR regex proteins:\t"
    cat $RxLR_Merged_AA | grep '>' | wc -l
    printf "\n"
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt  $SigP_Gff   RxLR_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff
    RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.gff
    RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.txt
    RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.aa
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.gff --db sigP_ORF_RxLR.db
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR.db --id sigP_ORF_RxLR --out sigP_ORF_RxLR_merged.db --gff > $RxLR_Merged_Gff
    cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
    $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
    printf "Merged RxLR-EER regex proteins:\t"
    cat $RxLR_Merged_AA | grep '>' | wc -l
    printf "\n"
    echo "$Strain done"
done
```

```
Discovar:
strain: SCRP249 species: P.rubi
the number of SigP genes is:    81145
the number of SigP-RxLR genes are:      2481
the number of SigP-RxLR-EER genes are:  311
Merged RxLR-EER regex proteins: 74

strain: SCRP324 species: P.rubi
the number of SigP genes is:    89423
the number of SigP-RxLR genes are:      2583
the number of SigP-RxLR-EER genes are:  304
Merged RxLR-EER regex proteins: 276

strain: SCRP333 species: P.rubi
the number of SigP genes is:    81434
the number of SigP-RxLR genes are:      2471
the number of SigP-RxLR-EER genes are:  298
Merged RxLR-EER regex proteins: 270

SPAdes:
strain: SCRP249 species: P.rubi
the number of SigP genes is:    82,268
the number of SigP-RxLR genes are:      2,616
the number of SigP-RxLR-EER genes are:  308
Merged RxLR regex proteins: 2,158
Merged RxLR-EER regex proteins: 277

strain: SCRP324 species: P.rubi
the number of SigP genes is:    82,082
the number of SigP-RxLR genes are:      2,613
the number of SigP-RxLR-EER genes are:  298
Merged RxLR regex proteins: 2,156
Merged RxLR-EER regex proteins: 270

strain: SCRP333 species: P.rubi
the number of SigP genes is:    82,132
the number of SigP-RxLR genes are:      2,594
the number of SigP-RxLR-EER genes are:  299
Merged RxLR regex proteins: 2,141
Merged RxLR-EER regex proteins: 272
```

####H.1) From ORF gene models - Hmm evidence of WY domains

Hmm models for the WY domain contained in many RxLRs were used to search ORFs predicted with atg.pl. These were run with the following commands:

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/P.rubi/*/*_all_secreted.fa)
do
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer
    HmmModel=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
    Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_ORF_WY_hmmer.txt
    hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
    echo "$Organism $Strain"
    cat $OutDir/$HmmResults | grep 'Initial search space'
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
    HmmFasta="$Strain"_ORF_WY_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
    Headers="$Strain"_ORF_WY_hmmer_headers.txt
    cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
    SigP_Merged_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_merged.gff
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_WY_hmmer.gff
    echo "$Strain done"
done
```

```
Discovar:
P.rubi SCRP249
Initial search space (Z):             81,145  [actual number of targets]
Domain search space  (domZ):             416  [number of targets reported over threshold]
P.rubi SCRP324
Initial search space (Z):             89,423  [actual number of targets]
Domain search space  (domZ):             419  [number of targets reported over threshold]
P.rubi SCRP333
Initial search space (Z):             81,434  [actual number of targets]
Domain search space  (domZ):             387  [number of targets reported over threshold]

SPAdes:
P.rubi SCRP249
Initial search space (Z):             82,268  [actual number of targets]
Domain search space  (domZ):             402  [number of targets reported over threshold]
P.rubi SCRP324
Initial search space (Z):             82,082  [actual number of targets]
Domain search space  (domZ):             410  [number of targets reported over threshold]
P.rubi SCRP333
Initial search space (Z):             82,132  [actual number of targets]
Domain search space  (domZ):             388  [number of targets reported over threshold]
```

####H.2) From ORF gene models - Hmm evidence of RxLR effectors

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/P.rubi/*/*_all_secreted.fa)
do
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer
    HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
    Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_ORF_RxLR_hmmer_unmerged.txt
    hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
    echo "$Organism $Strain"
    cat $OutDir/$HmmResults | grep 'Initial search space'
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
    HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
    Headers="$Strain"_ORF_RxLR_hmmer_headers_unmerged.txt
    cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
    SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3
    RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.gff
    RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.txt
    RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_hmm_merged.aa
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_hmmer_unmerged.gff3 --db sigP_ORF_RxLR_hmm.db
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR_hmm.db --id sigP_ORF_RxLR_hmm --out sigP_ORF_RxLR_hmm_merged.db --gff > $RxLR_Merged_Gff
    cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
    ORF_fasta=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
    $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
    printf "Merged RxLR-EER Hmm proteins:\t"
    cat $RxLR_Merged_AA | grep '>' | wc -l
    echo "$Strain done"
done
```

```
Discovar:
P.rubi SCRP249
Initial search space (Z):             81,145  [actual number of targets]
Domain search space  (domZ):             703  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   50
P.rubi SCRP324
Initial search space (Z):             89,423  [actual number of targets]
Domain search space  (domZ):             696  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   217
P.rubi SCRP333
Initial search space (Z):             81,434  [actual number of targets]
Domain search space  (domZ):             668  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   208

SPAdes:
P.rubi SCRP249
Initial search space (Z):             82,268  [actual number of targets]
Domain search space  (domZ):             684  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   211
P.rubi SCRP324
Initial search space (Z):             82,082  [actual number of targets]
Domain search space  (domZ):             682  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   212
P.rubi SCRP333
Initial search space (Z):             82,132  [actual number of targets]
Domain search space  (domZ):             664  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   206
```

####Combining RxLRs from Regex and hmm searches

The total RxLRs are

```bash
#Without EER
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_ORF_RxLR_regex_merged.txt)
do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF.gff3)
    Proteome=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
    HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmm_merged.txt)
    echo "$Organism - $Strain"
    echo "Number of RxLRs identified by Regex:"
    cat $RegexRxLR | sort | uniq | wc -l
    echo "Number of RxLRs identified by Hmm:"
    cat $HmmRxLR | sort | uniq | wc -l
    echo "Number of RxLRs in combined dataset:"
    cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
    OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
    mkdir -p $OutDir
    cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_total_ORF_RxLR_headers.txt $Gff ORF_RxLR Name Augustus > $OutDir/"$Strain"_total_ORF_RxLR.gff
    echo "Number of genes in the extracted gff file:"
    cat $OutDir/"$Strain"_total_ORF_RxLR.gff | grep -w 'gene' | wc -l
    echo ""
    echo "$Strain done without EER"
done

#With EER
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*_ORF_RxLR_EER_regex_merged.txt)
do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain"_ORF.gff3)
    Proteome=$(ls gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa)
    HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain/"$Strain"_ORF_RxLR_hmm_merged.txt)
    echo "$Organism - $Strain"
    echo "Number of RxLRs identified by Regex:"
    cat $RegexRxLR | sort | uniq | wc -l
    echo "Number of RxLRs identified by Hmm:"
    cat $HmmRxLR | sort | uniq | wc -l
    echo "Number of RxLRs in combined dataset:"
    cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
    OutDir=analysis/RxLR_effectors/combined_evidence/$Organism/$Strain
    mkdir -p $OutDir
    cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_EER_headers.txt
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_total_ORF_RxLR_EER_headers.txt $Gff ORF_RxLR Name Augustus > $OutDir/"$Strain"_total_ORF_RxLR_EER.gff
    echo "Number of genes in the extracted gff file:"
    cat $OutDir/"$Strain"_total_ORF_RxLR_EER.gff | grep -w 'gene' | wc -l
    echo ""
    echo "$Strain done with EER"
done
```

```
With EER:
Discovar:
P.rubi - SCRP249
Number of RxLRs identified by Regex:
74
Number of RxLRs identified by Hmm:
50
Number of RxLRs in combined dataset:
84
Number of genes in the extracted gff file:
84

P.rubi - SCRP324
Number of RxLRs identified by Regex:
276
Number of RxLRs identified by Hmm:
217
Number of RxLRs in combined dataset:
309
Number of genes in the extracted gff file:
309

P.rubi - SCRP333
Number of RxLRs identified by Regex:
270
Number of RxLRs identified by Hmm:
208
Number of RxLRs in combined dataset:
298
Number of genes in the extracted gff file:
298

SPAdes:
P.rubi - SCRP249
Number of RxLRs identified by Regex:
277
Number of RxLRs identified by Hmm:
211
Number of RxLRs in combined dataset:
306
Number of genes in the extracted gff file:
306

P.rubi - SCRP324
Number of RxLRs identified by Regex:
270
Number of RxLRs identified by Hmm:
212
Number of RxLRs in combined dataset:
302
Number of genes in the extracted gff file:
302

P.rubi - SCRP333
Number of RxLRs identified by Regex:
272
Number of RxLRs identified by Hmm:
206
Number of RxLRs in combined dataset:
299
Number of genes in the extracted gff file:
299

Without EER:
P.rubi - SCRP249
Number of RxLRs identified by Regex:
2,158
Number of RxLRs identified by Hmm:
211
Number of RxLRs in combined dataset:
2,172
Number of genes in the extracted gff file:
2,172

P.rubi - SCRP324
Number of RxLRs identified by Regex:
2,156
Number of RxLRs identified by Hmm:
212
Number of RxLRs in combined dataset:
2,171
Number of genes in the extracted gff file:
2,171

P.rubi - SCRP333
Number of RxLRs identified by Regex:
2,141
Number of RxLRs identified by Hmm:
206
Number of RxLRs in combined dataset:
2,153
Number of genes in the extracted gff file:
2,153
```

####Analysis of RxLR effectors - merger of Augustus / published genes with ORFs

Intersection between the coodinates of putative RxLRs from gene models and ORFs were identified to determine the total number of RxLRs predicted in these genomes.

The RxLR effectors from both Gene models and ORF finding approaches were combined into a single file.

```bash
#Without EER
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/*/*)
do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$MergeDir/"$Strain"_total_RxLR.gff
    AugTxt=$MergeDir/"$Strain"_total_RxLR_headers.txt
    AugFa=$(ls gene_pred/final/"$Species"/"$Strain"/final/final_genes_appended_renamed.pep.fasta)

    ORFGff=$(ls $MergeDir/"$Strain"_total_ORF_RxLR.gff)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
    ORFsTxt=$(ls $MergeDir/"$Strain"_total_ORF_RxLR_headers.txt)

    ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_RxLR_motif_hmm.gff
    AugInORFs=$MergeDir/"$Strain"_AugInORFs_RxLR_motif_hmm.gff
    ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_RxLR_motif_hmm.gff
    AugUniq=$MergeDir/"$Strain"_Aug_Uniq_RxLR_motif_hmm.gff
    TotalRxLRsTxt=$MergeDir/"$Strain"_Total_RxLR_motif_hmm.txt
    TotalRxLRsGff=$MergeDir/"$Strain"_Total_RxLR_motif_hmm.gff

    bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
    bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
    bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
    bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq

    echo "$Species - $Strain"
    echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
    cat $ORFsInAug | grep -w 'gene' | wc -l
    echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
    cat $AugInORFs | grep -w 'gene' | wc -l
    echo "The number of RxLRs unique to ORF models:"
    cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l
    echo "The number of RxLRs unique to Augustus models:"
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
    echo "The total number of putative RxLRs are:"
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalRxLRsTxt
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
    cat $ORFsUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f3 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
    cat $TotalRxLRsTxt | wc -l
    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff

    RxLRsFa=$MergeDir/"$Strain"_final_RxLR.fa
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalRxLRsTxt > $RxLRsFa
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalRxLRsTxt >> $RxLRsFa
    echo "The number of sequences extracted is"
    cat $RxLRsFa | grep '>' | wc -l
    echo "$Strain done without EER"
done

#With EER
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/*/*)
do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$MergeDir/"$Strain"_total_RxLR_EER.gff
    AugTxt=$MergeDir/"$Strain"_total_RxLR_EER_headers.txt
    AugFa=$(ls gene_pred/final/"$Species"/"$Strain"/final/final_genes_appended_renamed.pep.fasta)

    ORFGff=$(ls $MergeDir/"$Strain"_total_ORF_RxLR_EER.gff)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
    ORFsTxt=$(ls $MergeDir/"$Strain"_total_ORF_RxLR_EER_headers.txt)

    ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_RxLR_EER_motif_hmm.gff
    AugInORFs=$MergeDir/"$Strain"_AugInORFs_RxLR_EER_motif_hmm.gff
    ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff
    AugUniq=$MergeDir/"$Strain"_Aug_Uniq_RxLR_EER_motif_hmm.gff
    TotalRxLRsTxt=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.txt
    TotalRxLRsGff=$MergeDir/"$Strain"_Total_RxLR_EER_motif_hmm.gff

    bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
    bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
    bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
    bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq

    echo "$Species - $Strain"
    echo "The number of ORF RxLRs overlapping Augustus RxLRs:"
    cat $ORFsInAug | grep -w 'gene' | wc -l
    echo "The number of Augustus RxLRs overlapping ORF RxLRs:"
    cat $AugInORFs | grep -w 'gene' | wc -l
    echo "The number of RxLRs unique to ORF models:"
    cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l
    echo "The number of RxLRs unique to Augustus models:"
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
    echo "The total number of putative RxLRs are:"
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalRxLRsTxt
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
    cat $ORFsUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f3 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
    cat $TotalRxLRsTxt | wc -l
    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff

    RxLRsFa=$MergeDir/"$Strain"_final_RxLR_EER.fa
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalRxLRsTxt > $RxLRsFa
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalRxLRsTxt >> $RxLRsFa
    echo "The number of sequences extracted is"
    cat $RxLRsFa | grep '>' | wc -l
    echo "$Strain done without EER"
done
```

```
With EER:
Discovar:
P.rubi - SCRP249
The number of ORF RxLRs overlapping Augustus RxLRs:
50
The number of Augustus RxLRs overlapping ORF RxLRs:
50
The number of RxLRs unique to ORF models:
34
The number of RxLRs unique to Augustus models:
195
The total number of putative RxLRs are:
279
The number of sequences extracted is
276
P.rubi - SCRP324
The number of ORF RxLRs overlapping Augustus RxLRs:
211
The number of Augustus RxLRs overlapping ORF RxLRs:
210
The number of RxLRs unique to ORF models:
98
The number of RxLRs unique to Augustus models:
36
The total number of putative RxLRs are:
344
The number of sequences extracted is
339
P.rubi - SCRP333
The number of ORF RxLRs overlapping Augustus RxLRs:
196
The number of Augustus RxLRs overlapping ORF RxLRs:
195
The number of RxLRs unique to ORF models:
102
The number of RxLRs unique to Augustus models:
34
The total number of putative RxLRs are:
331
The number of sequences extracted is
329

SPAdes:
P.rubi - SCRP249
The number of ORF RxLRs overlapping Augustus RxLRs:
203
The number of Augustus RxLRs overlapping ORF RxLRs:
202
The number of RxLRs unique to ORF models:
103
The number of RxLRs unique to Augustus models:
36
The total number of putative RxLRs are:
341
The number of sequences extracted is
341
P.rubi - SCRP324
The number of ORF RxLRs overlapping Augustus RxLRs:
201
The number of Augustus RxLRs overlapping ORF RxLRs:
199
The number of RxLRs unique to ORF models:
101
The number of RxLRs unique to Augustus models:
31
The total number of putative RxLRs are:
331
The number of sequences extracted is
331
P.rubi - SCRP333
The number of ORF RxLRs overlapping Augustus RxLRs:
209
The number of Augustus RxLRs overlapping ORF RxLRs:
208
The number of RxLRs unique to ORF models:
90
The number of RxLRs unique to Augustus models:
36
The total number of putative RxLRs are:
334
The number of sequences extracted is
334

Without EER:
P.rubi - SCRP249
The number of ORF RxLRs overlapping Augustus RxLRs:
321
The number of Augustus RxLRs overlapping ORF RxLRs:
316
The number of RxLRs unique to ORF models:
1,851
The number of RxLRs unique to Augustus models:
93
The total number of putative RxLRs are:
2,261
The number of sequences extracted is
2,261
P.rubi - SCRP324
The number of ORF RxLRs overlapping Augustus RxLRs:
321
The number of Augustus RxLRs overlapping ORF RxLRs:
315
The number of RxLRs unique to ORF models:
1,850
The number of RxLRs unique to Augustus models:
82
The total number of putative RxLRs are:
2,247
The number of sequences extracted is
2,247
P.rubi - SCRP333
The number of ORF RxLRs overlapping Augustus RxLRs:
330
The number of Augustus RxLRs overlapping ORF RxLRs:
323
The number of RxLRs unique to ORF models:
1,823
The number of RxLRs unique to Augustus models:
88
The total number of putative RxLRs are:
2,234
The number of sequences extracted is
2,234
```

####I) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers in ORF gene models. This was done with the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/P.rubi/*/*.aa_cat.fa)
do
    # Setting variables
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/CRN_effectors/hmmer_CRN/$Organism/$Strain
    mkdir -p $OutDir
    # Hmmer variables
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer
    HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
    # Searches for LFLAK domain
    LFLAK_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm
    HmmResultsLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.txt
    hmmsearch -T 0 $LFLAK_hmm $Proteome > $OutDir/$HmmResultsLFLAK
    echo "Searching for LFLAK domains in: $Organism $Strain"
    cat $OutDir/$HmmResultsLFLAK | grep 'Initial search space'
    cat $OutDir/$HmmResultsLFLAK | grep 'number of targets reported over threshold'
    HmmFastaLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsLFLAK $Proteome > $OutDir/$HmmFastaLFLAK
    # Searches for DWL domain
    DWL_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm
    HmmResultsDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.txt
    hmmsearch -T 0 $DWL_hmm $Proteome > $OutDir/$HmmResultsDWL
    echo "Searching for DWL domains in: $Organism $Strain"
    cat $OutDir/$HmmResultsDWL | grep 'Initial search space'
    cat $OutDir/$HmmResultsDWL | grep 'number of targets reported over threshold'
    HmmFastaDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsDWL $Proteome > $OutDir/$HmmFastaDWL
    # Identify ORFs found by both models
    CommonHeaders=$OutDir/"$Strain"_ORF_CRN_DWL_LFLAK_unmerged_headers.txt
    cat $OutDir/$HmmFastaLFLAK $OutDir/$HmmFastaDWL | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $CommonHeaders
    echo "The number of CRNs common to both models are:"
    cat $CommonHeaders | wc -l
    # The sequences will be merged based upon the strength of their DWL domain score
    # For this reason headers as they appear in the DWL fasta file were extracted
    Headers=$OutDir/"$Strain"_CRN_hmmer_unmerged_headers.txt
    cat $OutDir/$HmmFastaDWL | grep '>' | grep -w -f $CommonHeaders | tr -d '>' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | tr -d '-' | sed 's/hmm_score/HMM_score/g' > $Headers
    # As we are dealing with JGI and Broad sequences, some features need formatting:
    ORF_Gff=$(ls gene_pred/ORF_finder/$Organism/$Strain/*_ORF.gff3)
    # Gff features were extracted for each header
    CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $Headers $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
    # Gff features were merged based upon the DWL hmm score
    DbDir=analysis/databases/$Organism/$Strain
    mkdir -p $DbDir
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $CRN_unmerged_Gff --db $DbDir/CRN_ORF.db
    CRN_Merged_Gff=$OutDir/"$Strain"_CRN_merged_hmmer.gff3
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp $DbDir/CRN_ORF.db --id LFLAK_DWL_CRN --out $DbDir/CRN_ORF_merged.db --gff > $CRN_Merged_Gff
    # Final results are reported:
    echo "Number of CRN ORFs after merging:"
    cat $CRN_Merged_Gff | grep 'gene' | wc -l
    echo "$Strain done"
done
```

```
Discovar:
Searching for LFLAK domains in: P.rubi SCRP249
Initial search space (Z):            716,788  [actual number of targets]
Domain search space  (domZ):             311  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP249
Initial search space (Z):            716,788  [actual number of targets]
Domain search space  (domZ):             512  [number of targets reported over threshold]
The number of CRNs common to both models are:
233
Number of CRN ORFs after merging:
141
Searching for LFLAK domains in: P.rubi SCRP324
Initial search space (Z):            795,934  [actual number of targets]
Domain search space  (domZ):             316  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP324
Initial search space (Z):            795,934  [actual number of targets]
Domain search space  (domZ):             517  [number of targets reported over threshold]
The number of CRNs common to both models are:
234
Number of CRN ORFs after merging:
140
Searching for LFLAK domains in: P.rubi SCRP333
Initial search space (Z):            720,807  [actual number of targets]
Domain search space  (domZ):             301  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP333
Initial search space (Z):            720,807  [actual number of targets]
Domain search space  (domZ):             468  [number of targets reported over threshold]
The number of CRNs common to both models are:
222
Number of CRN ORFs after merging:
132

SPAdes:
Searching for LFLAK domains in: P.rubi SCRP249
Initial search space (Z):             645,852  [actual number of targets]
Domain search space  (domZ):             301  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP249
Initial search space (Z):             645,852  [actual number of targets]
Domain search space  (domZ):             452  [number of targets reported over threshold]
The number of CRNs common to both models are:
226
Number of CRN ORFs after merging:
131
Searching for LFLAK domains in: P.rubi SCRP324
Initial search space (Z):             647,415  [actual number of targets]
Domain search space  (domZ):             294  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP324
Initial search space (Z):             647,415  [actual number of targets]
Domain search space  (domZ):             440  [number of targets reported over threshold]
The number of CRNs common to both models are:
217
Number of CRN ORFs after merging:
126
Searching for LFLAK domains in: P.rubi SCRP333
Initial search space (Z):             646,159  [actual number of targets]
Domain search space  (domZ):             293  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP333
Initial search space (Z):             646,159  [actual number of targets]
Domain search space  (domZ):             438  [number of targets reported over threshold]
The number of CRNs common to both models are:
220
Number of CRN ORFs after merging:
129
```

####Merge CRNs from ORF fragments and Augustus gene models

```bash
for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/P.rubi/*)
do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$(ls $MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff)
    AugFa=$(ls gene_pred/final/"$Species"/"$Strain"/final/final_genes_appended_renamed.pep.fasta)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
    ORFGff=$MergeDir/"$Strain"_CRN_merged_hmmer.gff3
    ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_CRN_hmmer.bed
    AugInORFs=$MergeDir/"$Strain"_AugInORFs_CRN_hmmer.bed
    ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_CRN_hmmer.bed
    AugUniq=$MergeDir/"$Strain"_Aug_Uniq_CRN_hmmer.bed
    TotalCRNsTxt=$MergeDir/"$Strain"_final_CRN.txt
    TotalCRNsGff=$MergeDir/"$Strain"_final_CRN.gff
    TotalCRNsHeaders=$MergeDir/"$Strain"_Total_CRN_headers.txt
    bedtools intersect -wa -u -a $ORFGff -b $AugGff > $ORFsInAug
    bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
    bedtools intersect -v -wa -a $ORFGff -b $AugGff > $ORFsUniq
    bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq
    echo "$Species - $Strain"

    echo "The number of ORF CRNs overlapping Augustus CRNs:"
    cat $ORFsInAug | grep -w -e 'transcript' -e 'mRNA' | wc -l
    echo "The number of Augustus CRNs overlapping ORF CRNs:"
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA' | wc -l
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalCRNsTxt
    echo "The number of CRNs unique to ORF models:"
    cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' | wc -l
    cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt
    echo "The number of CRNs unique to Augustus models:"
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt

    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalCRNsTxt > $TotalCRNsGff

    CRNsFa=$MergeDir/"$Strain"_final_CRN.fa
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalCRNsTxt > $CRNsFa
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalCRNsTxt >> $CRNsFa
    echo "The number of sequences extracted is"
    cat $CRNsFa | grep '>' | wc -l
    echo "$Strain done"
done
```

```
Discovar:
P.rubi - SCRP249
The number of ORF CRNs overlapping Augustus CRNs:
132
The number of Augustus CRNs overlapping ORF CRNs:
132
The number of CRNs unique to ORF models:
9
The number of CRNs unique to Augustus models:
3
The number of sequences extracted is
132
P.rubi - SCRP324
The number of ORF CRNs overlapping Augustus CRNs:
132
The number of Augustus CRNs overlapping ORF CRNs:
132
The number of CRNs unique to ORF models:
8
The number of CRNs unique to Augustus models:
3
The number of sequences extracted is
133
P.rubi - SCRP333
The number of ORF CRNs overlapping Augustus CRNs:
121
The number of Augustus CRNs overlapping ORF CRNs:
121
The number of CRNs unique to ORF models:
11
The number of CRNs unique to Augustus models:
2
The number of sequences extracted is
129

SPAdes:
P.rubi - SCRP249
The number of ORF CRNs overlapping Augustus CRNs:
122
The number of Augustus CRNs overlapping ORF CRNs:
122
The number of CRNs unique to ORF models:
9
The number of CRNs unique to Augustus models:
7
The number of sequences extracted is
138
P.rubi - SCRP324
The number of ORF CRNs overlapping Augustus CRNs:
66
The number of Augustus CRNs overlapping ORF CRNs:
68
The number of CRNs unique to ORF models:
60
The number of CRNs unique to Augustus models:
2
The number of sequences extracted is
130
P.rubi - SCRP333
The number of ORF CRNs overlapping Augustus CRNs:
121
The number of Augustus CRNs overlapping ORF CRNs:
123
The number of CRNs unique to ORF models:
8
The number of CRNs unique to Augustus models:
6
The number of sequences extracted is
137
```

####J) From ORF fragments - ApoplastP prediction of apoplastic effectors

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted_merged.aa)
do
    Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
    BaseName="$Organism"_"$Strain"_ApoplastP_ORF
    echo "$Organism - $Strain"
    OutDir=analysis/ApoplastP/$Organism/$Strain
    ProgDir=/home/adamst/git_repos/tools/seq_tools/apoplastic_effectors
    qsub $ProgDir/pred_apoplastP.sh $Secretome $BaseName $OutDir
done
```

TEMPORARY COMMANDS FOR RUNNING ON HEAD NODE

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*_all_secreted.fa)
do
    Strain=$(echo $Secretome | rev | cut -f2 -d "/" | rev)
    Organism=$(echo $Secretome | rev | cut -f3 -d "/" | rev)
    echo "$Organism - $Strain"
    BaseName="$Organism"_"$Strain"_ApoplastP_ORF
    OutDir=analysis/ApoplastP/$Organism/$Strain
    mkdir -p $OutDir
    ApoplastP.py -o "$OutDir"/"$BaseName".txt -A "$OutDir"/"$BaseName"_unmerged.fa -i $Secretome
done
```

The number of proteins predicted as being apoplastic effectors were summarised using the following commands

Gffutils seems only to search the same strand as the feature

Occassionally there is an error in the script, but this does not happen when repeating the commands

```bash
for File in $(ls analysis/ApoplastP/*/*/*_ApoplastP_ORF_unmerged.fa)
do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo $File | sed 's/_ApoplastP_ORF_unmerged.fa/_ApoplastP_ORF_unmerged_headers.txt/g')
    Gff=$(echo $File | sed 's/_ApoplastP_ORF_unmerged.fa/_ApoplastP_ORF_unmerged.gff3/g')
    echo "Creating Headers file"
    cat $File | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*$//' | tr -d ' ' > $Headers
    SigP_Gff=gene_pred/combined_sigP_ORF/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
    ORF_fasta=$(ls gene_pred/ORF_finder/*/$Strain/"$Strain".aa_cat.fa)
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $Headers $SigP_Gff ApoplastP_ORF Name Augustus > $Gff
    Apo_Merged_Gff_1_strand=analysis/ApoplastP/$Organism/$Strain/"$Strain"_ApoplastP_ORF_merged_1_strand.gff
    Apo_Merged_txt=analysis/ApoplastP/$Organism/$Strain/"$Strain"_ApoplastP_ORF_merged_headers.txt
    Apo_Merged_AA=analysis/ApoplastP/$Organism/$Strain/"$Strain"_ApoplastP_ORF_merged.aa
    Apo_Merged_Gff=analysis/ApoplastP/$Organism/$Strain/"$Strain"_ApoplastP_ORF_merged.gff
    Apo_Merged_Gff_no_strands=analysis/ApoplastP/$Organism/$Strain/"$Strain"_ApoplastP_ORF_merged_no_strands.gff
    # Remove features on the same strand
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $Gff --db sigP_ORF_ApoP.db
    echo "First database made"
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_ApoP.db --id sigP_ORF_ApoplastP --out sigP_ORF_ApoP_merged.db --gff > $Apo_Merged_Gff_1_strand
    echo "Overlapping features on the same strand removed"
    # Modify Gff to remove all strand specification
    Gff_no_strands=analysis/ApoplastP/$Organism/$Strain/"$Strain"_ORF_ApoplastP_no_strands.gff3
    cat $Gff | sed 's/+/./g' | sed 's/-/./g' > $Gff_no_strands
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $Gff_no_strands --db sigP_ORF_ApoP_mod.db
    echo "Second database made"
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_ApoP_mod.db --id sigP_ORF_ApoplastP --out sigP_ORF_ApoP_merged_mod.db --gff > $Apo_Merged_Gff_no_strands
    echo "Overlapping features on opposite strands removed"
    cat $Apo_Merged_Gff_no_strands | grep 'transcript' | rev | cut -f1 -d '=' | rev > $Apo_Merged_txt
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    # Easy reformatting
    Apo_Merged_Gff_temp=analysis/ApoplastP/$Organism/$Strain/"$Strain"_ApoplastP_ORF_merged_unformatted.gff
    $ProgDir/gene_list_to_gff.pl $Apo_Merged_txt $Gff ApoplastP_ORF Name Augustus > $Apo_Merged_Gff_temp
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $Apo_Merged_Gff_temp --db sigP_ORF_ApoP_reformat.db
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    echo "Reformatting begun"
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_ApoP_reformat.db --id sigP_ORF_ApoplastP --out sigP_ORF_ApoP_merged_reformat.db --gff > $Apo_Merged_Gff
    echo "Reformatting finished"
    echo "The number of genes predicted as Apoplastic effectors is:"
    cat $Apo_Merged_txt | wc -l
    echo "The number of genes extracted to the GFF is:"
    cat $Apo_Merged_Gff | grep 'transcript' | wc -l
    ORF_fasta=gene_pred/ORF_finder/$Organism/$Strain/"$Strain".aa_cat.fa
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $Apo_Merged_txt > $Apo_Merged_AA
    echo "The number of genes extracted to the FASTA is:"
    cat $Apo_Merged_AA | grep '>' | wc -l
done
```

```
P.rubi - SCRP249
The number of genes predicted as Apoplastic effectors is:
10,341
The number of genes extracted to the GFF is:
10,341
The number of genes extracted to the FASTA is:
10,341

P.rubi - SCRP324
The number of genes predicted as Apoplastic effectors is:
10,280
The number of genes extracted to the GFF is:
10,280
The number of genes extracted to the FASTA is:
10,280

P.rubi - SCRP333
The number of genes predicted as Apoplastic effectors is:
10,354
The number of genes extracted to the GFF is:
10,354
The number of genes extracted to the FASTA is:
10,354
```

####Merge apoplastic effectors from augustus models and ORF fragments

```bash
for MergeDir in $(ls -d analysis/ApoplastP/*/*)
do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    AugGff=$MergeDir/"$Species"_"$Strain"_ApoplastP.gff
    AugTxt=$MergeDir/"$Species"_"$Strain"_ApoplastP_headers.txt
    AugFa=$(ls gene_pred/final/"$Species"/"$Strain"/final/final_genes_appended_renamed.pep.fasta)
    AugRxLR=$(ls analysis/RxLR_effectors/combined_evidence/$Species/$Strain/"$Strain"_total_RxLR.gff)
    AugCRN=$(ls analysis/CRN_effectors/hmmer_CRN/$Species/$Strain/"$Strain"_pub_CRN_LFLAK_DWL.gff)

    ORFGff=$(ls $MergeDir/"$Strain"_ApoplastP_ORF_merged.gff)
    ORFsFa=$(ls gene_pred/ORF_finder/*/"$Strain"/"$Strain".aa_cat.fa)
    ORFsTxt=$(ls $MergeDir/"$Strain"_ApoplastP_ORF_merged_headers.txt)

    ORFsInAug=$MergeDir/"$Strain"_ORFsInAug_ApoplastP.gff
    AugInORFs=$MergeDir/"$Strain"_AugInORFs_ApoplastP.gff
    ORFsUniq=$MergeDir/"$Strain"_ORFsUniq_ApoplastP.gff
    AugUniq=$MergeDir/"$Strain"_Aug_Uniq_ApoplastP.gff
    TotalApoPTxt=$MergeDir/"$Strain"_Total_ApoplastP.txt
    TotalApoPGff=$MergeDir/"$Strain"_Total_ApoplastP.gff

    bedtools intersect -wa -u -a $ORFGff -b $AugGff $AugRxLR $AugCRN > $ORFsInAug
    bedtools intersect -wa -u -a $AugGff -b $ORFGff > $AugInORFs
    bedtools intersect -v -wa -a $ORFGff -b $AugGff $AugRxLR $AugCRN > $ORFsUniq
    bedtools intersect -v -wa -a $AugGff -b $ORFGff > $AugUniq

    echo "$Species - $Strain"
    echo "The number of ORF apoplastic effectors overlapping Augustus effectors:"
    cat $ORFsInAug | grep -w 'gene' | wc -l
    echo "The number of Augustus effectors overlapping ORF apoplastic effectors:"
    cat $AugInORFs | grep -w 'gene' | wc -l
    echo "The number of apoplastic effectors unique to ORF models:"
    cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l
    echo "The number of apoplastic effectors unique to Augustus models:"
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l
    echo "The total number of putative apoplastic effectors are:"
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalApoPTxt
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalApoPTxt
    cat $ORFsUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' >> $TotalApoPTxt
    cat $TotalApoPTxt | wc -l
    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalApoPTxt > $TotalApoPGff

    ApoPFa=$MergeDir/"$Strain"_final_ApoplastP.fa
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalApoPTxt > $ApoPFa
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalApoPTxt >> $ApoPFa
    echo "The number of sequences extracted is"
    cat $ApoPFa | grep '>' | wc -l
done
```

```
P.rubi - SCRP249
The number of ORF apoplastic effectors overlapping Augustus effectors:
958
The number of Augustus effectors overlapping ORF apoplastic effectors:
0
The number of apoplastic effectors unique to ORF models:
9,383
The number of apoplastic effectors unique to Augustus models:
142
The total number of putative apoplastic effectors are:
10,400
The number of sequences extracted is
10,400

P.rubi - SCRP324
The number of ORF apoplastic effectors overlapping Augustus apoplastic effectors:
979
The number of Augustus apoplastic effectors overlapping ORF apoplastic effectors:
896
The number of apoplastic effectors unique to ORF models:
9,301
The number of apoplastic effectors unique to Augustus models:
169
The total number of putative apoplastic effectors are:
10,373
The number of sequences extracted is
10,373

P.rubi - SCRP333
The number of ORF apoplastic effectors overlapping Augustus apoplastic effectors:
999
The number of Augustus apoplastic effectors overlapping ORF apoplastic effectors:
0
The number of apoplastic effectors unique to ORF models:
9,355
The number of apoplastic effectors unique to Augustus models:
171
The total number of putative apoplastic effectors are:
10,445
The number of sequences extracted is
10,445
```

<!-- Due to an unknown error, the softmasked files for SCRP249 and SCRP324 do not read into the hash table in the add_ORF_features.pl script. Wrapping the unmasked file every 60 characters provides an assembly file that does work.

```bash
for Strain in SCRP249 SCRP324
do
    Assembly=repeat_masked/P.rubi/$Strain/deconseq_Paen_repmask/"$Strain"_contigs_unmasked.fa
    fold -w 60 $Assembly > repeat_masked/P.rubi/$Strain/deconseq_Paen_repmask/"$Strain"_contigs_unmasked_wrapped.fa
done
``` -->

## Before making a combined GFF, several other intersects are needed

Due to high hit rate of ApoP, and the fact that some RxLRs are also hit by ApoP,
Extra intersects are needed to avoid genes being present twice.

```bash
for MergeDir in $(ls -d analysis/ApoplastP/*/*)
do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    RxLR_ORFs=$(ls analysis/RxLR_effectors/combined_evidence/P.rubi/$Strain/"$Strain"_ORFsUniq_RxLR_motif_hmm.gff)
    RxLR_EER_ORFs=$(ls analysis/RxLR_effectors/combined_evidence/P.rubi/$Strain/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff)
    CRN_ORFs=$(ls analysis/CRN_effectors/hmmer_CRN/P.rubi/$Strain/"$Strain"_ORFsUniq_CRN_hmmer.bed)
    ApoP_ORFs=$(ls analysis/ApoplastP/P.rubi/$Strain/"$Strain"_ORFsUniq_ApoplastP.gff)
    RxLR_No_ApoP=$MergeDir/"$Strain"_RxLR_No_ApoP_ORFs.gff
    RxLR_Plus_ApoP=$MergeDir/"$Strain"_RxLR_Plus_ApoP_ORFs.gff
    CRN_No_ApoP=$MergeDir/"$Strain"_CRN_No_ApoP_ORFs.gff
    CRN_Plus_ApoP=$MergeDir/"$Strain"_CRN_Plus_ApoP_ORFs.gff
    ApoP_No_RxLR_CRN=$MergeDir/"$Strain"_ApoP_No_RxLR_CRN_ORFs.gff
    RxLR_EER_No_ApoP=$MergeDir/"$Strain"_RxLR_EER_No_ApoP_ORFs.gff
    RxLR_EER_Plus_ApoP=$MergeDir/"$Strain"_RxLR_EER_Plus_ApoP_ORFs.gff
    ApoP_No_RxLR_EER_CRN=$MergeDir/"$Strain"_ApoP_No_RxLR_EER_CRN_ORFs.gff

    bedtools intersect -v -wa -a $RxLR_ORFs -b $ApoP_ORFs > $RxLR_No_ApoP
    bedtools intersect -wa -u -a $RxLR_ORFs -b $ApoP_ORFs > $RxLR_Plus_ApoP
    bedtools intersect -v -wa -a $CRN_ORFs -b $ApoP_ORFs > $CRN_No_ApoP
    bedtools intersect -wa -u -a $CRN_ORFs -b $ApoP_ORFs > $CRN_Plus_ApoP
    bedtools intersect -v -wa -a $ApoP_ORFs -b $RxLR_ORFs $CRN_ORFs > $ApoP_No_RxLR_CRN
    bedtools intersect -v -wa -a $RxLR_EER_ORFs -b $ApoP_ORFs > $RxLR_EER_No_ApoP
    bedtools intersect -wa -u -a $RxLR_EER_ORFs -b $ApoP_ORFs > $RxLR_EER_Plus_ApoP
    bedtools intersect -v -wa -a $ApoP_ORFs -b $RxLR_EER_ORFs $CRN_ORFs > $ApoP_No_RxLR_EER_CRN

    RxLR_No_ApoP_modified=$MergeDir/"$Strain"_RxLR_No_ApoP_ORFs_modified.gff
    RxLR_Plus_ApoP_modified=$MergeDir/"$Strain"_RxLR_Plus_ApoP_ORFs_modified.gff
    CRN_No_ApoP_modified=$MergeDir/"$Strain"_CRN_No_ApoP_ORFs_modified.gff
    CRN_Plus_ApoP_modified=$MergeDir/"$Strain"_CRN_Plus_ApoP_ORFs_modified.gff
    RxLR_EER_No_ApoP_modified=$MergeDir/"$Strain"_RxLR_EER_No_ApoP_ORFs_modified.gff
    RxLR_EER_Plus_ApoP_modified=$MergeDir/"$Strain"_RxLR_EER_Plus_ApoP_ORFs_modified.gff

    cat $RxLR_No_ApoP | sed -r 's/ORF_RxLR/RxLR_ORF/g' > $RxLR_No_ApoP_modified
    cat $RxLR_Plus_ApoP | sed -r 's/ORF_RxLR/RxLR_ORF_+_ApoplastP_ORF/g' > $RxLR_Plus_ApoP_modified
    cat $CRN_No_ApoP | sed -r 's/CRN_HMM/CRN_ORF/g' > $CRN_No_ApoP_modified
    cat $CRN_Plus_ApoP | sed -r 's/CRN_HMM/CRN_ORF_+_ApoplastP_ORF/g' > $CRN_Plus_ApoP_modified
    cat $RxLR_EER_No_ApoP | sed -r 's/ORF_RxLR/RxLR_ORF/g' > $RxLR_EER_No_ApoP_modified
    cat $RxLR_EER_Plus_ApoP | sed -r 's/ORF_RxLR/RxLR_ORF_+_ApoplastP_ORF/g' > $RxLR_EER_Plus_ApoP_modified

    mv $RxLR_No_ApoP_modified $RxLR_No_ApoP
    mv $RxLR_Plus_ApoP_modified $RxLR_Plus_ApoP
    mv $CRN_No_ApoP_modified $CRN_No_ApoP
    mv $CRN_Plus_ApoP_modified $CRN_Plus_ApoP
    mv $RxLR_EER_No_ApoP_modified $RxLR_EER_No_ApoP
    mv $RxLR_EER_Plus_ApoP_modified $RxLR_EER_Plus_ApoP

    echo "$Species - $Strain"
    echo "The number of RxLRs not hit by ApoplastP is:"
    cat $RxLR_No_ApoP | grep -w 'gene' | wc -l
    echo "The number of RxLRs that also hit ApoplastP is:"
    cat $RxLR_Plus_ApoP | grep -w 'gene' | wc -l
    echo "The number of CRNs not hit by ApoplastP is:"
    cat $CRN_No_ApoP | grep -w 'gene' | wc -l
    echo "The number of CRNs that also hit ApoplastP is:"
    cat $CRN_Plus_ApoP | grep -w 'gene' | wc -l
    echo "The number of ApoplastP hits that are not RxLRs or CRNs is:"
    cat $ApoP_No_RxLR_CRN | grep -w 'gene' | wc -l
    echo "The number of RxLRs, with an EER domain, not hit by ApoplastP is:"
    cat $RxLR_EER_No_ApoP | grep -w 'gene' | wc -l
    echo "The number of RxLRs, with an EER domain, also hit by ApoplastP is:"
    cat $RxLR_EER_Plus_ApoP | grep -w 'gene' | wc -l
    echo "The number of ApoplastP hits that are not RxLRs, with an EER domain, or CRNs is:"
    cat $ApoP_No_RxLR_EER_CRN | grep -w 'gene' | wc -l
    printf "\n"
done
```

```
P.rubi - SCRP249
The number of RxLRs not hit by ApoplastP is:
1,287
The number of RxLRs that also hit ApoplastP is:
564
The number of CRNs not hit by ApoplastP is:
8
The number of CRNs that also hit ApoplastP is:
1
The number of ApoplastP hits that are not RxLRs or CRNs is:
8,817
The number of RxLRs, with an EER domain, not hit by ApoplastP is:
93
The number of RxLRs, with an EER domain, also hit by ApoplastP is:
10
The number of ApoplastP hits that are not RxLRs, with an EER domain, or CRNs is:
9,371

P.rubi - SCRP324
The number of RxLRs not hit by ApoplastP is:
1,294
The number of RxLRs that also hit ApoplastP is:
556
The number of CRNs not hit by ApoplastP is:
56
The number of CRNs that also hit ApoplastP is:
4
The number of ApoplastP hits that are not RxLRs or CRNs is:
8,738
The number of RxLRs, with an EER domain, not hit by ApoplastP is:
89
The number of RxLRs, with an EER domain, also hit by ApoplastP is:
12
The number of ApoplastP hits that are not RxLRs, with an EER domain, or CRNs is:
9,283

P.rubi - SCRP333
The number of RxLRs not hit by ApoplastP is:
1,262
The number of RxLRs that also hit ApoplastP is:
561
The number of CRNs not hit by ApoplastP is:
8
The number of CRNs that also hit ApoplastP is:
0
The number of ApoplastP hits that are not RxLRs or CRNs is:
8,793
The number of RxLRs, with an EER domain, not hit by ApoplastP is:
80
The number of RxLRs, with an EER domain, also hit by ApoplastP is:
10
The number of ApoplastP hits that are not RxLRs, with an EER domain, or CRNs is:
9,344
```

## Making a combined file of Braker and CodingQuary genes with additional ORF effector candidates

First checks features file isn't empty with a series of if statements

```bash
for num in 1
do
    #Without EER discrimination and with ApoplastP
    for GeneGff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3)
    do
        Strain=$(echo $GeneGff | rev | cut -d '/' -f3 | rev)
        Organism=$(echo $GeneGff | rev | cut -d '/' -f4 | rev)
        MergeDir=$(ls -d analysis/ApoplastP/$Organism/$Strain)
        echo "$Strain - $Organism"
        RxLR_No_ApoP=$MergeDir/"$Strain"_RxLR_No_ApoP_ORFs.gff
        RxLR_Plus_ApoP=$MergeDir/"$Strain"_RxLR_Plus_ApoP_ORFs.gff
        CRN_No_ApoP=$MergeDir/"$Strain"_CRN_No_ApoP_ORFs.gff
        CRN_Plus_ApoP=$MergeDir/"$Strain"_CRN_Plus_ApoP_ORFs.gff
        ApoP_No_RxLR_CRN=$MergeDir/"$Strain"_ApoP_No_RxLR_CRN_ORFs.gff
        if [ -f repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked.fa)
            echo $Assembly
        elif [ -f repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked.fa)
            echo $Assembly
        elif [ -f repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked.fa)
            echo $Assembly
        else
            Assembly=$(ls repeat_masked/quiver_results/polished/filtered_contigs_repmask/*_softmasked.fa)
            echo $Assembly
        fi
        OutDir=gene_pred/annotation/P.rubi/$Strain
        mkdir -p $OutDir
        cat $GeneGff  \
        	| sed 's/\(\tCDS\t.*\)transcript_id "\(.*\)"; gene_id.*/\1ID=\2.CDS; Parent=\2/g' \
        	| sed 's/\(\exon\t.*\)transcript_id "\(.*\)"; gene_id.*/\1ID=\2.exon; Parent=\2/g' \
        	| sed 's/transcript_id "/ID=/g' | sed 's/";/;/g' | sed 's/ gene_id "/Parent=/g' \
        	| sed -r "s/\tg/\tID=g/g" | sed 's/ID=gene/gene/g' | sed -r "s/;$//g" \
        	| sed "s/\ttranscript\t.*ID=\(.*\).t.*$/\0;Parent=\1/" \
        	> $OutDir/"$Strain"_genes_incl_ORFeffectors.gff3
        # cat $GeneGff > $OutDir/10300_genes_incl_ORFeffectors.gff3
        ProgDir=/home/adamst/git_repos/scripts/phytophthora/10300_analysis
        $ProgDir/gff_name2id.py --gff $RxLR_No_ApoP > $OutDir/ORF_RxLR_No_ApoP_parsed.gff3
        $ProgDir/gff_name2id.py --gff $RxLR_Plus_ApoP > $OutDir/ORF_RxLR_Plus_ApoP_parsed.gff3
        $ProgDir/gff_name2id.py --gff $CRN_No_ApoP > $OutDir/ORF_CRN_No_ApoP_parsed.gff3
        $ProgDir/gff_name2id.py --gff $CRN_Plus_ApoP > $OutDir/ORF_CRN_Plus_ApoP_parsed.gff3
        $ProgDir/gff_name2id.py --gff $ApoP_No_RxLR_CRN > $OutDir/ORF_ApoP_No_RxLR_CRN_parsed.gff3

        ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
        if [ -s $OutDir/ORF_RxLR_No_ApoP_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_RxLR_No_ApoP_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors.gff3
        fi
        if [ -s $OutDir/ORF_RxLR_Plus_ApoP_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_RxLR_Plus_ApoP_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors.gff3
        fi
        if [ -s $OutDir/ORF_CRN_No_ApoP_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_CRN_No_ApoP_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors.gff3
        fi
        if [ -s $OutDir/ORF_CRN_Plus_ApoP_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_CRN_Plus_ApoP_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors.gff3
        fi
        if [ -s $OutDir/ORF_ApoP_No_RxLR_CRN_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_ApoP_No_RxLR_CRN_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors.gff3
        fi
        # Make gene models from gff files.
        ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
        $ProgDir/gff2fasta.pl $Assembly $OutDir/"$Strain"_genes_incl_ORFeffectors.gff3 $OutDir/"$Strain"_genes_incl_ORFeffectors
        # Note - these fasta files have not been validated - do not use
    done
    echo "Gff created with all RxLRs and ApoplastP hits"

    #With EER discrimination and with ApoplastP
    for GeneGff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3)
    do
        Strain=$(echo $GeneGff | rev | cut -d '/' -f3 | rev)
        Organism=$(echo $GeneGff | rev | cut -d '/' -f4 | rev)
        MergeDir=$(ls -d analysis/ApoplastP/$Organism/$Strain)
        echo "$Strain - $Organism"
        RxLR_EER_No_ApoP=$MergeDir/"$Strain"_RxLR_EER_No_ApoP_ORFs.gff
        RxLR_EER_Plus_ApoP=$MergeDir/"$Strain"_RxLR_EER_Plus_ApoP_ORFs.gff
        CRN_No_ApoP=$MergeDir/"$Strain"_CRN_No_ApoP_ORFs.gff
        CRN_Plus_ApoP=$MergeDir/"$Strain"_CRN_Plus_ApoP_ORFs.gff
        ApoP_No_RxLR_EER_CRN=$MergeDir/"$Strain"_ApoP_No_RxLR_EER_CRN_ORFs.gff
        if [ -f repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked.fa)
            echo $Assembly
        elif [ -f repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked.fa)
            echo $Assembly
        elif [ -f repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked.fa)
            echo $Assembly
        else
            Assembly=$(ls repeat_masked/quiver_results/polished/filtered_contigs_repmask/*_softmasked.fa)
            echo $Assembly
        fi
        OutDir=gene_pred/annotation/P.rubi/$Strain
        mkdir -p $OutDir
        cat $GeneGff  \
        	| sed 's/\(\tCDS\t.*\)transcript_id "\(.*\)"; gene_id.*/\1ID=\2.CDS; Parent=\2/g' \
        	| sed 's/\(\exon\t.*\)transcript_id "\(.*\)"; gene_id.*/\1ID=\2.exon; Parent=\2/g' \
        	| sed 's/transcript_id "/ID=/g' | sed 's/";/;/g' | sed 's/ gene_id "/Parent=/g' \
        	| sed -r "s/\tg/\tID=g/g" | sed 's/ID=gene/gene/g' | sed -r "s/;$//g" \
        	| sed "s/\ttranscript\t.*ID=\(.*\).t.*$/\0;Parent=\1/" \
        	> $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative.gff3
        # cat $GeneGff > $OutDir/10300_genes_incl_ORFeffectors.gff3
        ProgDir=/home/adamst/git_repos/scripts/phytophthora/10300_analysis
        $ProgDir/gff_name2id.py --gff $RxLR_EER_No_ApoP > $OutDir/ORF_RxLR_EER_No_ApoP_parsed.gff3
        $ProgDir/gff_name2id.py --gff $RxLR_EER_Plus_ApoP > $OutDir/ORF_RxLR_EER_Plus_ApoP_parsed.gff3
        $ProgDir/gff_name2id.py --gff $CRN_No_ApoP > $OutDir/ORF_CRN_No_ApoP_parsed.gff3
        $ProgDir/gff_name2id.py --gff $CRN_Plus_ApoP > $OutDir/ORF_CRN_Plus_ApoP_parsed.gff3
        $ProgDir/gff_name2id.py --gff $ApoP_No_RxLR_EER_CRN > $OutDir/ORF_ApoP_No_RxLR_EER_CRN_parsed.gff3

        ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
        if [ -s $OutDir/ORF_RxLR_EER_No_ApoP_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_RxLR_EER_No_ApoP_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative.gff3
        fi
        if [ -s $OutDir/ORF_RxLR_EER_Plus_ApoP_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_RxLR_EER_Plus_ApoP_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative.gff3
        fi
        if [ -s $OutDir/ORF_CRN_No_ApoP_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_CRN_No_ApoP_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative.gff3
        fi
        if [ -s $OutDir/ORF_CRN_Plus_ApoP_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_CRN_Plus_ApoP_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative.gff3
        fi
        if [ -s $OutDir/ORF_ApoP_No_RxLR_EER_CRN_parsed.gff3 ]
        then
            $ProgDir/add_ORF_features.pl $OutDir/ORF_ApoP_No_RxLR_EER_CRN_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative.gff3
        fi
        # Make gene models from gff files.
        ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
        $ProgDir/gff2fasta.pl $Assembly $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative.gff3 $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative
        # Note - these fasta files have not been validated - do not use
    done
    echo "Gff created with high confidence RxLRs and ApoplastP hits"

    #Without EER discrimination and without ApoplastP
    for GeneGff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3)
    do
        Strain=$(echo $GeneGff | rev | cut -d '/' -f3 | rev)
        Organism=$(echo $GeneGff | rev | cut -d '/' -f4 | rev)
        echo "$Strain - $Organism"
        GffOrfRxLR=$(ls analysis/RxLR_effectors/combined_evidence/P.rubi/$Strain/"$Strain"_ORFsUniq_RxLR_motif_hmm.gff)
        GffOrfCRN=$(ls analysis/CRN_effectors/hmmer_CRN/P.rubi/$Strain/"$Strain"_ORFsUniq_CRN_hmmer.bed)
        if [ -f repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked.fa)
            echo $Assembly
        elif [ -f repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked.fa)
            echo $Assembly
        elif [ -f repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked.fa)
            echo $Assembly
        else
            Assembly=$(ls repeat_masked/quiver_results/polished/filtered_contigs_repmask/*_softmasked.fa)
            echo $Assembly
        fi
        OutDir=gene_pred/annotation/P.rubi/$Strain
        mkdir -p $OutDir
        cat $GeneGff  \
        	| sed 's/\(\tCDS\t.*\)transcript_id "\(.*\)"; gene_id.*/\1ID=\2.CDS; Parent=\2/g' \
        	| sed 's/\(\exon\t.*\)transcript_id "\(.*\)"; gene_id.*/\1ID=\2.exon; Parent=\2/g' \
        	| sed 's/transcript_id "/ID=/g' | sed 's/";/;/g' | sed 's/ gene_id "/Parent=/g' \
        	| sed -r "s/\tg/\tID=g/g" | sed 's/ID=gene/gene/g' | sed -r "s/;$//g" \
        	| sed "s/\ttranscript\t.*ID=\(.*\).t.*$/\0;Parent=\1/" \
        	> $OutDir/"$Strain"_genes_incl_ORFeffectors_noApoP.gff3
        # cat $GeneGff > $OutDir/10300_genes_incl_ORFeffectors.gff3
        ProgDir=/home/adamst/git_repos/scripts/phytophthora/10300_analysis
        $ProgDir/gff_name2id.py --gff $GffOrfRxLR > $OutDir/ORF_RxLR_parsed.gff3
        $ProgDir/gff_name2id.py --gff $GffOrfCRN > $OutDir/ORF_CRN_parsed.gff3

        ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
        $ProgDir/add_ORF_features.pl $OutDir/ORF_RxLR_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors_noApoP.gff3
        $ProgDir/add_ORF_features.pl $OutDir/ORF_CRN_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors_noApoP.gff3
        # Make gene models from gff files.
        ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
        $ProgDir/gff2fasta.pl $Assembly $OutDir/"$Strain"_genes_incl_ORFeffectors_noApoP.gff3 $OutDir/"$Strain"_genes_incl_ORFeffectors_noApoP
        # Note - these fasta files have not been validated - do not use
    done
    echo "Gff created with all RxLRs and without ApoplastP hits"

    #With EER discrimination and without ApoplastP
    for GeneGff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3)
    do
        Strain=$(echo $GeneGff | rev | cut -d '/' -f3 | rev)
        Organism=$(echo $GeneGff | rev | cut -d '/' -f4 | rev)
        echo "$Strain - $Organism"
        GffOrfRxLR=$(ls analysis/RxLR_effectors/combined_evidence/P.rubi/$Strain/"$Strain"_ORFsUniq_RxLR_EER_motif_hmm.gff)
        GffOrfCRN=$(ls analysis/CRN_effectors/hmmer_CRN/P.rubi/$Strain/"$Strain"_ORFsUniq_CRN_hmmer.bed)
        if [ -f repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/manual_edits_repmask/*_softmasked.fa)
            echo $Assembly
        elif [ -f repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked.fa)
            echo $Assembly
        elif [ -f repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked.fa ]
        then
            Assembly=$(ls repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/*_softmasked.fa)
            echo $Assembly
        else
            Assembly=$(ls repeat_masked/quiver_results/polished/filtered_contigs_repmask/*_softmasked.fa)
            echo $Assembly
        fi
        OutDir=gene_pred/annotation/P.rubi/$Strain
        mkdir -p $OutDir
        cat $GeneGff  \
        	| sed 's/\(\tCDS\t.*\)transcript_id "\(.*\)"; gene_id.*/\1ID=\2.CDS; Parent=\2/g' \
        	| sed 's/\(\exon\t.*\)transcript_id "\(.*\)"; gene_id.*/\1ID=\2.exon; Parent=\2/g' \
        	| sed 's/transcript_id "/ID=/g' | sed 's/";/;/g' | sed 's/ gene_id "/Parent=/g' \
        	| sed -r "s/\tg/\tID=g/g" | sed 's/ID=gene/gene/g' | sed -r "s/;$//g" \
        	| sed "s/\ttranscript\t.*ID=\(.*\).t.*$/\0;Parent=\1/" \
        	> $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative_noApoP.gff3
        # cat $GeneGff > $OutDir/10300_genes_incl_ORFeffectors.gff3
        ProgDir=/home/adamst/git_repos/scripts/phytophthora/10300_analysis
        $ProgDir/gff_name2id.py --gff $GffOrfRxLR > $OutDir/ORF_RxLR_EER_parsed.gff3
        $ProgDir/gff_name2id.py --gff $GffOrfCRN > $OutDir/ORF_CRN_parsed.gff3

        ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
        $ProgDir/add_ORF_features.pl $OutDir/ORF_RxLR_EER_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative_noApoP.gff3
        $ProgDir/add_ORF_features.pl $OutDir/ORF_CRN_parsed.gff3 $Assembly >> $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative_noApoP.gff3
        # Make gene models from gff files.
        ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
        $ProgDir/gff2fasta.pl $Assembly $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative_noApoP.gff3 $OutDir/"$Strain"_genes_incl_ORFeffectors_conservative_noApoP
        # Note - these fasta files have not been validated - do not use
    done
    echo "Gff created with high confidence RxLRs and without ApoplastP hits"
done
```

Additional Gff clean up performed in the phytophthora_fragariae repository in
Annotation&Short_Read_submissions.md

## Functional annotation

### A)Interproscan

Interproscan was used to give gene models functional annotations.

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/interproscan/
for Genes in $(ls gene_pred/annotation/P.rubi/*/*_genes_incl_ORFeffectors_renamed.pep.fasta)
do
    $ProgDir/sub_interproscan_TA.sh $Genes
done
```

Following this, split files were combined as follows:

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/interproscan
for Proteome in $(ls gene_pred/annotation/P.rubi/*/*_genes_incl_ORFeffectors_renamed.pep.fasta)
do
    Strain=$(echo $Proteome | rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Proteome | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    echo $Strain
    InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
    $ProgDir/append_interpro_TA.sh $Proteome $InterProRaw
done
```

### B)Swissprot

```bash
for Proteome in $(ls gene_pred/annotation/P.rubi/*/*_genes_incl_ORFeffectors_renamed.pep.fasta)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    OutDir=gene_pred/swissprot/$Organism/$Strain/greedy
    SwissDbDir=../../uniprot/swissprot
    SwissDbName=uniprot_sprot
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/swissprot
    Jobs=$(qstat | grep 'sub_swiss' | wc -l)
    while [ $Jobs -gt 5 ]
    do
        sleep 1
        printf "."
        Jobs=$(qstat | grep 'sub_swiss' | wc -l)
    done
    qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```

### C)Identify genes with transmembrane domains

WARNING: This has a high false positive rate - modified from the base script to allow greedy/conservative

```bash
for Proteome in $(ls gene_pred/annotation/P.rubi/*/*_genes_incl_ORFeffectors_renamed.pep.fasta)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/transmembrane_helices
    Jobs=$(qstat | grep 'submit_TM' | wc -l)
    while [ $Jobs -gt 5 ]
    do
        sleep 1
        printf "."
        Jobs=$(qstat | grep 'submit_TM' | wc -l)
    done
    qsub $ProgDir/submit_TMHMM.sh $Proteome
done
```

Summarise numbers of TM Proteins

```bash
for TM in $(ls gene_pred/trans_mem/P.rubi/*/*/*_TM_genes_pos.txt)
do
    Strain=$(echo $TM | rev | cut -f3 -d '/' | rev)
    echo "$Strain"
    echo "The number of proteins scoring positive for a transmembrane helix is:"
    cat $TM | wc -l
    echo ""
done
```

```
SCRP249 - conservative_noApoP
The number of proteins scoring positive for a transmembrane helix is:
4,180

SCRP249 - conservative
The number of proteins scoring positive for a transmembrane helix is:
5,427

SCRP249 - greedy_noApoP
The number of proteins scoring positive for a transmembrane helix is:
4,510

SCRP249 - greedy
The number of proteins scoring positive for a transmembrane helix is:
5,757

SCRP324 - conservative_noApoP
The number of proteins scoring positive for a transmembrane helix is:
4,272

SCRP324 - conservative
The number of proteins scoring positive for a transmembrane helix is:
5,529

SCRP324 - greedy_noApoP
The number of proteins scoring positive for a transmembrane helix is:
4,610

SCRP324 - greedy
The number of proteins scoring positive for a transmembrane helix is:
5,867

SCRP333 - conservative_noApoP
The number of proteins scoring positive for a transmembrane helix is:
4,213

SCRP333 - conservative
The number of proteins scoring positive for a transmembrane helix is:
5,442

SCRP333 - greedy_noApoP
The number of proteins scoring positive for a transmembrane helix is:
4,541

SCRP333 - greedy
The number of proteins scoring positive for a transmembrane helix is:
5,770
```

Create a headers file

```bash
for PosFile in $(ls gene_pred/trans_mem/*/*/*/*_TM_genes_pos.txt)
do
    TmHeaders=$(echo $PosFile | sed 's/.txt/_headers.txt/g')
    cat $PosFile | cut -f1 > $TmHeaders
done
```

### D)Identify genes with GPI anchors

Proteins were identified by submitting the combined protein file to webserver at http://gpi.unibe.ch

Output directory made

```bash
for Proteome in $(ls gene_pred/annotation/P.rubi/*/*_genes_incl_ORFeffectors_renamed.pep.fasta)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/GPIsom/$Organism/$Strain/greedy
    mkdir -p $OutDir
done
```

Split fasta files into files no larger than 39,000 sequences

```bash
#Greedy
for file in $(ls gene_pred/annotation/P.rubi/*/*effectors_renamed.pep.fasta)
do
    Strain=$(echo $file | rev | cut -f2 -d "/" | rev)
    OutDir=split_files/P.rubi/$Strain
    mkdir -p $OutDir
    awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%39000==0){file=sprintf("greedy%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $file
    for output in $(ls *.fa)
    do
    mv $output $OutDir/.
    done
done
```

Results were parsed to the file

```bash
nano gene_pred/GPIsom/P.rubi/SCRP249/greedy/GPI_pos.fa
nano gene_pred/GPIsom/P.rubi/SCRP324/greedy/GPI_pos.fa
nano gene_pred/GPIsom/P.rubi/SCRP333/greedy/GPI_pos.fa
```

Create a file just listing gene names

```bash
for PosFile in $(ls gene_pred/GPIsom/*/*/*/GPI_pos.fa)
do
    GPIHeaders=$(echo $PosFile | sed 's/.fa/.txt/g')
    cat $PosFile | grep -e ">" | cut -f1 -d ' ' | sed 's/>//g' > $GPIHeaders
done
```

Summarise numbers of GPI Proteins

```bash
for GPI in $(ls gene_pred/GPIsom/P.rubi/*/*/*.txt)
do
    Strain=$(echo $GPI | rev | cut -f3 -d '/' | rev)
    echo "$Strain"
    echo "The number of proteins scoring positive for being GPI anchored is:"
    cat $GPI | wc -l
    echo ""
done
```

Only greedy run due to the webserver crashing and dependencies for a local install no longer being available.

```
SCRP249 - greedy
The number of proteins scoring positive for being GPI anchored is:
2,096

SCRP324 - greedy
The number of proteins scoring positive for being GPI anchored is:
1,436

SCRP333 - greedy
The number of proteins scoring positive for being GPI anchored is:
1,367
```

## Genes with transcription factor annotations

Identify TFs & TRs in gene set
Run in a qlogin session with a screen session

```bash
# greedy ApoP
screen -a
qlogin
cd /home/groups/harrisonlab/project_files/phytophthora_rubi/
for Proteome in $(ls gene_pred/annotation/P.rubi/*/*_genes_incl_ORFeffectors_renamed.pep.fasta)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer
    HMMModel=/home/adamst/TF_TR_Buitrago-Florez.hmm
    OutDir=analysis/transcription_factors/$Organism/$Strain/greedy
    mkdir -p $OutDir
    HMMResults="$Strain"_TF_TR_hmmer.txt
    hmmsearch -T 0 $HMMModel $Proteome > $OutDir/$HMMResults
    echo "$Organism - $Strain:"
    ProgDir=/home/adamst/git_repos/tools/seq_tools/transcription_factors
    python $ProgDir/parse_hmmer_TF_TR.py --HMM_results $OutDir/$HMMResults --outdir $OutDir --Isolate $Strain
    cat $OutDir/*TF_TR_Headers.txt | wc -l
    HMMFasta="$Strain"_TF_TR_hmmer.fa
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/*TF_TR_Headers.txt > $OutDir/"$Strain"_TF_TR_hmmer.fa
    Gff=$(ls gene_pred/annotation/$Organism/$Strain/*_genes_incl_ORFeffectors_renamed.gff3)
    cat $Gff | grep -w -f $OutDir/*TF_TR_Headers.txt > $OutDir/"$Strain"_TF_TR.gff3
done
```

## Rename predicted effectors to allow further analysis with gene models on NCBI

```bash
for Isolate in SCRP249 SCRP324 SCRP333
do
    Species=P.rubi
    RxLR_IDs=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_motif_hmm.txt
    RxLR_EER_IDs=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_EER_motif_hmm.txt
    CRN_IDs=analysis/CRN_effectors/hmmer_CRN/$Species/$Isolate/"$Isolate"_final_CRN.txt
    ApoP_IDs=analysis/ApoplastP/$Species/$Isolate/"$Isolate"_Total_ApoplastP.txt
    SigP2_AA=gene_pred/final_sigP/$Species/$Isolate/"$Isolate"_aug_sp.aa
    SigP2_Headers=gene_pred/final_sigP/$Species/$Isolate/"$Isolate"_aug_sp_Headers.txt
    SigP2_ORF_AA=gene_pred/ORF_sigP/$Species/$Isolate/"$Isolate"_aug_sp.aa
    SigP2_ORF_Headers=gene_pred/ORF_sigP/$Species/$Isolate/"$Isolate"_aug_sp_Headers.txt
    SigP3_AA=gene_pred/final_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp.aa
    SigP3_Headers=gene_pred/final_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp_Headers.txt
    SigP3_ORF_AA=gene_pred/ORF_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp.aa
    SigP3_ORF_Headers=gene_pred/ORF_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp_Headers.txt
    SigP4_AA=gene_pred/final_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp.aa
    SigP4_Headers=gene_pred/final_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp_Headers.txt
    SigP4_ORF_AA=gene_pred/ORF_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp.aa
    SigP4_ORF_Headers=gene_pred/ORF_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp_Headers.txt
    Phobius_AA=analysis/phobius_CQ/$Species/$Isolate/"$Isolate"_phobius.fa
    Phobius_Headers=analysis/phobius_CQ/$Species/$Isolate/"$Isolate"_phobius_Headers.txt
    Phobius_ORF_AA=analysis/phobius_ORF/$Species/$Isolate/"$Isolate"_phobius.fa
    Phobius_ORF_Headers=analysis/phobius_ORF/$Species/$Isolate/"$Isolate"_phobius_Headers.txt
    Aug_Secreted=gene_pred/combined_sigP_CQ/$Species/$Isolate/"$Isolate"_secreted.txt
    ORF_Secreted=gene_pred/combined_sigP_ORF/$Species/$Isolate/"$Isolate"_all_secreted_merged.txt
    cat $SigP2_AA | grep '>' | cut -f1 | tr -d '>' > $SigP2_Headers
    cat $SigP2_ORF_AA | grep '>' | cut -f1 | tr -d '>' > $SigP2_ORF_Headers
    cat $SigP3_AA | grep '>' | cut -f1 | tr -d '>' > $SigP3_Headers
    cat $SigP3_ORF_AA | grep '>' | cut -f1 | tr -d '>' > $SigP3_ORF_Headers
    cat $SigP4_AA | grep '>' | cut -f1 | tr -d '>' > $SigP4_Headers
    cat $SigP4_ORF_AA | grep '>' | cut -f1 | tr -d '>' > $SigP4_ORF_Headers
    cat $Phobius_AA | grep '>' | cut -f1 | tr -d '>' > $Phobius_Headers
    cat $Phobius_ORF_AA | grep '>' | cut -f1 | tr -d '>' > $Phobius_ORF_Headers
    RxLR_IDs_renamed=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_motif_hmm_renamed.txt
    RxLR_EER_IDs_renamed=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_EER_motif_hmm_renamed.txt
    CRN_IDs_renamed=analysis/CRN_effectors/hmmer_CRN/$Species/$Isolate/"$Isolate"_final_CRN_renamed.txt
    ApoP_IDs_renamed=analysis/ApoplastP/$Species/$Isolate/"$Isolate"_Total_ApoplastP_renamed.txt
    SigP2_Headers_renamed=gene_pred/final_sigP/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP2_ORF_Headers_renamed=gene_pred/ORF_sigP/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP3_Headers_renamed=gene_pred/final_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP3_ORF_Headers_renamed=gene_pred/ORF_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP4_Headers_renamed=gene_pred/final_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP4_ORF_Headers_renamed=gene_pred/ORF_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    Phobius_Headers_renamed=analysis/phobius_CQ/$Species/$Isolate/"$Isolate"_phobius_Headers_renamed.txt
    Phobius_ORF_Headers_renamed=analysis/phobius_ORF/$Species/$Isolate/"$Isolate"_phobius_ORF_Headers_renamed.txt
    Aug_Secreted_renamed=gene_pred/combined_sigP_CQ/$Species/$Isolate/"$Isolate"_secreted_renamed.txt
    ORF_Secreted_renamed=gene_pred/combined_sigP_ORF/$Species/$Isolate/"$Isolate"_all_secreted_merged_renamed.txt
    Log_file=gene_pred/annotation/$Species/$Isolate/Renaming_log.log

    ProgDir=/home/adamst/git_repos/scripts/phytophthora_fragariae
    python $ProgDir/Rename_Effectors.py --Feature_File $RxLR_IDs --Conversion_File $Log_file --Out_File $RxLR_IDs_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $RxLR_EER_IDs --Conversion_File $Log_file --Out_File $RxLR_EER_IDs_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $CRN_IDs --Conversion_File $Log_file --Out_File $CRN_IDs_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $ApoP_IDs --Conversion_File $Log_file --Out_File $ApoP_IDs_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $SigP2_Headers --Conversion_File $Log_file --Out_File $SigP2_Headers_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $SigP2_ORF_Headers --Conversion_File $Log_file --Out_File $SigP2_ORF_Headers_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $SigP3_Headers --Conversion_File $Log_file --Out_File $SigP3_Headers_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $SigP3_ORF_Headers --Conversion_File $Log_file --Out_File $SigP3_ORF_Headers_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $SigP4_Headers --Conversion_File $Log_file --Out_File $SigP4_Headers_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $SigP4_ORF_Headers --Conversion_File $Log_file --Out_File $SigP4_ORF_Headers_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $Phobius_Headers --Conversion_File $Log_file --Out_File $Phobius_Headers_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $Phobius_ORF_Headers --Conversion_File $Log_file --Out_File $Phobius_ORF_Headers_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $Aug_Secreted --Conversion_File $Log_file --Out_File $Aug_Secreted_renamed
    python $ProgDir/Rename_Effectors.py --Feature_File $ORF_Secreted --Conversion_File $Log_file --Out_File $ORF_Secreted_renamed
done
```

### Extract fasta files with renamed gene IDs for features

```bash
for Isolate in SCRP249 SCRP324 SCRP333
do
    Species=P.rubi
    RxLR_IDs_renamed=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_motif_hmm_renamed.txt
    RxLR_EER_IDs_renamed=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_EER_motif_hmm_renamed.txt
    CRN_IDs_renamed=analysis/CRN_effectors/hmmer_CRN/$Species/$Isolate/"$Isolate"_final_CRN_renamed.txt
    ApoP_IDs_renamed=analysis/ApoplastP/$Species/$Isolate/"$Isolate"_Total_ApoplastP_renamed.txt
    SigP2_Headers_renamed=gene_pred/final_sigP/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP2_ORF_Headers_renamed=gene_pred/ORF_sigP/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP3_Headers_renamed=gene_pred/final_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP3_ORF_Headers_renamed=gene_pred/ORF_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP4_Headers_renamed=gene_pred/final_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    SigP4_ORF_Headers_renamed=gene_pred/ORF_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.txt
    Phobius_Headers_renamed=analysis/phobius_CQ/$Species/$Isolate/"$Isolate"_phobius_Headers_renamed.txt
    Phobius_ORF_Headers_renamed=analysis/phobius_ORF/$Species/$Isolate/"$Isolate"_phobius_ORF_Headers_renamed.txt
    Aug_Secreted_renamed=gene_pred/combined_sigP_CQ/$Species/$Isolate/"$Isolate"_secreted_renamed.txt
    ORF_Secreted_renamed=gene_pred/combined_sigP_ORF/$Species/$Isolate/"$Isolate"_all_secreted_merged_renamed.txt

    Proteome=gene_pred/annotation/$Species/$Isolate/"$Isolate"_genes_incl_ORFeffectors_renamed.pep.fasta

    RxLR_renamed_seqs=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_motif_hmm_renamed.fa
    RxLR_EER_renamed_seqs=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_EER_motif_hmm_renamed.fa
    CRN_renamed_seqs=analysis/CRN_effectors/hmmer_CRN/$Species/$Isolate/"$Isolate"_final_CRN_renamed.fa
    ApoP_renamed_seqs=analysis/ApoplastP/$Species/$Isolate/"$Isolate"_Total_ApoplastP_renamed.fa
    SigP2_renamed_seqs=gene_pred/final_sigP/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.fa
    SigP2_ORF_renamed_seqs=gene_pred/ORF_sigP/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.fa
    SigP3_renamed_seqs=gene_pred/final_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.fa
    SigP3_ORF_renamed_seqs=gene_pred/ORF_signalp-3.0/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.fa
    SigP4_renamed_seqs=gene_pred/final_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.fa
    SigP4_ORF_renamed_seqs=gene_pred/ORF_signalp-4.1/$Species/$Isolate/"$Isolate"_aug_sp_Headers_renamed.fa
    Phobius_renamed_seqs=analysis/phobius_CQ/$Species/$Isolate/"$Isolate"_phobius_Headers_renamed.fa
    Phobius_ORF_renamed_seqs=analysis/phobius_ORF/$Species/$Isolate/"$Isolate"_phobius_ORF_Headers_renamed.fa
    Aug_Secreted_renamed_seqs=gene_pred/combined_sigP_CQ/$Species/$Isolate/"$Isolate"_secreted_renamed.fa
    ORF_Secreted_renamed_seqs=gene_pred/combined_sigP_ORF/$Species/$Isolate/"$Isolate"_all_secreted_merged_renamed.fa

    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $RxLR_IDs_renamed > $RxLR_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $RxLR_EER_IDs_renamed > $RxLR_EER_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $CRN_IDs_renamed > $CRN_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $ApoP_IDs_renamed > $ApoP_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $SigP2_Headers_renamed > $SigP2_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $SigP2_ORF_Headers_renamed > $SigP2_ORF_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $SigP3_Headers_renamed > $SigP3_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $SigP3_ORF_Headers_renamed > $SigP3_ORF_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $SigP4_Headers_renamed > $SigP4_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $SigP4_ORF_Headers_renamed > $SigP4_ORF_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $Phobius_Headers_renamed > $Phobius_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $Phobius_ORF_Headers_renamed > $Phobius_ORF_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $Aug_Secreted_renamed > $Aug_Secreted_renamed_seqs
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $ORF_Secreted_renamed > $ORF_Secreted_renamed_seqs
done
```

#### Summarise final numbers of each class of gene

```bash
for Isolate in SCRP249 SCRP324 SCRP333
do
    Organism=P.rubi
    echo "$Organism - $Isolate"
    echo "The number of predicted RxLRs is:"
    RxLR_IDs=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_motif_hmm_renamed.txt
    cat $RxLR_IDs | wc -l
    echo "The number of predicted RxLRs with fasta sequences is:"
    RxLR_seqs=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_motif_hmm_renamed.fa
    cat $RxLR_seqs | grep '>' | wc -l
    echo "The number of high-confidence predicted RxLRs is:"
    RxLR_EER_IDs=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_EER_motif_hmm_renamed.txt
    cat $RxLR_EER_IDs | wc -l
    echo "The number of high-confidence predicted RxLRs with fasta sequences is:"
    RxLR_EER_seqs=analysis/RxLR_effectors/combined_evidence/$Species/$Isolate/"$Isolate"_Total_RxLR_EER_motif_hmm_renamed.fa
    cat $RxLR_EER_seqs | grep '>' | wc -l
    echo "The number of predicted CRN effectors is:"
    CRN_IDs=analysis/CRN_effectors/hmmer_CRN/$Species/$Isolate/"$Isolate"_final_CRN_renamed.txt
    cat $CRN_IDs | wc -l
    echo "The number of predicted CRN effectors with fasta sequences is:"
    CRN_seqs=analysis/CRN_effectors/hmmer_CRN/$Species/$Isolate/"$Isolate"_final_CRN_renamed.fa
    cat $CRN_seqs | grep '>' | wc -l
    echo "The number of predicted Apoplastic effectors is:"
    ApoP_IDs=analysis/ApoplastP/$Species/$Isolate/"$Isolate"_Total_ApoplastP_renamed.txt
    cat $ApoP_IDs | wc -l
    echo "The number of predicted Apoplastic effectors with fasta sequences is:"
    ApoP_seqs=analysis/ApoplastP/$Species/$Isolate/"$Isolate"_Total_ApoplastP_renamed.fa
    cat $ApoP_seqs | grep '>' | wc -l
    echo "The number of predicted secreted proteins is:"
    Aug_Secreted=gene_pred/combined_sigP_CQ/$Species/$Isolate/"$Isolate"_secreted_renamed.txt
    ORF_Secreted=gene_pred/combined_sigP_ORF/$Species/$Isolate/"$Isolate"_all_secreted_merged_renamed.txt
    All_Secreted=gene_pred/combined_sigP_CQ/$Species/$Isolate/"$Isolate"_all_secreted_Aug_ORF.txt
    cat $Aug_Secreted $ORF_Secreted > $All_Secreted
    cat $All_Secreted | wc -l
    echo "The number of predicted secreted proteins with fasta sequences is:"
    Proteome=gene_pred/annotation/$Species/$Isolate/"$Isolate"_genes_incl_ORFeffectors_renamed.pep.fasta
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    All_Secreted_seqs=gene_pred/combined_sigP_CQ/$Species/$Isolate/"$Isolate"_all_secreted_Aug_ORF.fa
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $All_Secreted > $All_Secreted_seqs
    cat $All_Secreted_seqs | grep '>' | wc -l
done
```

Further downstream analysis done in OrthoFinder.md and SNP_analysis.md in
phytophthora_fragariae repository
