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

```bash
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

### QUAST used to summarise assembly statistics

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Strain in SCRP249 SCRP324 SCRP333
do
    for Assembly in $(ls assembly/spades/P.rubi/$Strain/deconseq_Paen/*_500bp_filtered_renamed.fasta)
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

##Discovar

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
SCRP333: 10260 **

#Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assemblies were used to perform repeatmasking

for discovar assembly:

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
```

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

The number of bases masked by transposonPSI and Repeatmasker were summarised using the following commands:

```bash
for Assembler in discovar spades
do
    for RepDir in $(ls -d repeat_masked/$Assembler/P.*/*/filtered_contigs_repmask)
    do
        Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
        Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
        RepMaskGff=$(ls $RepDir/"$Strain"_contigs_hardmasked.gff)
        TransPSIGff=$(ls $RepDir/"$Strain"_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
        printf "$Assembler\n"
        printf "$Organism\t$Strain\n"
        printf "The number of bases masked by RepeatMasker:\t"
        sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
        printf "The number of bases masked by TransposonPSI:\t"
        sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
        printf "The total number of masked bases are:\t"
        cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    done
done
```

** discovar
P.rubi	SCRP249
The number of bases masked by RepeatMasker:	31215924
The number of bases masked by TransposonPSI:	7577505
The total number of masked bases are:	33062283
discovar
P.rubi	SCRP324
The number of bases masked by RepeatMasker:	31949318
The number of bases masked by TransposonPSI:	7607213
The total number of masked bases are:	33859484
discovar
P.rubi	SCRP333
The number of bases masked by RepeatMasker:	31597333
The number of bases masked by TransposonPSI:	7644529
The total number of masked bases are:	33426973
spades
P.rubi	SCRP249
The number of bases masked by RepeatMasker:	23482028
The number of bases masked by TransposonPSI:	5953026
The total number of masked bases are:	25419316
spades
P.rubi	SCRP324
The number of bases masked by RepeatMasker:	23613201
The number of bases masked by TransposonPSI:	5940852
The total number of masked bases are:	25417379
spades
P.rubi	SCRP333
The number of bases masked by RepeatMasker:	23158566
The number of bases masked by TransposonPSI:	5961557
The total number of masked bases are:	25028389 **

#Merging RepeatMasker and TransposonPSI outputs

```bash
for File in $(ls -d repeat_masked/*/P.*/*/filtered_contigs_repmask/*_contigs_softmasked.fa)
do
    OutDir=$(dirname $File)
    TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    Assembler=$(echo $File | rev | cut -f5 -d '/' | rev)
    OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
    bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
    echo $Assembler
    echo "$OutFile"
    echo "Number of masked bases:"
    cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
```

#Gene Prediction
Gene prediction followed three steps: Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified. Gene model training - Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline Gene prediction - Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

##Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
for Assembly in $(ls assembly/spades/P.rubi/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta)
do
    Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    echo "$Strain"
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/busco
    BuscoDB=Eukaryotic
    OutDir=assembly/spades/P.rubi/$Strain/deconseq_Paen/
    qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
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
for Assembly in $(ls repeat_masked/P.rubi/*/deconseq_Paen_repmask/"$Strain"_contigs_softmasked_repeatmasker_TPSI_appended.fa)
do
    Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for FileF in $(ls qc_rna/qc_rna/raw_rna/consortium/P.rubi/F/*_trim.fq.gz)
    do
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        while [ $Jobs -gt 1 ]
        do
            sleep 1m
            printf "."
            Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
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
for Assembly in $(ls repeat_masked/spades/P.rubi/*/deconseq_Paen_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
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
    mkdir -p alignment/$Assembler/$Organism/$Strain/concatenated
    samtools merge -f alignment/star/$Organism/$Strain/concatenated/concatenated.bam \
    alignment/star/$Organism/$Strain/4671V8/accepted_hits.bam \
    alignment/star/$Organism/$Strain/Pr4671PB/accepted_hits.bam \
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
for Assembly in $(ls repeat_masked/*/*/deconseq_Paen_repmask/*_contigs_unmasked.fa)
do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
    mkdir -p $OutDir
    AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
    ProgDir=/home/adamst/git_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Secondly, genes were predicted using CodingQuarry:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquarry/$Organism/$Strain
    CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```

Then, additional transcripts were added to Braker1 gene models, when CodingQuarry genes were predicted in regions of the genome, not containing Braker1 gene models:

```bash
for BrakerGff in $(ls gene_pred/braker/P.*/*_braker/*/augustus.gff3)
do
    Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g' | sed 's/_braker_pacbio//g')
    Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Assembly=$(ls repeat_masked/$Organism/$Strain/deconseq_Paen_repmask/"$Strain"_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    CodingQuarryGff=gene_pred/codingquarry/$Organism/$Strain/out/PredictedPass.gff3
    PGNGff=gene_pred/codingquarry/$Organism/$Strain/out/PGN_predictedPass.gff3
    AddDir=gene_pred/codingquarry/$Organism/$Strain/additional
    FinalDir=gene_pred/codingquarry/$Organism/$Strain/final
    AddGenesList=$AddDir/additional_genes.txt
    AddGenesGff=$AddDir/additional_genes.gff
    FinalGff=$AddDir/combined_genes.gff
    mkdir -p $AddDir
    mkdir -p $FinalDir

    bedtools intersect -v -a $CodingQuarryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
    bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuarryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
    $ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    # GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/additional/additional_genes.gff
    # GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/out/PredictedPass.gff3

    $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuarry.gff3
    $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuarry_unspliced.gff3
    $ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuarry_unspliced.gff3 > $FinalDir/final_genes_CodingQuarry.gff3
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

The final number of genes per isolate was observed using:

```bash
for DirPath in $(ls -d gene_pred/codingquarry/P.*/*/final)
do
    echo $DirPath
    echo Braker:
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l
    echo CodingQuarry:
    cat $DirPath/final_genes_CodingQuarry.pep.fasta | grep '>' | wc -l
    echo Total:
    cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l
    echo ""
done
```

```
discovar - P.rubi - SCRP249
Braker:
33629
CodingQuarry:
2490
Total:
36119

discovar - P.rubi - SCRP324
Braker:
39285
CodingQuarry:
2909
Total:
42194

discovar - P.rubi - SCRP333
Braker:
34276
CodingQuarry:
3236
Total:
37512

spades - P.rubi - SCRP249
Braker:
29919
CodingQuarry:
1947
Total:
31866

spades - P.rubi - SCRP324
Braker:
34668
CodingQuarry:
2353
Total:
37021

spades - P.rubi - SCRP333
Braker:
29669
CodingQuarry:
2069
Total:
31738
```

#Gene prediction 2 - atg.pl prediction of ORFs

Open reading frame predictions were made using the atg.pl script as part of the path_pipe.sh pipeline. This pipeline also identifies open reading frames containing Signal peptide sequences and RxLRs. This pipeline was run with the following commands:


```bash
ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
for Assembler in discovar spades
do
    for Genome in $(ls repeat_masked/$Assembler/P.*/*/filtered_contigs_repmask/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
    do
        echo "$Genome"
        qsub $ProgDir/run_ORF_finder2.sh $Genome
    done
done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were corrected using the following commands:

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
for OrfGff in $(ls gene_pred/ORF_finder/*/P.*/*/*_ORF.gff | grep -v 'atg')
do
    echo "$OrfGff"
    OrfGffMod=$(echo $OrfGff | sed 's/.gff/.gff3/g')
    $ProgDir/gff_corrector.pl $OrfGff > $OrfGffMod
done
```

The final number of genes per isolate were determined using

```bash
for Assembler in discovar spades
do
    echo $Assembler
    for Strain in SCRP249 SCRP324 SCRP333
    do
        for DirPath in $(ls -d gene_pred/ORF_finder/$Assembler/*/$Strain)
        do
            echo $DirPath
            cat $DirPath/"$Strain".aa_cat.fa | grep '>' | wc -l
            echo ""
        done
    done
done
```

```
discovar

SCRP249
716788

SCRP324
795934

SCRP333
720807

spades

SCRP249
645937

SCRP324
725179

SCRP333
646243
```

#Genomic analysis

##RxLR genes

Putative RxLR genes were identified within Augustus gene models using a number of approaches:

A) From Augustus gene models - Signal peptide & RxLR motif
B) From Augustus gene models - Hmm evidence of WY domains
C) From Augustus gene models - Hmm evidence of RxLR effectors
D) From Augustus gene models - Hmm evidence of CRN effectors
E) From ORF fragments - Signal peptide & RxLR motif
F) From ORF fragments - Hmm evidence of WY domains
G) From ORF fragments - Hmm evidence of RxLR effectors

##A) From Augustus gene models - Signal peptide & RxLR motif

Required programs:

SigP
biopython

####A.1) Signal peptide prediction using SignalP 2.0

Proteins that were predicted to contain signal peptides were identified using the following commands:

```bash
for Assembler in discovar spades
do
    for Strain in SCRP249 SCRP324 SCRP333
    do
        for Proteome in $(ls gene_pred/braker/$Assembler/*/"$Strain"_braker/*/augustus.aa)
        do
            SplitfileDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
            ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
            Organism=P.rubi
            SplitDir=gene_pred/braker_split/$Assembler/$Organism/$Strain
            mkdir -p $SplitDir
            BaseName="$Organism""_$Strain"_braker
            $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
            for File in $(ls $SplitDir/*_braker_*)
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
                qsub $ProgDir/pred_sigP.sh $File
                qsub $ProgDir/pred_sigP.sh $File signalp-4.1
            done
        done
    done
done
```

The batch files of predicted secreted proteins needed to be combined into a single file for each strain. This was done with the following commands:

```bash
for Assembler in discovar spades
do
    for Strain in SCRP249 SCRP324 SCRP333
    do
        for SplitDir in $(ls -d gene_pred/braker_split/$Assembler/P.*/$Strain)
        do
            Organism=P.fragariae
            echo "$Organism - $Strain"
            InStringAA=''
            InStringNeg=''
            InStringTab=''
            InStringTxt=''
            for SigpDir in $(ls -d gene_pred/"$Assembler"_sig* | cut -f2 -d'/')
            do
                for GRP in $(ls -l $SplitDir/*_braker_*.fa | rev | cut -d '_' -f1 | rev | sort -n)
                do  
                    InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_$GRP""_sp.aa"
                    InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_$GRP""_sp_neg.aa"
                    InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_$GRP""_sp.tab"
                    InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_$GRP""_sp.txt"
                done
                cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
                cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
                tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
                cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
            done
        done
    done
done
```

####B.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
for Assembler in discovar spades
do
    for Strain in SCRP249 SCRP324 SCRP333
    do
        for Proteome in $(ls gene_pred/braker/$Assembler/*/"$Strain"_braker/*/augustus.aa)
        do
            Organism=P.rubi
            echo "$Organism - $Strain"
            OutDir=analysis/phobius/$Assembler/$Organism/$Strain
            mkdir -p $OutDir
            phobius.pl $Proteome > $OutDir/"$Strain"_phobius.txt
            ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
            $ProgDir/phobius_parser.py --inp_fasta $Proteome --phobius_txt $OutDir/"$Strain"_phobius.txt --out_fasta $OutDir/"$Strain"_phobius.fa
        done
    done
done
```

Secreted proteins from different sources were combined into a single file:

```bash
for Assembler in discovar spades
do
    echo $Assembler
    for Strain in SCRP249 SCRP324 SCRP333
    do
        for Proteome in $(ls gene_pred/braker/$Assembler/*/"$Strain"_braker/*/augustus.aa)
        do
            Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
            echo "$Organism - $Strain"
            OutDir=gene_pred/combined_sigP/$Assembler/$Organism/$Strain
            mkdir -p $OutDir
            echo "The following number of sequences were predicted as secreted:"
            cat gene_pred/braker_sig*/$Assembler/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Assembler/$Organism/$Strain/"$Strain"_phobius.fa > $OutDir/"$Strain"_all_secreted.fa
            cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | wc -l
            echo "This represented the following number of unique genes:"
            cat gene_pred/braker_sig*/$Assembler/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Assembler/$Organism/$Strain/"$Strain"_phobius.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
            ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
            $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
            cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
        done
    done
done
```

```
discovar

P.rubi - SCRP249
The following number of sequences were predicted as secreted:
10696
This represented the following number of unique genes:
3748
P.rubi - SCRP324
The following number of sequences were predicted as secreted:
12476
This represented the following number of unique genes:
4434
P.rubi - SCRP333
The following number of sequences were predicted as secreted:
10397
This represented the following number of unique genes:
3714

spades

P.rubi - SCRP249
The following number of sequences were predicted as secreted:
10174
This represented the following number of unique genes:
3546
P.rubi - SCRP324
The following number of sequences were predicted as secreted:
12309
This represented the following number of unique genes:
4347
P.rubi - SCRP333
The following number of sequences were predicted as secreted:
9722
This represented the following number of unique genes:
3451
```

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP/*/*/*/*_all_secreted.fa)
do
    Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev)
    Assembler=$(echo $Secretome | rev |  cut -d '/' -f4 | rev)
    Proteome=$(ls gene_pred/braker/"$Assembler"/$Organism/"$Strain"_braker/*/augustus.aa)
    Gff=$(ls gene_pred/braker/"$Assembler"/$Organism/"$Strain"_braker/*/augustus_extracted.gff)
    OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Assembler"/"$Organism"/"$Strain"
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
discovar

strain: SCRP249	species: P.rubi
the total number of SigP gene is:	10696
the number of unique SigP gene is:	3748
the number of SigP-RxLR genes are:	379
the number of SigP-RxLR-EER genes are:	200

strain: SCRP324	species: P.rubi
the total number of SigP gene is:	12476
the number of unique SigP gene is:	4434
the number of SigP-RxLR genes are:	377
the number of SigP-RxLR-EER genes are:	192

strain: SCRP333	species: P.rubi
the total number of SigP gene is:	10397
the number of unique SigP gene is:	3714
the number of SigP-RxLR genes are:	354
the number of SigP-RxLR-EER genes are:	179

spades

strain: SCRP249	species: P.rubi
the total number of SigP gene is:	10174
the number of unique SigP gene is:	3546
the number of SigP-RxLR genes are:	346
the number of SigP-RxLR-EER genes are:	176

strain: SCRP324	species: P.rubi
the total number of SigP gene is:	12309
the number of unique SigP gene is:	4347
the number of SigP-RxLR genes are:	378
the number of SigP-RxLR-EER genes are:	191

strain: SCRP333	species: P.rubi
the total number of SigP gene is:	9722
the number of unique SigP gene is:	3451
the number of SigP-RxLR genes are:	341
the number of SigP-RxLR-EER genes are:	171
```

####G) From Secreted gene models - Hmm evidence of RxLR effectors

```bash
for Assembler in discovar spades
do
    echo $Assembler
    for Strain in SCRP249 SCRP324 SCRP333
    do
        for Proteome in $(ls gene_pred/braker/$Assembler/*/"$Strain"_braker/*/augustus.aa)
        do
            ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer
            HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
            Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
            OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Assembler/$Organism/$Strain
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
            Gff=$(ls gene_pred/braker/$Assembler/$Organism/"$Strain"_braker/*/augustus_extracted.gff)
            cat $Gff | grep -w -f $OutDir/$Headers > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
        done
    done
done
```

```
discovar

P.rubi SCRP249
Initial search space (Z):              36943  [actual number of targets]
Domain search space  (domZ):             202  [number of targets reported over threshold]
P.rubi SCRP324
Initial search space (Z):              42604  [actual number of targets]
Domain search space  (domZ):             207  [number of targets reported over threshold]
P.rubi SCRP333
Initial search space (Z):              36843  [actual number of targets]
Domain search space  (domZ):             189  [number of targets reported over threshold]

spades

P.rubi SCRP249
Initial search space (Z):              32541  [actual number of targets]
Domain search space  (domZ):             195  [number of targets reported over threshold]
P.rubi SCRP324
Initial search space (Z):              38842  [actual number of targets]
Domain search space  (domZ):             202  [number of targets reported over threshold]
P.rubi SCRP333
Initial search space (Z):              32562  [actual number of targets]
Domain search space  (domZ):             188  [number of targets reported over threshold]
```

####F) Combining RxLRs from Regex and hmm searches

```bash
for Assembler in discovar spades
do
    echo $Assembler
    for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/$Assembler/P.rubi/*/*_RxLR_EER_regex.txt)
    do
        Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
        Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
        Gff=$(ls gene_pred/braker/$Assembler/$Organism/"$Strain"_braker/*/augustus_extracted.gff)
        Proteome=$(ls gene_pred/braker/$Assembler/$Organism/"$Strain"_braker/*/augustus.aa)
        HmmRxLR=analysis/RxLR_effectors/hmmer_RxLR/$Assembler/$Organism/$Strain/*_RxLR_hmmer_headers.txt
        echo "$Organism - $Strain"
        echo "Number of RxLRs identified by Regex:"
        cat $RegexRxLR | sort | uniq | wc -l
        echo "Number of RxLRs identified by Hmm:"
        cat $HmmRxLR | sort | uniq | wc -l
        echo "Number of RxLRs in combined dataset:"
        cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l
        # echo "Number of RxLRs in both datasets:"
        # cat $RegexRxLR $HmmRxLR | sort | uniq -d | wc -l
        # echo "Extracting RxLRs from datasets"
        OutDir=analysis/RxLR_effectors/combined_evidence/$Assembler/$Organism/$Strain
        mkdir -p $OutDir
        cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_RxLR_headers.txt
        cat $Gff | grep -w -f $OutDir/"$Strain"_total_RxLR_headers.txt > $OutDir/"$Strain"_total_RxLR.gff
        echo "Number of genes in the extracted gff file:"
        cat $OutDir/"$Strain"_total_RxLR.gff | grep -w 'gene' | wc -l
        echo ""
    done
done
```

```
discovar

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

spades

P.rubi - SCRP249
Number of RxLRs identified by Regex:
176
Number of RxLRs identified by Hmm:
195
Number of RxLRs in combined dataset:
231
Number of genes in the extracted gff file:
231

P.rubi - SCRP324
Number of RxLRs identified by Regex:
191
Number of RxLRs identified by Hmm:
202
Number of RxLRs in combined dataset:
242
Number of genes in the extracted gff file:
242

P.rubi - SCRP333
Number of RxLRs identified by Regex:
171
Number of RxLRs identified by Hmm:
188
Number of RxLRs in combined dataset:
224
Number of genes in the extracted gff file:
224
```

####D) From Augustus gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers in Augustus gene models. This was done with the following commands:

```bash
HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
LFLAK_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm)
DWL_hmm=$(ls $HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm)
for Assembler in discovar spades
do
    echo $Assembler
    for Strain in SCRP249 SCRP324 SCRP333
    do
        for Proteome in $(ls gene_pred/braker/$Assembler/P.rubi/"$Strain"_braker/*/augustus.aa)
        do
            Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
            OutDir=analysis/CRN_effectors/hmmer_CRN/$Assembler/$Organism/$Strain
            mkdir -p $OutDir
            echo "$Organism - $Strain"
            # Run hmm searches LFLAK domains
            CrinklerProts_LFLAK=$OutDir/"$Strain"_pub_CRN_LFLAK_hmm.txt
            hmmsearch -T0 $LFLAK_hmm $Proteome > $CrinklerProts_LFLAK
            cat $CrinklerProts_LFLAK | grep 'Initial search space'
            cat $CrinklerProts_LFLAK | grep 'number of targets reported over threshold'
            ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
            $ProgDir/hmmer2fasta.pl $CrinklerProts_LFLAK $Proteome > $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa
            # Run hmm searches DWL domains
            CrinklerProts_DWL=$OutDir/"$Strain"_pub_CRN_DWL_hmm.txt
            hmmsearch -T0 $DWL_hmm $Proteome > $CrinklerProts_DWL
            cat $CrinklerProts_DWL | grep 'Initial search space'
            cat $CrinklerProts_DWL | grep 'number of targets reported over threshold'
            ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
            $ProgDir/hmmer2fasta.pl $CrinklerProts_DWL $Proteome > $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa
            echo "Identify the genes detected in both models:"
            cat $OutDir/"$Strain"_pub_CRN_LFLAK_hmm.fa $OutDir/"$Strain"_pub_CRN_DWL_hmm.fa | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt
            cat $OutDir/"$Strain"_pub_CRN_LFLAK_DWL.txt | wc -l
        done
    done
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

spades

P.rubi - SCRP249
Initial search space (Z):              32541  [actual number of targets]
Domain search space  (domZ):             149  [number of targets reported over threshold]
Initial search space (Z):              32541  [actual number of targets]
Domain search space  (domZ):             137  [number of targets reported over threshold]
Identify the genes detected in both models:     126

P.rubi - SCRP324
Initial search space (Z):              38842  [actual number of targets]
Domain search space  (domZ):             154  [number of targets reported over threshold]
Initial search space (Z):              38842  [actual number of targets]
Domain search space  (domZ):             140  [number of targets reported over threshold]
Identify the genes detected in both models:     124

P.rubi - SCRP333
Initial search space (Z):              32562  [actual number of targets]
Domain search space  (domZ):             141  [number of targets reported over threshold]
Initial search space (Z):              32562  [actual number of targets]
Domain search space  (domZ):             136  [number of targets reported over threshold]
Identify the genes detected in both models:     121
```

Extract gff annotations for Crinklers:

```bash
for Assembler in discovar spades
do
    for CRNlist in $(ls analysis/CRN_effectors/hmmer_CRN/$Assembler/P.rubi/*/*_pub_CRN_LFLAK_DWL.txt)
    do
        Strain=$(echo $CRNlist | rev | cut -f2 -d '/' | rev)
        Organism=$(echo $CRNlist | rev | cut -f3 -d '/' | rev)
        OutName=$(echo $CRNlist | sed 's/.txt/.gff/g')
        echo "$Organism - $Strain"
        Gff=$(ls gene_pred/braker/$Assembler/$Organism/"$Strain"_braker/*/augustus_extracted.gff)
        cat $CRNlist | sed -r 's/\.t.$//g' > tmp.txt
        cat $Gff | grep -w -f tmp.txt > $OutName
        rm tmp.txt
    done
done
```

####E) From ORF gene models - Signal peptide & RxLR motif

Required programs:

SigP
Phobius
biopython

#####E.1) Prediction using SignalP

Proteins that were predicted to contain signal peptides were identified using the following commands:

```bash
for Assembler in discovar spades
do
    for Proteome in $(ls gene_pred/ORF_finder/$Assembler/P.*/*/*.aa_cat.fa)
    do
        SplitfileDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
        ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation/signal_peptides
        Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
        Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
        echo "$Organism - $Strain"
        SplitDir=gene_pred/ORF_split/$Assembler/$Organism/$Strain
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
            qsub $ProgDir/pred_sigP.sh $File
            qsub $ProgDir/pred_sigP.sh $File signalp-4.1
        done
    done
done
```

The batch files of predicted secreted proteins needed to be combined into a single file for each strain. This was done with the following commands:

```bash
for Assembler in discovar spades
do
    for SplitDir in $(ls -d gene_pred/ORF_split/$Assembler/*/*)
    do
        Strain=$(echo $SplitDir | cut -d '/' -f5)
        Organism=$(echo $SplitDir | cut -d '/' -f4)
        echo "$Organism - $Strain"
        InStringAA=''
        InStringNeg=''
        InStringTab=''
        InStringTxt=''
        for SigpDir in $(ls -d gene_pred/"$Assembler"_sig* | cut -f2 -d'/')
        do
            for GRP in $(ls -l $SplitDir/*_ORF_*.fa | rev | cut -d '_' -f1 | rev | sort -n)
            do
                InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.aa"
                InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp_neg.aa"
                InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.tab"
                InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_preds_$GRP""_sp.txt"
            done
            SigPDir=gene_pred/ORF_sig*/$Assembler/$Organism/$Strain
            cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
            cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
            tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
            cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
        done
    done
done
```

E.2) Prediction using Phobius

Secreted proteins were also predicted using Phobius

```bash
for Proteome in $(ls gene_pred/ORF_finder/*/P.*/*/*.aa_cat.fa)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Assembler=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/phobius/$Assembler/$Organism/$Strain
    mkdir -p $OutDir
    phobius.pl $Proteome > $OutDir/"$Strain"_phobius_ORF.txt
    cat $OutDir/"$Strain"_phobius_ORF.txt | grep -B1 'SIGNAL' | grep 'ID' | sed s'/ID.*g/g/g' > $OutDir/"$Strain"_phobius_headers_ORF.txt
done
```

Because of the way ORF_finder predicts proteins, phobius predictions cannot be used downstream as there is no way to remove overlapping features.

Secreted proteins from different sources were combined into a single file:

```bash
for Proteome in $(ls gene_pred/ORF_finder/*/P.*/*/*.aa_cat.fa)
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
    cat $OutDir/"$Strain"_all_secreted.fa | grep '>' | tr -d '>' | tr -d ' ' | sed "s/HMM_score\t/HMM_score=\t/g" > $OutDir/"$Strain"_all_secreted_headers.txt
    cat $OutDir/"$Strain"_all_secreted_headers.txt | wc -l
    echo "This represented the following number of unique genes:"
    # cat gene_pred/"$Assembler"_sig*/$Organism/$Strain/*_aug_sp.aa analysis/phobius/$Assembler/$Organism/$Strain/"$Strain"_phobius_ORF.fa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    cat gene_pred/"$Assembler"_sig*/$Organism/$Strain/*_aug_sp.aa | grep '>' | cut -f1 | tr -d ' >' | sort -g | uniq > $OutDir/"$Strain"_secreted.txt
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Proteome --headers $OutDir/"$Strain"_secreted.txt > $OutDir/"$Strain"_secreted.fa
    cat $OutDir/"$Strain"_secreted.fa | grep '>' | wc -l
done
```

For spades assembly of SCRP333, the production of all_secreted_headers.txt required the removal of one tr step shown below:

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
```

```
discovar - P.rubi - SCRP249
The following number of sequences were predicted as secreted:
81145
This represented the following number of unique genes:
38987
discovar - P.rubi - SCRP324
The following number of sequences were predicted as secreted:
89423
This represented the following number of unique genes:
42890
discovar - P.rubi - SCRP333
The following number of sequences were predicted as secreted:
81434
This represented the following number of unique genes:
39145
spades - P.rubi - SCRP249
The following number of sequences were predicted as secreted:
74633
This represented the following number of unique genes:
35754
spades - P.rubi - SCRP324
The following number of sequences were predicted as secreted:
83177
This represented the following number of unique genes:
39827
spades - P.rubi - SCRP333
The following number of sequences were predicted as secreted:
74477
This represented the following number of unique genes:
35686
```

E.3) Prediction of RxLRs

Names of ORFs containing signal peptides were extracted from fasta files. This included information on the position and hmm score of RxLRs.

```bash
for FastaFile in $(ls gene_pred/combined_sigP_ORF/*/*/*/*_all_secreted.fa)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Assembler=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SigP_headers=gene_pred/combined_sigP_ORF/$Assembler/$Organism/$Strain/"$Strain"_all_secreted_headers.txt
    cat $FastaFile | grep '>' | sed -r 's/>//g' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | sed 's/--//g' > $SigP_headers
done
```

Due to errors in spades SCRP333 fasta file, script was modified to look the same as those produced above

```bash
for FastaFile in $(ls gene_pred/combined_sigP_ORF/spades/*/SCRP333/SCRP333_all_secreted.fa)
do
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Assembler=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SigP_headers=gene_pred/combined_sigP_ORF/$Assembler/$Organism/$Strain/"$Strain"_all_secreted_headers.txt
    cat $FastaFile | grep '>' | sed -r 's/>//g' | sed "s/HMM_score\t/HMM_score=\t/g" | sed -r 's/\s+/\t/g' > $SigP_headers
done
```

Due to the nature of predicting ORFs, some features overlapped with one another. A single ORF was selected from each set of overlapped ORFs. This was was selected on the basis of its SignalP Hmm score. Biopython was used to identify overlaps and identify the ORF with the best signalP score.

```bash
for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*/*_ORF.gff3)
do
    Organism=$(echo $ORF_Gff | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $ORF_Gff | rev | cut -d '/' -f2 | rev)
    Assembler=$(echo $ORF_Gff | rev | cut -d '/' -f4 | rev)
    OutDir=$(ls -d gene_pred/combined_sigP_ORF/$Assembler/$Organism/$Strain)
    echo "$Assembler - $Organism - $Strain"
    SigP_headers=$(ls $OutDir/"$Strain"_all_secreted_headers.txt)
    SigP_fasta=$(ls $OutDir/"$Strain"_all_secreted.fa)
    ORF_fasta=$(ls gene_pred/ORF_finder/$Assembler/$Organism/$Strain/"$Strain".aa_cat.fa)

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
    echo "$Assembler $Strain complete"
done
```

The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identify RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*/*_all_secreted.fa)
do
    ProgDir=/home/adamst/git_repos/tools/pathogen/RxLR_effectors
    Strain=$(echo $Secretome | rev | cut -d '/' -f2 | rev)
    Organism=$(echo $Secretome | rev |  cut -d '/' -f3 | rev)
    Assembler=$(echo $Secretome | rev | cut -d '/' -f4 | rev)
    OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Assembler"/"$Organism"/"$Strain"
    mkdir -p $OutDir
    printf "\nassembler: $Assembler\tstrain: $Strain\tspecies: $Organism\n" >> report.txt
    printf "the number of SigP genes is:\t" >> report.txt
    cat $Secretome | grep '>' | wc -l >> report.txt
    printf "the number of SigP-RxLR genes are:\t" >> report.txt
    $ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa
    cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt
    cat $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l >> report.txt
    printf "the number of SigP-RxLR-EER genes are:\t" >> report.txt
    cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.fa | grep '>' | grep 'EER_motif_start' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' '> $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt
    cat $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt | tr -d ' ' | sort | uniq | wc -l >> report.txt
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    SigP_Gff=gene_pred/combined_sigP_ORF/$Assembler/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
    ORF_fasta=$(ls gene_pred/ORF_finder/$Assembler/$Organism/$Strain/"$Strain".aa_cat.fa)
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_ORF_RxLR_EER_regex_unmerged.txt  $SigP_Gff   RxLR_EER_regex_finder.py Name Augustus > $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff
    RxLR_Merged_Gff=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.gff
    RxLR_Merged_txt=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.txt
    RxLR_Merged_AA=$OutDir/"$Strain"_ORF_RxLR_EER_regex_merged.aa
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $OutDir/"$Strain"_ORF_RxLR_regex_unmerged.gff --db sigP_ORF_RxLR.db
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp sigP_ORF_RxLR.db --id sigP_ORF_RxLR --out sigP_ORF_RxLR_merged.db --gff > $RxLR_Merged_Gff
    cat $RxLR_Merged_Gff | grep 'transcript' | rev | cut -f1 -d '=' | rev > $RxLR_Merged_txt
    $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
    printf "Merged RxLR-EER regex proteins:\t" >> report.txt
    cat $RxLR_Merged_AA | grep '>' | wc -l >> report.txt
    printf "\n" >> report.txt
    echo "$Assembler $Strain done"
done
```

```
assembler: discovar     strain: SCRP249 species: P.rubi
the number of SigP genes is:    81145
the number of SigP-RxLR genes are:      2481
the number of SigP-RxLR-EER genes are:  311
Merged RxLR-EER regex proteins: 74

assembler: discovar     strain: SCRP324 species: P.rubi
the number of SigP genes is:    89423
the number of SigP-RxLR genes are:      2583
the number of SigP-RxLR-EER genes are:  304
Merged RxLR-EER regex proteins: 276

assembler: discovar     strain: SCRP333 species: P.rubi
the number of SigP genes is:    81434
the number of SigP-RxLR genes are:      2471
the number of SigP-RxLR-EER genes are:  298
Merged RxLR-EER regex proteins: 270

assembler: spades       strain: SCRP249 species: P.rubi
the number of SigP genes is:    74633
the number of SigP-RxLR genes are:      2303
the number of SigP-RxLR-EER genes are:  300
Merged RxLR-EER regex proteins: 270

assembler: spades       strain: SCRP324 species: P.rubi
the number of SigP genes is:    83177
the number of SigP-RxLR genes are:      2402
the number of SigP-RxLR-EER genes are:  294
Merged RxLR-EER regex proteins: 266

assembler: spades       strain: SCRP333 species: P.rubi
the number of SigP genes is:    74477
the number of SigP-RxLR genes are:      2279
the number of SigP-RxLR-EER genes are:  290
Merged RxLR-EER regex proteins: 264
```

E5) From ORF gene models - Hmm evidence of WY domains

Hmm models for the WY domain contained in many RxLRs were used to search ORFs predicted with atg.pl. These were run with the following commands:

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*/*_all_secreted.fa)
do
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
    Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
    Assembler=$(echo $Secretome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/RxLR_effectors/hmmer_WY/$Assembler/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_ORF_WY_hmmer.txt
    hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
    echo "$Assembler $Organism $Strain" >> report.txt
    cat $OutDir/$HmmResults | grep 'Initial search space' >> report.txt
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold' >> report.txt
    HmmFasta="$Strain"_ORF_WY_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
    Headers="$Strain"_ORF_WY_hmmer_headers.txt
    cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
    SigP_Merged_Gff=gene_pred/combined_sigP_ORF/$Assembler/$Organism/$Strain/"$Strain"_all_secreted_merged.gff
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/$Headers $SigP_Merged_Gff $HmmModel Name Augustus > $OutDir/"$Strain"_ORF_WY_hmmer.gff
    echo "$Assembler $Strain done"
done
```

```
discovar P.rubi SCRP249
Initial search space (Z):              81145  [actual number of targets]
Domain search space  (domZ):             416  [number of targets reported over threshold]
discovar P.rubi SCRP324
Initial search space (Z):              89423  [actual number of targets]
Domain search space  (domZ):             419  [number of targets reported over threshold]
discovar P.rubi SCRP333
Initial search space (Z):              81434  [actual number of targets]
Domain search space  (domZ):             387  [number of targets reported over threshold]
spades P.rubi SCRP249
Initial search space (Z):              74633  [actual number of targets]
Domain search space  (domZ):             404  [number of targets reported over threshold]
spades P.rubi SCRP324
Initial search space (Z):              83177  [actual number of targets]
Domain search space  (domZ):             413  [number of targets reported over threshold]
spades P.rubi SCRP333
Initial search space (Z):              74477  [actual number of targets]
Domain search space  (domZ):             390  [number of targets reported over threshold]
```

E6) From ORF gene models - Hmm evidence of RxLR effectors

```bash
for Secretome in $(ls gene_pred/combined_sigP_ORF/*/*/*/*_all_secreted.fa)
do
    ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
    HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
    Strain=$(echo $Secretome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
    Assembler=$(echo $Secretome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Assembler/$Organism/$Strain
    mkdir -p $OutDir
    HmmResults="$Strain"_ORF_RxLR_hmmer_unmerged.txt
    hmmsearch -T 0 $HmmModel $Secretome > $OutDir/$HmmResults
    echo "$Assembler $Organism $Strain" >> report.txt
    cat $OutDir/$HmmResults | grep 'Initial search space' >> report.txt
    cat $OutDir/$HmmResults | grep 'number of targets reported over threshold' >> report.txt
    HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Secretome > $OutDir/$HmmFasta
    Headers="$Strain"_ORF_RxLR_hmmer_headers_unmerged.txt
    cat $OutDir/$HmmFasta | grep '>' | cut -f1 | tr -d '>' | sed -r 's/\.t.*//' | tr -d ' ' > $OutDir/$Headers
    SigP_Gff=gene_pred/combined_sigP_ORF/$Assembler/$Organism/$Strain/"$Strain"_all_secreted_unmerged.gff
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
    ORF_fasta=$(ls gene_pred/ORF_finder/$Assembler/$Organism/$Strain/"$Strain".aa_cat.fa)
    $ProgDir/extract_from_fasta.py --fasta $ORF_fasta --headers $RxLR_Merged_txt > $RxLR_Merged_AA
    printf "Merged RxLR-EER Hmm proteins:\t" >> report.txt
    cat $RxLR_Merged_AA | grep '>' | wc -l >> report.txt
    echo "$Assembler $Strain done"
done
```

```
discovar P.rubi SCRP249
Initial search space (Z):              81145  [actual number of targets]
Domain search space  (domZ):             703  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   50
discovar P.rubi SCRP324
Initial search space (Z):              89423  [actual number of targets]
Domain search space  (domZ):             696  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   217
discovar P.rubi SCRP333
Initial search space (Z):              81434  [actual number of targets]
Domain search space  (domZ):             668  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   208
spades P.rubi SCRP249
Initial search space (Z):              74633  [actual number of targets]
Domain search space  (domZ):             682  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   211
spades P.rubi SCRP324
Initial search space (Z):              83177  [actual number of targets]
Domain search space  (domZ):             681  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   212
spades P.rubi SCRP333
Initial search space (Z):              74477  [actual number of targets]
Domain search space  (domZ):             662  [number of targets reported over threshold]
Merged RxLR-EER Hmm proteins:   206
```

E7) Combining RxLRs from Regex and hmm searches

The total RxLRs are

```bash
for RegexRxLR in $(ls analysis/RxLR_effectors/RxLR_EER_regex_finder/*/*/*/*_ORF_RxLR_EER_regex_merged.txt)
do
    Organism=$(echo $RegexRxLR | rev |  cut -d '/' -f3 | rev)
    Strain=$(echo $RegexRxLR | rev | cut -d '/' -f2 | rev)
    Assembler=$(echo $RegexRxLR | rev | cut -d '/' -f4 | rev)
    Gff=$(ls gene_pred/ORF_finder/$Assembler/$Organism/$Strain/"$Strain"_ORF.gff3)
    Proteome=$(ls gene_pred/ORF_finder/$Assembler/$Organism/$Strain/"$Strain".aa_cat.fa)
    HmmRxLR=$(ls analysis/RxLR_effectors/hmmer_RxLR/$Assembler/$Organism/$Strain/"$Strain"_ORF_RxLR_hmm_merged.txt)
    echo "$Assembler - $Organism - $Strain" >> report.txt
    echo "Number of RxLRs identified by Regex:" >> report.txt
    cat $RegexRxLR | sort | uniq | wc -l >> report.txt
    echo "Number of RxLRs identified by Hmm:" >> report.txt
    cat $HmmRxLR | sort | uniq | wc -l >> report.txt
    echo "Number of RxLRs in combined dataset:" >> report.txt
    cat $RegexRxLR $HmmRxLR | sort | uniq | wc -l >> report.txt
    OutDir=analysis/RxLR_effectors/combined_evidence/$Assembler/$Organism/$Strain
    mkdir -p $OutDir
    cat $RegexRxLR $HmmRxLR | sort | uniq > $OutDir/"$Strain"_total_ORF_RxLR_headers.txt
    ProgDir=/home/adamst/git_repos/tools/seq_tools/feature_annotation
    $ProgDir/gene_list_to_gff.pl $OutDir/"$Strain"_total_ORF_RxLR_headers.txt $Gff ORF_RxLR Name Augustus > $OutDir/"$Strain"_total_ORF_RxLR.gff
    echo "Number of genes in the extracted gff file:" >> report.txt
    cat $OutDir/"$Strain"_total_ORF_RxLR.gff | grep -w 'gene' | wc -l >> report.txt
    echo "" >> report.txt
    echo "$Assembler $Strain done"
done
```

```
discovar - P.rubi - SCRP249
Number of RxLRs identified by Regex:
74
Number of RxLRs identified by Hmm:
50
Number of RxLRs in combined dataset:
84
Number of genes in the extracted gff file:
84

discovar - P.rubi - SCRP324
Number of RxLRs identified by Regex:
276
Number of RxLRs identified by Hmm:
217
Number of RxLRs in combined dataset:
309
Number of genes in the extracted gff file:
309

discovar - P.rubi - SCRP333
Number of RxLRs identified by Regex:
270
Number of RxLRs identified by Hmm:
208
Number of RxLRs in combined dataset:
298
Number of genes in the extracted gff file:
298

spades - P.rubi - SCRP249
Number of RxLRs identified by Regex:
270
Number of RxLRs identified by Hmm:
211
Number of RxLRs in combined dataset:
300
Number of genes in the extracted gff file:
300

spades - P.rubi - SCRP324
Number of RxLRs identified by Regex:
266
Number of RxLRs identified by Hmm:
212
Number of RxLRs in combined dataset:
299
Number of genes in the extracted gff file:
299

spades - P.rubi - SCRP333
Number of RxLRs identified by Regex:
264
Number of RxLRs identified by Hmm:
206
Number of RxLRs in combined dataset:
292
Number of genes in the extracted gff file:
292
```

4.2.c Analysis of RxLR effectors - merger of Augustus / published genes with ORFs

Intersection between the coodinates of putative RxLRs from gene models and ORFs were identified to determine the total number of RxLRs predicted in these genomes.

The RxLR effectors from both Gene models and ORF finding approaches were combined into a single file.

```bash
for MergeDir in $(ls -d analysis/RxLR_effectors/combined_evidence/*/*/*)
do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    Assembler=$(echo "$MergeDir" | rev | cut -f3 -d '/' | rev)
    AugGff=$MergeDir/"$Strain"_total_RxLR.gff
    AugTxt=$MergeDir/"$Strain"_total_RxLR_headers.txt
    AugFa=$(ls gene_pred/codingquary/"$Assembler"/"$Species"/"$Strain"/final/final_genes_combined.pep.fasta)

    ORFGff=$(ls $MergeDir/"$Strain"_total_ORF_RxLR.gff)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Assembler"/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
    ORFsTxt=$(ls $MergeDir/"$Strain"_total_ORF_RxLR_headers.txt)

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

    echo "$Assembler - $Species - $Strain" >> report.txt
    echo "The number of ORF RxLRs overlapping Augustus RxLRs:" >> report.txt
    cat $ORFsInAug | grep -w 'gene' | wc -l >> report.txt
    echo "The number of Augustus RxLRs overlapping ORF RxLRs:" >> report.txt
    cat $AugInORFs | grep -w 'gene' | wc -l >> report.txt
    echo "The number of RxLRs unique to ORF models:" >> report.txt
    cat $ORFsUniq | grep -w 'transcript' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' | wc -l >> report.txt
    echo "The number of RxLRs unique to Augustus models:" >> report.txt
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l >> report.txt
    echo "The total number of putative RxLRs are:" >> report.txt
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalRxLRsTxt
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
    cat $ORFsUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f3 -d ';' | cut -f2 -d '=' >> $TotalRxLRsTxt
    cat $TotalRxLRsTxt | wc -l >> report.txt
    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalRxLRsTxt > $TotalRxLRsGff

    RxLRsFa=$MergeDir/"$Strain"_final_RxLR_EER.fa
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalRxLRsTxt > $RxLRsFa
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalRxLRsTxt >> $RxLRsFa
    echo "The number of sequences extracted is" >> report.txt
    cat $RxLRsFa | grep '>' | wc -l >> report.txt
    echo "$Assembler $Strain done"
done
```

```
discovar - P.rubi - SCRP249
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
discovar - P.rubi - SCRP324
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
discovar - P.rubi - SCRP333
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
spades - P.rubi - SCRP249
The number of ORF RxLRs overlapping Augustus RxLRs:
198
The number of Augustus RxLRs overlapping ORF RxLRs:
197
The number of RxLRs unique to ORF models:
102
The number of RxLRs unique to Augustus models:
34
The total number of putative RxLRs are:
333
The number of sequences extracted is
317
spades - P.rubi - SCRP324
The number of ORF RxLRs overlapping Augustus RxLRs:
212
The number of Augustus RxLRs overlapping ORF RxLRs:
211
The number of RxLRs unique to ORF models:
87
The number of RxLRs unique to Augustus models:
31
The total number of putative RxLRs are:
329
The number of sequences extracted is
291
spades - P.rubi - SCRP333
The number of ORF RxLRs overlapping Augustus RxLRs:
188
The number of Augustus RxLRs overlapping ORF RxLRs:
187
The number of RxLRs unique to ORF models:
104
The number of RxLRs unique to Augustus models:
37
The total number of putative RxLRs are:
328
The number of sequences extracted is
304
```

H) From ORF gene models - Hmm evidence of CRN effectors

A hmm model relating to crinkler domains was used to identify putative crinklers in ORF gene models. This was done with the following commands:

```bash
for Proteome in $(ls gene_pred/ORF_finder/*/*/*/*.aa_cat.fa)
do
    # Setting variables
    Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Assembler=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=analysis/CRN_effectors/hmmer_CRN/$Assembler/$Organism/$Strain
    mkdir -p $OutDir
    # Hmmer variables
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/hmmer
    HmmDir=/home/groups/harrisonlab/project_files/idris/analysis/CRN_effectors/hmmer_models
    # Searches for LFLAK domain
    LFLAK_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_LFLAK.hmm
    HmmResultsLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.txt
    hmmsearch -T 0 $LFLAK_hmm $Proteome > $OutDir/$HmmResultsLFLAK
    echo "Searching for LFLAK domains in: $Organism $Strain" >> report.txt
    cat $OutDir/$HmmResultsLFLAK | grep 'Initial search space' >> report.txt
    cat $OutDir/$HmmResultsLFLAK | grep 'number of targets reported over threshold' >> report.txt
    HmmFastaLFLAK="$Strain"_ORF_CRN_LFLAK_unmerged_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsLFLAK $Proteome > $OutDir/$HmmFastaLFLAK
    # Searches for DWL domain
    DWL_hmm=$HmmDir/Pinf_Pram_Psoj_Pcap_DWL.hmm
    HmmResultsDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.txt
    hmmsearch -T 0 $DWL_hmm $Proteome > $OutDir/$HmmResultsDWL
    echo "Searching for DWL domains in: $Organism $Strain" >> report.txt
    cat $OutDir/$HmmResultsDWL | grep 'Initial search space' >> report.txt
    cat $OutDir/$HmmResultsDWL | grep 'number of targets reported over threshold' >> report.txt
    HmmFastaDWL="$Strain"_ORF_CRN_DWL_unmerged_hmmer.fa
    $ProgDir/hmmer2fasta.pl $OutDir/$HmmResultsDWL $Proteome > $OutDir/$HmmFastaDWL
    # Identify ORFs found by both models
    CommonHeaders=$OutDir/"$Strain"_ORF_CRN_DWL_LFLAK_unmerged_headers.txt
    cat $OutDir/$HmmFastaLFLAK $OutDir/$HmmFastaDWL | grep '>' | cut -f1 | tr -d '>' | sort | uniq -d > $CommonHeaders
    echo "The number of CRNs common to both models are:" >> report.txt
    cat $CommonHeaders | wc -l >> report.txt
    # The sequences will be merged based upon the strength of their DWL domain score
    # For this reason headers as they appear in the DWL fasta file were extracted
    Headers=$OutDir/"$Strain"_CRN_hmmer_unmerged_headers.txt
    cat $OutDir/$HmmFastaDWL | grep '>' | grep -w -f $CommonHeaders | tr -d '>' | sed -r 's/\s+/\t/g'| sed 's/=\t/=/g' | tr -d '-' | sed 's/hmm_score/HMM_score/g' > $Headers
    # As we are dealing with JGI and Broad sequences, some features need formatting:
    ORF_Gff=$(ls gene_pred/ORF_finder/$Assembler/$Organism/$Strain/*_ORF.gff3)
    # Gff features were extracted for each header
    CRN_unmerged_Gff=$OutDir/"$Strain"_CRN_unmerged_hmmer.gff3
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $Headers $ORF_Gff CRN_HMM Name > $CRN_unmerged_Gff
    # Gff features were merged based upon the DWL hmm score
    DbDir=analysis/databases/$Assembler/$Organism/$Strain
    mkdir -p $DbDir
    ProgDir=/home/adamst/git_repos/scripts/phytophthora/pathogen/merge_gff
    $ProgDir/make_gff_database.py --inp $CRN_unmerged_Gff --db $DbDir/CRN_ORF.db
    CRN_Merged_Gff=$OutDir/"$Strain"_CRN_merged_hmmer.gff3
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/merge_sigP_ORFs.py --inp $DbDir/CRN_ORF.db --id LFLAK_DWL_CRN --out $DbDir/CRN_ORF_merged.db --gff > $CRN_Merged_Gff
    # Final results are reported:
    echo "Number of CRN ORFs after merging:" >> report.txt
    cat $CRN_Merged_Gff | grep 'gene' | wc -l >> report.txt
    echo "$Assembler $Strain done"
done
```

```
Searching for LFLAK domains in: P.rubi SCRP249
Initial search space (Z):             716788  [actual number of targets]
Domain search space  (domZ):             311  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP249
Initial search space (Z):             716788  [actual number of targets]
Domain search space  (domZ):             512  [number of targets reported over threshold]
The number of CRNs common to both models are:
233
Number of CRN ORFs after merging:
141
Searching for LFLAK domains in: P.rubi SCRP324
Initial search space (Z):             795934  [actual number of targets]
Domain search space  (domZ):             316  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP324
Initial search space (Z):             795934  [actual number of targets]
Domain search space  (domZ):             517  [number of targets reported over threshold]
The number of CRNs common to both models are:
234
Number of CRN ORFs after merging:
140
Searching for LFLAK domains in: P.rubi SCRP333
Initial search space (Z):             720807  [actual number of targets]
Domain search space  (domZ):             301  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP333
Initial search space (Z):             720807  [actual number of targets]
Domain search space  (domZ):             468  [number of targets reported over threshold]
The number of CRNs common to both models are:
222
Number of CRN ORFs after merging:
132
Searching for LFLAK domains in: P.rubi SCRP249
Initial search space (Z):             645937  [actual number of targets]
Domain search space  (domZ):             301  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP249
Initial search space (Z):             645937  [actual number of targets]
Domain search space  (domZ):             452  [number of targets reported over threshold]
The number of CRNs common to both models are:
226
Number of CRN ORFs after merging:
131
Searching for LFLAK domains in: P.rubi SCRP324
Initial search space (Z):             725179  [actual number of targets]
Domain search space  (domZ):             294  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP324
Initial search space (Z):             725179  [actual number of targets]
Domain search space  (domZ):             440  [number of targets reported over threshold]
The number of CRNs common to both models are:
217
Number of CRN ORFs after merging:
126
Searching for LFLAK domains in: P.rubi SCRP333
Initial search space (Z):             646243  [actual number of targets]
Domain search space  (domZ):             293  [number of targets reported over threshold]
Searching for DWL domains in: P.rubi SCRP333
Initial search space (Z):             646243  [actual number of targets]
Domain search space  (domZ):             438  [number of targets reported over threshold]
The number of CRNs common to both models are:
220
Number of CRN ORFs after merging:
129
```

Extract crinklers from published gene models

```bash
for MergeDir in $(ls -d analysis/CRN_effectors/hmmer_CRN/*/*/*)
do
    Strain=$(echo "$MergeDir" | rev | cut -f1 -d '/' | rev)
    Species=$(echo "$MergeDir" | rev | cut -f2 -d '/' | rev)
    Assembler=$(echo "$MergeDir" | rev | cut -f3 -d '/' | rev)
    AugGff=$(ls $MergeDir/"$Strain"_pub_CRN_LFLAK_DWL.gff)
    AugFa=$(ls gene_pred/codingquary/"$Assembler"/"$Species"/"$Strain"/final/final_genes_combined.pep.fasta)
    ORFsFa=$(ls gene_pred/ORF_finder/"$Assembler"/"$Species"/"$Strain"/"$Strain".aa_cat.fa)
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
    echo "$Assembler - $Species - $Strain" >> report.txt

    echo "The number of ORF CRNs overlapping Augustus CRNs:" >> report.txt
    cat $ORFsInAug | grep -w -e 'transcript' -e 'mRNA' | wc -l >> report.txt
    echo "The number of Augustus CRNs overlapping ORF CRNs:" >> report.txt
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA' | wc -l >> report.txt
    cat $AugInORFs | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $TotalCRNsTxt
    echo "The number of CRNs unique to ORF models:" >> report.txt
    cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' | wc -l >> report.txt
    cat $ORFsUniq | grep -w 'transcript'| grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f4 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt
    echo "The number of CRNs unique to Augustus models:" >> report.txt
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA' | wc -l >> report.txt
    cat $AugUniq | grep -w -e 'transcript' -e 'mRNA'  | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' >> $TotalCRNsTxt

    cat $AugInORFs $AugUniq $ORFsUniq | grep -w -f $TotalCRNsTxt > $TotalCRNsGff

    CRNsFa=$MergeDir/"$Strain"_final_CRN.fa
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $AugFa --headers $TotalCRNsTxt > $CRNsFa
    $ProgDir/extract_from_fasta.py --fasta $ORFsFa --headers $TotalCRNsTxt >> $CRNsFa
    echo "The number of sequences extracted is" >> report.txt
    cat $CRNsFa | grep '>' | wc -l >> report.txt
    echo "$Assembler $Strain done"
done
```

```
discovar - P.rubi - SCRP249
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
discovar - P.rubi - SCRP324
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
discovar - P.rubi - SCRP333
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
spades - P.rubi - SCRP249
The number of ORF CRNs overlapping Augustus CRNs:
123
The number of Augustus CRNs overlapping ORF CRNs:
123
The number of CRNs unique to ORF models:
8
The number of CRNs unique to Augustus models:
3
The number of sequences extracted is
126
spades - P.rubi - SCRP324
The number of ORF CRNs overlapping Augustus CRNs:
120
The number of Augustus CRNs overlapping ORF CRNs:
120
The number of CRNs unique to ORF models:
6
The number of CRNs unique to Augustus models:
4
The number of sequences extracted is
117
spades - P.rubi - SCRP333
The number of ORF CRNs overlapping Augustus CRNs:
120
The number of Augustus CRNs overlapping ORF CRNs:
120
The number of CRNs unique to ORF models:
9
The number of CRNs unique to Augustus models:
1
The number of sequences extracted is
119
```
