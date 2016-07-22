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
for Strain in A4 Bc1 Bc23 Nov27 Nov5 Nov71 Nov77 Nov9 ONT3 SCRP245_v2
do
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    OutDir=$(ls -d assembly/spades/*/$Strain/filtered_contigs)
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
for Strain in A4 Bc1 Bc23 Nov27 Nov5 Nov71 Nov77 Nov9 ONT3 SCRP245_v2
do
    for Assembly in $(ls assembly/spades/*/$Strain/filtered_contigs/*_500bp_renamed.fasta)
    do
        Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
        Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
        OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
        qsub $ProgDir/sub_quast.sh $Assembly $OutDir
    done
done
```

** N50:
A4: 18245
BC-16: 437436
BC-23: 18227
NOV-27: 19406
NOV-5: 17887
NOV-77: 18909
ONT-3: 22074
SCRP245_v2: 20105
NOV-71: 20226
NOV-9: 21522
BC-1: 21834 **

** L50:
A4: 1116
BC-16: 59
BC-23: 1119
NOV-27: 1046
NOV-5: 1134
NOV-77: 1102
ONT-3: 917
SCRP245_v2: 994
NOV-71: 1016
NOV-9: 978
BC-1: 954 **

** Number of contigs > 1kb:
A4: 8660
BC-16: 406
BC-23: 8556
NOV-27: 8040
NOV-5: 8760
NOV-77: 8500
ONT-3: 8540
SCRP245_v2: 8584
NOV-71: 7885
NOV-9: 7655
BC-1: 7504 **

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

### Quast

```bash
for Strain in A4 Bc1 Bc23 Nov27 Nov5 Nov71 Nov77 Nov9 ONT3 SCRP245_v2
do
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    OutDir=$(ls -d assembly/spades/*/$Strain/filtered_contigs)
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
for Strain in A4 Bc1 Bc23 Nov27 Nov5 Nov71 Nov77 Nov9 ONT3 SCRP245_v2
do
    for Assembly in $(ls assembly/spades/*/$Strain/filtered_contigs/*_500bp_renamed.fasta)
    do
        Strain=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
        Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
        OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
        qsub $ProgDir/sub_quast.sh $Assembly $OutDir
    done
done
```

** N50:
A4: 18245
BC-16: 437436
BC-23: 18227
NOV-27: 19406
NOV-5: 17887
NOV-77: 18909
ONT-3: 22074
SCRP245_v2: 20105
NOV-71: 20226
NOV-9: 21522
BC-1: 21834 **

** L50:
A4: 1116
BC-16: 59
BC-23: 1119
NOV-27: 1046
NOV-5: 1134
NOV-77: 1102
ONT-3: 917
SCRP245_v2: 994
NOV-71: 1016
NOV-9: 978
BC-1: 954 **

** Number of contigs > 1kb:
A4: 8660
BC-16: 406
BC-23: 8556
NOV-27: 8040
NOV-5: 8760
NOV-77: 8500
ONT-3: 8540
SCRP245_v2: 8584
NOV-71: 7885
NOV-9: 7655
BC-1: 7504 **

#Repeatmasking

Repeat masking was performed and used the following programs: Repeatmasker Repeatmodeler

The best assemblies were used to perform repeatmasking

for BC-16 pacbio data:

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/merged_canu_spades/*/*/95m/filtered_contigs/Bc16_contigs_renamed.fasta)
do
    qsub $ProgDir/rep_modeling.sh $BestAss
    qsub $ProgDir/transposonPSI.sh $BestAss
done
```

for other isolates Illumina data:

```bash
for Strain in A4 Bc1 Bc23 Nov27 Nov5 Nov71 Nov77 Nov9 ONT3 SCRP245_v2
do
    for BestAss in $(ls assembly/spades/*/$Strain/filtered_contigs/*_500bp_renamed.fasta)
    do
        qsub $ProgDir/rep_modeling.sh $BestAss
        qsub $ProgDir/transposonPSI.sh $BestAss
    done
done   
```

The number of bases masked by transposonPSI and Repeatmasker were summarised using the following commands:

```bash
for RepDir in $(ls -d repeat_masked/P.*/*/filtered_contigs_repmask)
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

** A4
The number of bases masked by RepeatMasker:	24836372
The number of bases masked by TransposonPSI:	6237528
The total number of masked bases are:	26598776
Bc1
The number of bases masked by RepeatMasker:	24254593
The number of bases masked by TransposonPSI:	6219671
The total number of masked bases are:	26154357
Bc23
The number of bases masked by RepeatMasker:	23771588
The number of bases masked by TransposonPSI:	6101880
The total number of masked bases are:	25516134
Nov27
The number of bases masked by RepeatMasker:	24653573
The number of bases masked by TransposonPSI:	6209723
The total number of masked bases are:	26343538
Nov5
The number of bases masked by RepeatMasker:	24011096
The number of bases masked by TransposonPSI:	6242538
The total number of masked bases are:	25856769
Nov71
The number of bases masked by RepeatMasker:	24200190
The number of bases masked by TransposonPSI:	6080704
The total number of masked bases are:	25824977
Nov77
The number of bases masked by RepeatMasker:	24253868
The number of bases masked by TransposonPSI:	6250930
The total number of masked bases are:	26117699
Nov9
The number of bases masked by RepeatMasker:	24774161
The number of bases masked by TransposonPSI:	6290033
The total number of masked bases are:	26664169
ONT3
The number of bases masked by RepeatMasker:	25224812
The number of bases masked by TransposonPSI:	6238377
The total number of masked bases are:	26981713
SCRP245_v2
The number of bases masked by RepeatMasker:	23381847
The number of bases masked by TransposonPSI:	6037837
The total number of masked bases are:	25248164 **

#Gene Prediction
Gene prediction followed three steps: Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified. Gene model training - Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline Gene prediction - Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

##Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

```bash
ProgDir=/home/adamst/git_repos/tools/gene_prediction/cegma
for Genome in $(ls repeat_masked/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa)
do
    echo $Genome
    qsub $ProgDir/sub_cegma.sh $Genome dna
done
```

Outputs were summarised using the commands:

```bash
for File in $(ls gene_pred/cegma/*/*/*_dna_cegma.completeness_report)
do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Species=$(echo $File | rev | cut -f3 -d '/' | rev)
    printf "$Species\t$Strain\n"
    cat $File | head -n18 | tail -n+4;printf "\n"
done >> gene_pred/cegma/cegma_results_dna_summary.txt

less gene_pred/cegma/cegma_results_dna_summary.txt
```

** A4
Complete: 95.16%
Partial: 97.98%

Bc16
Complete: 94.35%
Partial: 96.37%

Bc1
Complete: 95.16%
Partial: 97.58%

Bc23
Complete: 95.16%
Partial: 97.58%

Nov27
Complete: 94.76%
Partial: 97.18%

Nov5
Complete: 94.76%
Partial: 97.18%

Nov71
Complete: 95.16%
Partial: 97.98%

Nov77
Complete: 94.76%
Partial: 97.18%

Nov9
Complete: 94.35%
Partial: 97.18%

ONT3
Complete: 95.16%
Partial: 97.18%

SCRP245_v2
Complete: 95.16%
Partial: 97.18% **

#Gene prediction
Gene prediction was performed for the P. fragariae genomes. Two gene prediction approaches were used:

Gene prediction using Braker1 and Prediction of all putative ORFs in the genome using the ORF finder (atg.pl) approach.

##Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to P. fragariae genomes.

qc of RNA seq data was performed as part of sequencing the 10300 genome:

```bash
FileF=qc_rna/raw_rna/genbank/P.cactorum/F/SRR1206032_trim.fq.gz
FileR=qc_rna/raw_rna/genbank/P.cactorum/R/SRR1206033_trim.fq.gz
```

#Aligning

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa)
do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNA in $(ls qc_rna/genbank/P.cactorum/10300/*/*_trim.fq.gz)
    do
        Timepoint=$(echo $RNA | rev | cut -f1 -d '/' | rev | sed 's/_trim.*//g')
        echo "$Timepoint"
        OutDir=alignment/$Organism/$Strain/$Timepoint
        ProgDir=/home/adamst/git_repos/tools/seq_tools/RNAseq
        qsub $ProgDir/tophat_alignment_unpaired.sh $Assembly $RNA $OutDir
    done
done
```

#Braker prediction

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa)
do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    mkdir -p alignment/$Organism/$Strain/concatenated
    samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
    alignment/$Organism/$Strain/SRR1206032/accepted_hits.bam \
    alignment/$Organism/$Strain/SRR1206033/accepted_hits.bam
    OutDir=gene_pred/braker/$Organism/"$Strain"_braker
    AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
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
for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa)
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
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa)
do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain
    CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
    ProgDir=/home/adamst/git_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```

Then, additional transcripts were added to Braker1 gene models, when CodingQuarry genes were predicted in regions of the genome, not containing Braker1 gene models:


    for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_pacbio/*/augustus.gff3 | grep -e 'Fus2'); do
        Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g')
        Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
        echo "$Organism - $Strain"
        # BrakerGff=gene_pred/braker/$Organism/$Strain/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
        Assembly=$(ls repeat_masked/$Organism/$Strain/*/"$Strain"_contigs_softmasked.fa)
        CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
        PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
        AddDir=gene_pred/codingquary/$Organism/$Strain/additional
        FinalDir=gene_pred/codingquary/$Organism/$Strain/final
        AddGenesList=$AddDir/additional_genes.txt
        AddGenesGff=$AddDir/additional_genes.gff
        FinalGff=$AddDir/combined_genes.gff
        mkdir -p $AddDir
        mkdir -p $FinalDir

        bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
        bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
        $ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
        $ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
        ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
        # GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/additional/additional_genes.gff
        # GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/out/PredictedPass.gff3

        $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
        $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
        cp $BrakerGff $FinalDir/final_genes_Braker.gff3
        $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
        cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
        cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
        cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
        cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

        GffBraker=$FinalDir/final_genes_CodingQuary.gff3
        GffQuary=$FinalDir/final_genes_Braker.gff3
        GffAppended=$FinalDir/final_genes_appended.gff3
        cat $GffBraker $GffQuary > $GffAppended

        # cat $BrakerGff $AddDir/additional_gene_parsed.gff3 | bedtools sort > $FinalGff
    done
The final number of genes per isolate was observed using:

for DirPath in $(ls -d gene_pred/codingquary/F.*/*/final | grep -w -e'Fus2'); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
