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

** N50:
SCRP249: 16613
SCRP324: 19133
SCRP333: 16946 **

** L50:
SCRP249: 1233
SCRP324: 1057
SCRP333: 1210 **

** Number of contigs > 1kb:
SCRP249: 9215
SCRP324: 9110
SCRP333: 8986 **

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
for Strain in SCRP249 SCRP324 SCRP333
do
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    OutDir=$(ls -d assembly/discovar/*/$Strain/assembly/a.final)
    AssFiltered=$OutDir/a.lines.fasta
    AssRenamed=$OutDir/a.lines.renamed.fasta
    echo $AssFiltered
    printf '.\t.\t.\t.\n' > editfile.tab
    $ProgDir/remove_contaminants.py --inp $AssFiltered --out $AssRenamed --coord_file editfile.tab
    rm editfile.tab
done
```

### QUAST used to summarise assembly statistics

Discovar assembly:

```bash
ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Strain in SCRP249 SCRP324 SCRP333
do
    for Assembly in $(ls assembly/discovar/*/$Strain/assembly/a.final/a.lines.renamed.fasta)
    do
        Organism=P.rubi
        OutDir=assembly/discovar/$Organism/$Strain/assembly/a.final/QUAST
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
    for BestAss in $(ls assembly/discovar/*/$Strain/assembly/a.final/a.lines.renamed.fasta)
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
    for BestAss in $(ls assembly/spades/*/$Strain/filtered_contigs/*_500bp_renamed.fasta)
    do
        qsub $ProgDir/rep_modeling.sh $BestAss
        qsub $ProgDir/transposonPSI.sh $BestAss
    done
done   
```

The number of bases masked by transposonPSI and Repeatmasker were summarised using the following commands:

for discovar assembly

```bash
for RepDir in $(ls -d repeat_masked/P.*/*/assembly/a.final_repmask)
do
    Strain=$(echo $RepDir | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $RepDir | rev | cut -f4 -d '/' | rev)  
    RepMaskGff=$(ls $RepDir/assembly_contigs_hardmasked.gff)
    TransPSIGff=$(ls $RepDir/assembly_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
    printf "$Organism\t$Strain\n"
    printf "The number of bases masked by RepeatMasker:\t"
    sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The number of bases masked by TransposonPSI:\t"
    sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
    printf "The total number of masked bases are:\t"
    cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done
```

** SCRP249
The number of bases masked by RepeatMasker:	53632316
The number of bases masked by TransposonPSI:	9072448
The total number of masked bases are:	55262319
SCRP324
The number of bases masked by RepeatMasker:	52598215
The number of bases masked by TransposonPSI:	9123758
The total number of masked bases are:	54313575
SCRP333
The number of bases masked by RepeatMasker:	49285458
The number of bases masked by TransposonPSI:	9172511
The total number of masked bases are:	51115236 **

#Merging RepeatMasker and TransposonPSI outputs

```bash
for File in $(ls repeat_masked/P.rubi/*/assembly/a.final_repmask/*_contigs_softmasked.fa)
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

for SPAdes assembly:

```bash
for RepDir in $(ls -d repeat_masked/spades/P.*/*/filtered_contigs_repmask)
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

** SCRP249
The number of bases masked by RepeatMasker:	23482028
The number of bases masked by TransposonPSI:	5953026
The total number of masked bases are:	25419316
SCRP324
The number of bases masked by RepeatMasker:	23613201
The number of bases masked by TransposonPSI:	5940852
The total number of masked bases are:	25417379
SCRP333
The number of bases masked by RepeatMasker:	23158566
The number of bases masked by TransposonPSI:	5961557
The total number of masked bases are:	25028389 **

#Merging RepeatMasker and TransposonPSI outputs

```bash
for File in $(ls repeat_masked/spades/*/*/filtered_contigs_repmask/*_contigs_softmasked.fa)
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
Gene prediction followed three steps: Pre-gene prediction - Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified. Gene model training - Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline Gene prediction - Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

##Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.

for discovar assemblies:

```bash
ProgDir=/home/adamst/git_repos/tools/gene_prediction/cegma
for Genome in $(ls repeat_masked/P.*/*/assembly/a.final_repmask/*_contigs_unmasked.fa)
do
    echo $Genome
    qsub $ProgDir/sub_cegma.sh $Genome dna
done
```

Outputs were summarised using the commands:

```bash
for Strain in SCRP249 SCRP324 SCRP333
do
for File in $(ls gene_pred/cegma/$Strain/*/*_dna_cegma.completeness_report)
do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Species=$(echo $File | rev | cut -f3 -d '/' | rev)
    printf "$Species\t$Strain\n"
    cat $File | head -n18 | tail -n+4;printf "\n"
done
done >> gene_pred/cegma/discovar/cegma_results_dna_summary.txt

less gene_pred/cegma/discovar/cegma_results_dna_summary.txt
```

** SCRP249
Complete: 93.55%
Partial: 96.37%

SCRP324:
Complete: 94.76%
Partial: 97.58%

SCRP333:
Complete: 94.35%
Partial: 96.77% **

for SPAdes assemblies:

```bash
ProgDir=/home/adamst/git_repos/tools/gene_prediction/cegma
for Genome in $(ls repeat_masked/spades/P.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa)
do
    echo $Genome
    qsub $ProgDir/sub_cegma.sh $Genome dna
done
```

Outputs were summarised using the commands:

```bash
for File in $(ls gene_pred/cegma/P.rubi/*/*_dna_cegma.completeness_report)
do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Species=$(echo $File | rev | cut -f3 -d '/' | rev)
    printf "$Species\t$Strain\n"
    cat $File | head -n18 | tail -n+4;printf "\n"
done >> gene_pred/cegma/spades/cegma_results_dna_summary.txt

less gene_pred/cegma/spades/cegma_results_dna_summary.txt
```

** SCRP249
Complete: 93.95%
Partial: 96.77%

SCRP324:
Complete: 95.56%
Partial: 97.98%

SCRP333:
Complete: 94.35%
Partial: 96.77% **

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
