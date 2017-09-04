#After being run through the Genbank contamination screen, seven of thirteen genomes have had contaminants identified from various sources.

The identity of the contaminants has been passed to Andy to allow tweaking of Deconseq for screening assemblies pre-submission.

##One of the assemblies had a number of hits to Bacillus spp., run my isolates through to see if there are any hits in the other isolates.

###Contigs were identified that had BLAST hits to non-phytophthora genomes

Illumina rubi genomes

```bash
for Assembly in $(ls ../phytophthora_rubi/assembly/spades/*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta)
do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    # Exclude_db="bact,virus,hsref"
    Exclude_db="bacillus"
    Good_db="phytoph"
    AssemblyDir=$(dirname $Assembly)
    # OutDir=$AssemblyDir/../deconseq
    OutDir=$AssemblyDir/../deconseq_Paen
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    qsub $ProgDir/sub_deconseq.sh $Assembly $Exclude_db $Good_db $OutDir
done
```

Illumina rubi genomes

```bash
# for File in $(ls assembly/spades/P.*/*/deconseq/log.txt); do
for File in $(ls ../phytophthora_rubi/assembly/spades/P.*/*/deconseq_Paen/log.txt)
do
    Name=$(echo $File | rev | cut -f3 -d '/' | rev)
    Good=$(cat $File |cut -f2 | head -n1 | tail -n1)
    Both=$(cat $File |cut -f2 | head -n2 | tail -n1)
    Bad=$(cat $File |cut -f2 | head -n3 | tail -n1)
    printf "$Name\t$Good\t$Both\t$Bad\n"
done
```

```
SCRP249	14020	3	0
SCRP324	13945	1	1
SCRP333	13678	1	0
```

Contaminant organisms identified by NCBI BLAST

```
Only contig identified as a contaminent here BLASTs to a Phytophthora genome, do not believe this. The false positive rate is too high.
```

Phytophthora rubi

```bash
for Assembly in $(ls assembly/spades/P.*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e 'SCRP324')
do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
    mkdir -p $NCBI_report_dir
done
```

##These files were used to correct assemblies

```bash
for Assembly in $(ls assembly/spades/P.*/*/deconseq_Paen/contigs_min_500bp_filtered_renamed.fasta | grep -e 'SCRP324')
do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/Contamination*.txt)
    if [[ $NCBI_report ]]
    then
        echo "Contamination report found"
    else
        NCBI_report=genome_submission/$Organism/$Strain/initial_submission/no_edits.txt
        printf "Exclude:\nSequence name, length, apparent source\n" > $NCBI_report
    fi
    OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
    mkdir -p $OutDir
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
    # $ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done
```

###Summarise assemblies using QUAST

```bash
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta)
do
    Kmer=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    # OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    OutDir=$(dirname $Assembly)
    ProgDir=/home/adamst/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

###Assembly statistics were summarised using

```bash
for File in $(ls assembly/spades/P.*/*/ncbi_edits/filtered_contigs/report.tsv)
do
    echo "$File" | rev | cut -f3 -d '/' | rev
    cat $File | cut -f2
done
```
