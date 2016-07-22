#Discovar assembly

*Phytophthora rubi* genomes were also assembled using Discovar. These commands were run via an interactive job within a screen session as opposed to via the grid engine to more efficiently use cluster resources.

All commands run in /home/groups/harrisonlab/project_files/phytophthora_rubi

```bash
qlogin -pe smp 24 -l virtual_free=15.7G -l h=blacklace01.blacklace
```

# Perform genome assembly of diploid organisms using discovar

Usage="sub_discovar.sh <F_read.fa> <R_read.fa> <output_directory>"

#Collect inputs

```bash
R1=qc_dna/paired/P.rubi/SCRP249/F/SCRP249_S1_L001_R1_001_trim.fq.gz
R2=qc_dna/paired/P.rubi/SCRP249/R/SCRP249_S1_L001_R2_001_trim.fq.gz
OutDir=assembly/discovar/P.rubi/SCRP249
CurPath=$PWD
WorkDir=tmp/discovar

echo "Running Discovar with the following reads:"
echo "$R1"
echo "$R2"
echo "Outputting to the following location:"
echo "$CurPath/$OutDir "
```

#Copy files

```bash
F_Read=$(basename $R1)
R_Read=$(basename $R2)

mkdir -p $WorkDir
cd $WorkDir

cp $CurPath/$R1 $WorkDir/$F_Read
cp $CurPath/$R2 $WorkDir/$R_Read
```

#Run Discovar

```bash
DiscovarDeNovo \
  READS=$F_Read,$R_Read \
  OUT_DIR=$WorkDir/assembly \
  MEMORY_CHECK=True
```

#Cleanup

```bash
mkdir -p $CurPath/$OutDir
rm $WorkDir/$F_Read
rm $WorkDir/$R_Read
cp -r $WorkDir/assembly $CurPath/$OutDir/.
```
