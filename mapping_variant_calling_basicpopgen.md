# GBS data processing for *Ericameria nauseosa nauseosa* from Pyramid lake.

This document presents our workflow and rationale for genotype inference from high throughput sequencing of reduced representation libraries (GBS, RADseq, ddRADseq, etc). Several canned software packages or computational workflows exist for handling this type of data. These methods rely on set thresholds for sequencing coverage depth per locus to call hard genotypes. The biggest cost of using these methods is throwing away much, if not most, of the data. Examples of such software/workflows include:

* [Stacks](http://catchenlab.life.illinois.edu/stacks/)
* [Ddocent](http://www.ddocent.com//)
* [ipyrad](https://ipyrad.readthedocs.io/en/master/)

See [Nielsen et al. 2011](./papers/Nielsen_etal_2011.pdf) and [Buerkle and Gompert 2013](./papers/Buerkle_Gompert_2013.pdf) for articulate thoughts about this.

## Table of Contents

 1. Species Information
 2. Directory Structure
 3. DNA Extraction
 4. Library Preparation
 5. Sequencing
 6. Data Cleaning
    1. Cleaning Contaminants
    2. Barcode Parsing
    3. Splitting fastqs
 7. Denovo Assembly
    1. Prepare Directories and Files
    2. Generate Unique Sequence Sets
    3. Assemble from Sequence Sets
 8. Mapping Reads with BWA
    1. Prepare Directories and Files
    2. Map, Sort, and Index
    3. Explanation of Script
 9. Build Mpileup with bcftools
    1. Prepare Directories
    2. Pileup, Call, and Filter
    3. Understanding bcftools Parameters
10. Convert bcf to vcf
    1. Identify SNPs
    2. Understanfing vcftools Parameters
11. Filtering
    1. Filtering on Individuals
    2. Filtering on Loci
12. Uncategorized
    1. Reference-based assembly (*T. podura*)
    2. Reheader vcf files?
13. Appendices
    1. Appendix 1: Useful Commands

## Focal Species: *Ericameria nauseosa nauseosa*

# Background and rationale here

Rubber rabbitbrush *Ericameria nauseosa* is a perennial shrub in the Aster family (Asteraceae), often co-occuring with another member, sagebrush (*Artemisia tridentata*). *E. nauseosa* is broadly distributed across desert, woodland, and arid montane habitats of the west from northern Mexico to southern Canada. It is a foundational plant species that thrives as a primary successional colonizer of dryland and distrubed sites and acts as a soil stabilizer. *E. nauseosa* is an important resource for birds, mammals, and insects. It hosts large diversitys and quantities of gall-forming insects and serves as a critical late-season pollen source for diverse pollinator communities.  

*E. nauseosa* displays a wide variety of phenotypic variation and has been classified as two recogonized subspecies, *E. n. nauseosa* (grey-stemmed) and *E. n. consimillis* (green-stemmed), and 22 varities; however, underlying genetic data do not support most of these variety taxanomic distictions. Regardless, phenotypic and genetic divergence has been identified across environmental gradients for the two subspecies. Phenotypic differences in morphology and chemistry are predicted to drive distinct insect and herbivore communities.

Insert information about sampling 

## Directory Structure - optional

* Create a project folder in your working directory called KRLA with `mkdir ERNA`
* The following directory structure should be made throughout the tutorial

```mermaid
flowchart TD;
    A(personal directory <br> /working/ebrewer/) --> B(species folder <br> /ebrewer/ERNA/)
    B --> C(assembly)
    B --> N (select_seqs)
    B --> D(bwa)
    B --> E(fastq)
    B --> F(scripts)
    E --> G(fastq files <br> e.g. *.fastq.gz)
    C --> H(de novo assembly <br> e.g. rf.*)
    C --> I(indexed assembly files <br> e.g. *.amb, *.ann, *.bwt, *.pac, *.sa)
    C --> J(alt_assemblies <br> OPTIONAL)
    J --> K(seq subset files <br> e.g. k4.i2.seqs)
    J --> L(assembly files <br> e.g. rf.4.2.92)
    D --> M(mapped/sorted reads <br> + index files <br> e.g. *.bam and *.bam.bai)
```

## DNA extraction

* Performed by AG Biotech in 12/2022
* DNA information with plate maps and IDs can be found at [insert doc path here]

## Library preparation

* Performed in Parchman Lab in 12/2022
* [add upstream path]KRLA\_RFseq\_mastermixcockatils.xlsx contains information about reagents used in library prep
* R/L and PCR for plates 1-6

## Sequencing

* 1 lane of S2 chemistry NovaSeq data at UTGSAF in 03/2023
* sequence results in /archive/parchman_lab/rawdata_to_backup/
## Data cleaning

### Cleaning contaminants

* **Goal:** Remove reads that are not from the target organims
* Performed 06/09/23

1. Create a fastq folder inside KRLA folder:

   ```sh
   mkdir fastq
   ```

2. Copy the sequencing results to your fastq folder:

   ```sh
   cd fastq
   cp /archive/parchman_lab/rawdata_to_backup/KRLA/KRLA_S1_L001_R1_001.fastq.gz .
   ```

3. Decompress the file with:

   ```sh
   gunzip KRLA_S1_L001_R1_001.fastq.gz
   ```

4. Count the reads **before** cleaning with:

   ```sh
   nohup grep -c "^@" KRLA_S1_L001_R1_001.fastq > KRLA_number_of_rawreads.txt &
   ```

   * `nohup` (no hang up) allows the command to keep running, even if you close the terminal, lose connection with the server, or log out.
   * `&` allows the command to run in the background.
   * `grep -c` will count the number of occurances of the pattern ("^@") in the file (KRLA_S1_L001_R1_001.fastq).
   * the result file (KRLA_number_of_rawreads.txt) will contain the initial number of reads.
   * **Reads before cleaning:** n

5. Run the cleaning script (note: the script will need to be edited to change the paths to the appropriate files):

   ```sh
   nano ../scripts/cleaning_bash_KRLA.sh # this is where you edit the path names
   conda deactivate # to make sure you are using the loaded modules
   module load fqutils/0.4.1
   module load bowtie2/2.2.5
   bash cleaning_bash_KRLA.sh &
   ```

   * Run this in the fastq directory
   * `module load` is needed to allow access to the modules used in the script
   * fqutils and bowtie2 are used in the cleaning script
   * cleaning\_bash\_KRLA.sh uses the tapioca pipeline's contamination analysis script to remove artifacts from:
      * Illumina oligos
      * PhiX
      * *E. coli*

6. Check that KRLA.clean.fastq has been created (from cleaning script), then remove the raw data file (from KRLA/fastq - a backup is stored at /archive/parchman\_lab/rawdata\_to\_backup/KRLA) :

   ```sh
   ls # look for KRLA.clean.fastq
   rm -rf KRLA_S1_L001_R1_001.fastq &
   ```

7. Count the reads **after** cleaning with:

   ```sh
   nohup cd-hit-est -i k4.i2.seqs -o rf.4.2.92 -M 0 -T 0 -c 0.92 &>/dev/null &
   ```

   ```sh
   nohup grep -c "^@" KRLA.clean.fastq > KRLA_clean_reads.txt &
   ```

   * reads after cleaning are typically ~75% of reads before cleaning
   * **Reads after cleaning:** n

### Barcode parsing

* **Goal:** Remove barcode and adapter sequences from reads, place this information (what individual the reads came from) in the read ID line

1. Run parsing script:

   ```sh
   cp /working/parchman/KRLA/parse_barcodes769.pl ../scripts
   cp /working/parchman/KRLA/KRLA_barcode_key.csv ../scripts
   nohup perl ../scripts/parse_barcodes769.pl ../scripts/KRLA_barcode_key.csv KRLA.clean.fastq A00 &>/dev/null &
   ```

   * Run this in the fastq directory
   * KRLA\_barcode\_key.csv provides the barcodes for each individual
   * KRLA.clean.fastq are the clean reads from contaminant cleaning
   * A00 is the sequencer ID (first 3 characters after @ in the fastq identifier)
   * &>/dev/null prevents display of stdout

2. View the parsing results (created by the parsing script):

   ```sh
   less parsereport_KRLA.clean.fastq
   ```

   * Results compare good and bad mids (molecular IDs), bad mids are usually <5% of total mids
   * The other removed results are typically insignificant
   * **Parse report:**
      * Good mids count: 1571963061
      * Bad mids count: 73508410
      * Number of seqs with potential MSE adapter in seq: 321195
      * Seqs that were too short after removing MSE and beyond: 428

### Splitting fastqs

* **Goal:** Sort the reads that are in the one .fastq file into multiple .fastq files, one for each individual

1. Create an IDs file:

   ```sh
   cut -f 3 -d "," KRLA_barcode_key.csv | grep "_" > KRLA_ids_noheader.txt
   ```

   * `-f 3` option looks at the 3rd field (column)
   * `-d ","` option specifies comma as delimiter
   * `cut` extracts that column and pipes it to `grep`
   * `grep` ensures the ID column extracted contains an underscore (this removes the header line / label)

2. Split the (cleaned and parsed) fastq files by individual ID:

   ```sh
   cp /working/parchman/KRLA/splitfastqs/splitFastq_universal_regex.pl ../scripts
   nohup perl ../scripts/splitFastq_universal_regex.pl KRLA_ids_noheader.txt parsed_KRLA.clean.fastq &>/dev/null &
   ```

   * Run this in the fastq directory
   * The perl script creates a separate fastq file for each individual's reads

3. Compress all of the resultant fastq files:

   ```sh
   nohup gzip *fastq &>/dev/null &
   ```
# E. Brewer starts here

# EB/TLP: Start ERNA reference assembly here

Emily: maybe write out some basic features or questions for discussion about the ERNA genome haplotype.

## Mapping reads to reference genome with `bwa`

### *Ericameria nauseosa nauseosa* reference genome being used for read mapping. Must be indexed for efficient mapping with `bwa`

- Mapping is being conducted in `/working/ebrewer/bwa_mapping`
- `.fastq` files for each of the 400 ERNA individuals are in `/working/ebrewer/fastq`

2. Indexing reference genome with `bwa`
- Builds a 'lookup system' that makes read mapping fast without scanning the entire genome each time
- Converts the raw reference fasta sequence into a compresed, searchable data structure using BWT and FM-index. Allows bwa to quickly and efficiently find where sequence reads match with the genome during alignment.

   ```sh
   module load bwa/0.7.17-r1188
   bwa index -p ernacon -a bwtsw ojincantatabio-uni4263-hap2-mb-hirise-zpxz8__01-06-2024__final_assembly.fasta &
   ```
   

### Map, sort and index

1. Load both `bwa` and `samtools` which are required for running the mapping script

   ```sh
   module load bwa/0.7.17-r1188
   module load samtools/1.10
   ```

2. Run the script `bwa_mem2sorted_bam.sh` using the following nohup settings. Here we are calling the bash script from within the bwa directory. The shell script needs editing (try `nano` for this) to specify the name of the reference used above in the indexing step. FOr this project that name is `ernacon`. Note the relative paths in call below:

   ```sh
   nohup bash ../scripts/bwa_mem2sorted_bam.sh 2> /dev/null &
   ```

   * Running the script in this way prevents the process from being interrupted (i.e. you can disconnect from the server while this runs) while also capturing progress print statements in `nohup.out`. You can re-login to the server and check the progress of mapping by going into the `bwa` directory and entering the following:

# work done to here 10/28 TLP

      ```sh
      tail -n 1 nohup.out
      ```

   * This step took **~6 hours** using **24 nodes** on ponderosa for **93 individuals** in the POMA dataset.

### Explanation of `bwa_mem2sorted_bam.sh`

The contents of the previous script is the following:

```sh
#!/bin/bash

ctr=1
fTotal=$(ls ../fastq/*.gz -1 | wc -l)

for file in ../fastq/*.gz
   do
   if test -f "$file"
   then
      fPrefix=$(echo "$file" | sed 's|.*/||' | cut -f1 -d ".")
      echo mapping and sorting individual "$ctr" of "$fTotal"
      bwa mem -t 24 POMA_ref "$file" | \
      samtools view -u -b - | \
      samtools sort -l0 -@24 -o "$fPrefix.sorted.bam"
      ((ctr++))
   fi
done
for sBam in *.sorted.bam
   do
   if test -f "$sBam"
   then
      samtools index -@24 "$sBam"
   fi
done
```

We should explain the steps that are happening here particularly any settings used with

* `bwa mem` - maps sequences to the reference, creating a .sam file
  * `-t` - number of threads used
* `samtools view` - converts .sam format to .bam format to save space
  * `-u` - output uncompressed data
  * `-b` - output in .bam format
* `samtools sort` - sorts .bam files by position on the reference
  * `-l` - output compression level (l0) = uncompressed
  * `-@` - number of threads used (@24 = 24 threads)
  * `-o` - output file name
* `samtools index` - indexes the .bam file for faster search, creating .bai files
  * `-@` - number of threads used (@24 = 24 threads)


## Build pileup and variant call with BCFTools

### Prepare Directoties

1. Make vcf directory

   ```sh
   cd ..
   mkdir vcf
   ```

2. Make list of bam files from, do this in the `bwa/` directory

```sh
ls *.sorted.bam | sed 's/\*//' > bam_list.txt
```

### Pileup, call, and filter

1. Lets do this in the `bwa/` directory, should now have all of the .bam and .bam.bai files as well as a copy of the reference genome (`ojincantatabio-uni4263-hap2-mb-hirise-zpxz8__01-06-2024__final_assembly.fasta`, and associated index files `ernacon*`). We will use `bcftools 1.9` and run the following:
   ```sh
   module load bcftools/1.9
   ```

   ```sh
   nohup bcftools mpileup -a DP,AD,INFO/AD -C 50 -d 250 -f ojincantatabio-uni4263-hap2-mb-hirise-zpxz8__01-06-2024__final_assembly.fasta -q 30 -Q 20 -I -b bam_list.txt -o ERNAdenovo.bcf 2> /dev/null &
   ```
### step-by-step
* `bcftools mpileup` -- command reads one or more BAM files and creates a BCF file describing the pileup (base information per genomic position)
   * summarizing evidence for SNPs/indels 
* `-a DP, AD, INFO/AD` -- includes extra annotations in the output: 
   * `DP`: read depth per position
   * `AD`: allele depth 
   * `INFO/AD`: adds allele depth on the INFO field 
* `-C 50` -- adjusts mapping quality for reads with excessive mismatches 
   * Helps reduce false positive in variant calling 
* `-d 250` -- maximum per-BAM depth
   * pileup will ignore sites with coverage deeper than 250 to prevent bias from high coverage regions
* `-f ojincantatabio-uni4263-hap2-mb-hirise-zpxz8__01-06-2024__final_assembly.fasta` -- reference genome file 
* `-q 30` -- minimum mapping quality to include a read (filters out poorly mapped reads)
* `-Q 20` -- minimum base quality for a base to be considered (filters out unreliable base calls)
* `-I` -- ignore indels, only report SNP information 
* `-b bam_list.txt` -- provides a text file listing all BAM files to process
* `-o ERNAdenovo.bcf` -- output BCF file

   ```sh
   nohup bcftools call -v -m -f GQ ERNAdenovo.bcf -O z -o ERNAdenovo.vcf.gz 2> /dev/null &
   ```
### step-by-step
* `bcftools call` -- calls variants (SNPs and indels) from BCF file
   * looks at pileup of reads at each genomic position and determines whether there is evidence for a variant relative to reference
* `ERNAdenovo.bcf` -- BCF file previously generated
* `-v` -- output only variant sites 
   * positions that match reference exactly are excluded from the output
* `-m` -- multiallelic calling model, can handle positions with more than one alternative allele
* `f GQ` -- specifies format tags to calculate for genotypes
   * `GQ` = Genotype Quality, a measure of confidence in called genotype
* `-0 z` -- output file format, compressed 
* `-o ERNAdenovo.vcf.gz` specifies name of output file



### Understanding bcftools parameters

* -C --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]
* -d --max-depth INT     max per-file depth; avoids excessive memory usage [250]
* -f --fasta-ref FILE    faidx indexed reference sequence file
* -q --min-MQ INT        skip alignments with mapQ smaller than INT [0]
* -Q --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
* -I --skip-indels       do not perform indel calling
* -b --bam-list FILE     list of input BAM filenames, one per line
* -O --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
* -o --output FILE       write output to FILE [standard output]

# TLP stopped here, talk to Seth.
## Preview vcf to be sure that things are looking good in terms of expected number of variants in the raw .vcf. 

Once the unfiltered vcf looks good, move to variant_filtering

   ```sh
   module load vcftools/0.1.14
   ```

Using vcftools to do the bare minimum filtering, which entails removing any indels, removing loci with more than 2 alleles, thinning by 100 (because our sampled genomic regions are 100 bases in length)

   ```sh
   vcftools \
   --remove-indels \
   --min-alleles 2 \
   --max-alleles 2 \
   --thin 100 \
   --remove-filtered-all \
   --recode \
   --recode-INFO-all \
   --gzvcf POMA.vcf.gz \
   --out 
   ```

```sh
vcftools --gzvcf POMA.vcf.gz --out POMA.vcf.gz --missing-indv

```

Below we are making a list of individuals that have too much missing data to move forward with. Here we taking individuals that have data at 50% or more of loci.

```sh
   awk '$5 > 0.5 {print $1}' POMA.imiss | tail -n +2 > indmiss50.txt
```

Running vcftools to make a vcf that contains only the individuals specified in indmiss50.txt made above. First step is making the list of individuals with TOO MUCH missing data.

```sh
   awk '$5 > 0.5 {print $1}' POMA.imiss | tail -n +2 > indmiss50.txt
```

Filter the full vcf using the `--exclude` argument and the indmiss50.txt file made above.

```sh
vcftools --gzvcf POMA.vcf.gz --exclude indmiss50.txt --maf 0.04 --max-meanDP 100 --min-meanDP 2 --minQ 20 --missing .7 --recode --recode-INFO-all --remove-filter-all --out POMA.04.maxdp100.mindp2.miss70.vcf
```


### Understanding vcftools Parameters

* -v --variants-only             output variant sites only
* -c --consensus-caller          the original calling method (conflicts with -m)
* -f --format-fields <list>      output format fields: GQ,GP (lowercase allowed) []
* -p --pval-threshold <float>    variant if P(ref|D)<FLOAT with -c [0.5]
* -P --prior <float>         mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
* -O --output-type <b|u|z|v>     output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]

## Filtering

### Filtering on Individuals

* **Coverage:** Also can be thought of as depth. See 3mapping. Calculated on bam files. Average read count per locus per individual.
* **Missing:** Proportion of missing data allowed across all loci for individual. Common and high in GBS/RADseq data. Kinda an issue all around. Many methods, including PCA (all ordination methods), require a complete matrix with no missing data. Additionally, PCA will cluster by missing data with individuals with higher missing data clustering closer to the center and get this "fan" effect. Can be the same for coverage too. This (among other reasons) is why people use a variance-covariance matrix of genetic data to do ordinations. Other methods involve imputation. This can be fancy and use phased haplotype data OR simply, when you z-score, (g - mean(g))/sd(g), your genotype data across each locus, you make all missing data equal to 0 or Mean (i.e., the global allele frequency). There's more to this standardization, see [Patterson et al. 2006](https://dx.plos.org/10.1371/journal.pgen.0020190) for more info. See PCAsim_ex in examples directory for showing all these issues. This is another reason to use entropy. Entropy is a hierarchical bayesian model so it gets an updated genotype estimate for each missing value based on genotype likelihoods across loci, individuals, and the allele frequency of the cluster/deme that individual assigns to.

### Filtering on Loci

* **Biallelic:** Only keep biallelic SNPs. Multiallelic SNPs are rare at the time scale we work (Citation??) and also, mathematical nightmare and we have enough data so just ignore. Everyone does unless deep time phylogenetics.
* **thin:** Keeps one locus within a specified range. Not 100% how it decides with one to keep. I think it's on quality or depth. This is a necessary step as loci in close physical are prone to sequencing error and linkage disequalibrium (LD) confounds many different population genetic parameters. For de novo reference assemblies, we thin to 100 as contigs/reads are ~92 bp in length. This keeps one locus per contig to control for LD and sequencing error, which is really common in pop gen and necessary for many analyses.
* **max-missing** = max proportion of missing data per locus
* **MAF** = minor allele frequency. Proportion of individuals a alternate allele needs to be present in order for that locus to be kept as a SNP. (e.g. maf = 0.02 for 250 individuals means that an alternate allele needs to be present in at least 5 individuals to be kept) Many papers have shown this is a big issue in clustering and demography (Citation). We do this a second time near the end if we removed individuals during missing data filtering.
* **Mean Depth:** Average allelic depth or read depth per locus. Too low could be sequencing error, too high could be PCR replication artifact (Citation).
* **Qual:** Locus Quality. You can look up the math. Usually above 20-30 is good but given our coverage and number of individuals, we can usually go way higher.
* **Fis:** Inbreeding coefficient. This is a contentous topic. This has to do with paralogs or paralogous loci. This is where loci map to multiple regions of the genome. Issues in highly repeative genomes. Usually leads to an excess of heterozygotes. Filtering on negative Fis can help. See these two McKinney papers [1](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12763), and [2](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12613). Katie and others in the lab use his package called HDPlot to deal with this.

## Uncategorized

### Reference-based assembly (*T. podura*)

1. Refernce file at: /working/parchman/tpodura/raw_ind_fastqs/Tpodura_consensus.fa
2. Make index for `bwa`

   ```sh
   module load bwa/0.7.17-r1188
   bwa index -p Tpodura_consensus -a bwtsw Tpodura_consensus.fa &
   ```

    * `bwa` wrapper, runbwa_memTLP.pl, modified to run the `mem` algorithim (rather than aln and samse), and used bwa 0.7.17-r1188. Parameter settings are described within the wrapper, more info with `bwa mem`

      ```sh
      module load bwa/0.7.17-r1188
      perl runbwa_memTLP.pl  *fastq &
      ```

### Reheader vcf files?

1. Make id file for reheadering

   ```sh
   ls *fastq > fastqs.txt
   sed -s "s/.fastq//" fastqs.txt > Tpod_ids_col.txt
   ```

2. Reheader vcf

   ```sh
   module load bcftools/1.9
   module load vcftools/0.1.14
   bcftools reheader -s Tpod_ids_col.txt tpod.vcf -o rehead_tpod.vcf
   ```

3. Initial filtering

   ```sh
   vcftools --vcf rehead_tpod.vcf --out variants_maf5_miss5 --remove-filtered-all --maf 0.03 --max-missing 0.5 --recode --thin 100
   ```

## Appendices

### Appendix 1: Reference of useful commands (as far as what Parchman lab people like)

* `parallel` - easiest way to run jobs in parallel across CPUs. See [GNU Parallel tutorial](https://www.gnu.org/software/parallel/parallel_tutorial.html)
* `-j` - max jobs to run in parallel
* `--no-notice` - eliminates excessive notices printing to shell when running jobs
* `:::` - followed by list of input files to process
* `nohup <command> &> /dev/null &` - a way to run a longer background process that won't be interupted
  * `nohup` - keeps process running even if shell or terminal is exited (i.e. long jobs don't get terminated)
  * `&` - process is put in background (have access to shell while process is running)
  * `/dev/null` - essentially a black hole to direct st. output from nohup into assuming you've already captured the output of interest
  * can do similar things with `screen` but `nohup` is simpler and enough for most of the use cases here
* `time <command>` - prints the time a process takes after it completes
  * will generate 3 different times, but "real" is what you're usually interested in
  * useful for testing pararmeters of parallelization and getting idea of how long different tasks in pipeline take
* `du -sch <files/dir/etc.> | tail -n 1` - way to see how much disk space a set of files is using, useful if a lot of temporary/intermediate files are being generated
* `htop` - monitor status of jobs and CPU usage (just google for details)
