#TEA
Welcome to the companion site for the Tea repeat analysis pipeline. Tea stands for **T**ransposable **E**lement **A**nalyzer. This site will serve as a source code repository for the Tea pipeline. If you would like to dowonload repeat library files and additional data related to the paper [Visit TEA website!](http://compbio.med.harvard.edu/Tea/)
# Installing TEA
Installation requires the following software to be present on the system: 

1. BWA
 * BWA 0.6.2 is required. 
2. Samtools
 * Samtools 0.1.5 or above. 
3. Meerkat
 * Download Meetkat and see its intall instruction at [Meerkat website](http://compbio.med.harvard.edu/Meerkat/). 
  * Make sure that Meerkat is correctly installed.  
4. TEA
 * Source code for the Tea pipeline at [GitHub](https://github.com/hastj7373/TEA).
 * CAP3 assembler.
 * R packages : Bioconductor, Rsamtools, IRanges and [spp](http://compbio.med.harvard.edu/Supplements/ChIP-seq/).
 * Download repeat library at [TEA website](http://compbio.med.harvard.edu/Tea/)(Donwload repeat.combined.div30.isize150.fa and move it to /TEA_installed/lib/assembly). 

# Running TEA

1. Fill out **conf_file**
 * Example
 ```
samtools_home = /cluster/data/scratch/share/samtools-0.1.19
bwa_home = /cluster/data/scratch/share/bwa-0.6.2
blast_home = /cluster/bio/bin
cap3_home = /cluster/data/scratch/share/CAP3/
fasta_location = /cluster/data/scratch/db/hg19/hg19.fasta
meerkat_home= /cluster/data/scratch/share/Meerkat/scripts
tea_home = /home/hcjung/Tea
tea_run = scripts/tea.pl
tea_assembly = lib/assembly
tmp_dir = /cluster/data/scratch/TEA_COAD/tmp
retro = repeat.combined.div30.isize150.fa
```
2. Copy of **conf_file** and **tea_run** into a directory with bam files that you want to run. 
 * Example
 
  ![Image of secondstep](https://github.com/hastj7373/TEA/blob/master/second_step.gif)

3. Execute **tea_run** 
 ```
 tea.run bam_input tea ra -c hg19 -p 5 -f -K
 ```
 In the current version, users need to input a bam file and set # of threads for the parameter *-p* (Keep other parameters like above). Note that users need to input a bam file without ".bam" (see the example below)
  * Example (see details by clicking 'View Raw' in [here](https://github.com/hastj7373/TEA/blob/master/Example.docx))
  
  ![Image of thridstep](https://github.com/hastj7373/TEA/blob/master/third_step.gif)
 

4. Output Results

TEA produces sample_name.germline.contig
