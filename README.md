## Data analysis code for 'The nucleotide messenger (p)ppGpp is an anti-inducer of the purine synthesis transcription regulator PurR in Bacillus' ##
## Introduction ##
This directory contains the information on how reads were handled. Raw data are available at GEO accession GSE185164. Have fun!  

## Pre-processing of ChIP data ##
For each fastq file, the following steps were taken to preprocess the data

Reads were concatenated from multiple files into a single large fastq.gz file

Here:  
arg0 is substituted for a samples prefix

Change directories into the directory containing the reads to be concatenated,
then run the following:  

```bash
cat *R1*.fastq.gz > arg0_all_R1.fastq.gz
```

QC of original data was conducted. In each directory containing concatenated reads, run:

```bash
mkdir fastqc_before
fastqc arg0_all_R1.fastq.gz -f fastq -o fastqc_before/
```

Illumina adapters were trimmed from reads.

Here:
arg0 is substituted for a samples prefix
arg1 represents the fastq for R1 of the paired end sequencing

We found that after removing adapters starting with "AGATCGGAAGAGC" we were still left with many reads with adapter starting with "CGGAAGAGCACAC". We also have many reads with "GGGGGGGGGGGGGGGGGGGGGG...", so we remove them too.  
NOTE: make sure you adjust the -threads argument to trimmomatic to be compatible with your processor's number of threads.  

In each directy with concatenated reads, run:

```bash
cutadapt --quality-base=33 -a AGATCGGAAGAGC -a CGGAAGAGCACAC -a GGGGGGGGGGGGGGG -n 3 -m 20 --match-read-wildcards -o arg0_R1_trim_unpaired.fastq.gz arg1 > arg0_cutadapt.log 2> arg0_cutadapt.err
```

To QC trimmed files, in the directory with trimmed files, run:

```bash
mkdir fastqc_after
fastqc arg0_R1_trim_unpaired.fastq.gz -f fastq -o fastqc_after/
```

Change directories into the directory containing your concatenated reads,
then run `mkdir original_data`, then run `mv *.fastq.gz original_data` to move
your reads into the original_data folder.  
Now run `cd original_data` to enter your data directory.  
Now run:  

```bash
mkdir fastqc_before
mkdir fastqc_after
```

You are now ready to do your QC of your reads.  

```bash
fastqc arg1 -f fastq -o fastqc_before/
cutadapt --quality-base=33 -a AGATCGGAAGAGC -a CGGAAGAGCACAC -a GGGGGGGGGGGGGGG -n 3 -m 20 --mask-adapter --match-read-wildcards -o arg0_R1_cutadapt.fastq.gz arg1 > arg0_cutadapt.log 2> arg0_cutadapt.err
java -jar /path/to/trimmomatic-0.39.jar SE -threads 10 -phred33 arg0_R1_cutadapt.fastq.gz arg0_R1_trim_unpaired.fastq.gz TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 > arg0_trim.log 2> arg0_trim.err
fastqc arg0_R1_trim_paired.fastq.gz -f fastq -o fastqc_after/
```

HINT: you can loop through your sample directories.  
For example:  

```bash
for dir in Sample_*; do (cd $dir; cd original_data &&  echo $dir"_all_R1.fastq.gz"; cutadapt --quality-base=33 -a AGATCGGAAGAGC -a CGGAAGAGCACAC -a GGGGGGGGGGGGGGG -n 3 -m 20 --mask-adapter --match-read-wildcards -o $dir"_R1_cutadapt.fastq.gz" $dir"_all_R1.fastq.gz" > $dir"_cutadapt.log" 2> $dir"_cutadapt.err"); done)
```

In the above line of code, each time through the loop, the directory's name is prepended to your read file suffix using `$dir`.  
The `echo $dir"_all_R1.fastq.gz"` just serves the role of reporting to you which file your loop is working on at a given moment.  
This is useful, since it lets you know your loop isn't just doing nothing.  

## Alignment of processed reads ##
Next, each pair of processed reads was aligned with bowtie2.

Here  
arg0 is substituted for a samples prefix  
arg1 represents the trimmed fastq for R1 

For aligning p0 reads to NCIB 3610 genome (no pBS32)

```bash
bowtie2 -x /home/wanglab/Users_local/Jeremy/Sequencing/RefGenomes/NCIB_3610_CP020102 \
        -U arg1 -q --end-to-end --very-sensitive -p 8 \
        --phred33 2> arg0_bow.log | samtools view -bSh - | samtools sort -@ 8 -o arg0_sort.bam -
samtools index arg0_sort.bam
samtools view -bh -F 2820 -q 30 arg0_sort.bam > arg0_sort_filtered.bam
```

For aligning WT reads to NCIB 3610 genome + pBS32 plasmid:

```bash
bowtie2 -x /mnt/jadelab/lab/current/NGS/RefGenomes/NCIB_3610_CP020102_pBS32_CP020103 -U arg1 -q --end-to-end --very-sensitive -p 6 --phred33 2> arg0_bow.log | samtools view -bSh - | samtools sort -@ 8 -o arg0_sort.bam -
samtools index arg0_sort.bam
samtools view -bh -F 2820 -q 30 arg0_sort.bam > arg0_sort_filtered.bam
samtools index arg0_sort_filtered.bam
```

When executed in directory containing all Sample_* directories, this loop allows quick indexing of *_sort_filtered.bam files:

```bash
for dir in Sample_*; do (cd $dir && samtools index $dir"_sort_filtered.bam"); done
```

## Mapping of coverage and preparation for bootstrapping ##
Next, reads were filtered with samtools and made into a sampling object using
the custom script `bootstrap_sam_file.py`'s `parse` option.
'bootstrap_sam_file.py' is available at https://github.com/mikewolfe/2018_Lrp_ChIP.git

Here  
arg0 is the sample prefix  
arg1 is the input bam file  
To understand the flags indicated by the -F option passed to samtools, go to this [website](https://broadinstitute.github.io/picard/explain-flags.html "SAM flags explained").  

```bash
samtools view arg0_sort_filtered.bam CP020102 | python3 ~/src/2018_Lrp_ChIP/ChIP_analysis/bootstrap_sam_file.py parse - arg0.ob 2> arg0_sampler.err 
```

When executed in directory containing all Sample_* directories, this loop will accomplish this for each sample:

```bash
for dir in Sample_1290*; do (cd $dir && samtools view $dir"_sort_filtered.bam" CP020102 | python3 ~/src/2018_Lrp_ChIP/ChIP_analysis/bootstrap_sam_file.py parse - $dir".ob" 2> $dir"_sampler.err"); done
```

## Obtaining summary tracks at 10 bp resolution

To obtain bootstrapped read counts we will use `bootstrap_sam_file.py sample`.
'bootstrap_sam_file.py' is available at https://github.com/mikewolfe/2018_Lrp_ChIP.git

```bash
for dir in Sample_1290*;
    do (cd $dir && python3 ~/src/2018_Lrp_ChIP/ChIP_analysis/bootstrap_sam_file.py sample $dir".ob.npy" $dir"_sampled" --reference_genome /PATH --num_samples 100 --resolution 10 --seed 1234 2> $dir"_bootstrap.err");
done
```

Now we'll summarize the bootstrapped samples by calculating counts per million at each position for each bootstrap, and divide ChIP by input to get enrichment

```bash
python3 summarize_all_sample_counts.py
```

## Smoothing counts per million

To smooth, use 'NGS_smoothing.py' with 50bp windows
Exported as 'smoothed_data.csv'

## Calculating enrichment, with either smoothed or unsmoothed data

### For unsmoothed data:
Enrichment was calculated in R ('PurR_analysis.R' file in the PurR_ChIP_analysis repository).
Briefly, for unsmoothed bootstrapped data, the median of the bootstraps at each 10bp window was determined.

A csv was created with information on each sample: PurR_chip_and_input_pairs.csv
This csv also contained information on whether each sample was ChIP or input as well as
which ChIP samples were paired with which input samples.

Using the function get_ChIP_enrichment found in 'helperFunctions.R', ChIP enrichment relative to input
was determined by dividing the median cpm in ChIP by the median cpm in input at each 10bp window.

The enrichment was converted to log2 and infinite values were excluded.
The final dataset was saved as 'summary_df.RData'

### For smoothed data:
The csv containing the smoothed data was loaded into R.
Various manipulations were necessary to get the dataset into the proper format.
Otherwise, the enrichment was calculated as above.
The function get_ChIP_enrichment in 'helperFunctions.R' was used to calculate enrichment just as above.


##Obtaining ChIP peaks
To obtain ChIP peaks, the irreproducible discovery rate (IDR) was first calculated for the log2 enrichment of each sample. Log2 enrichment for each sample was positioned in a n by m dataframe (n rows for each 10bp summary track and m columns for each replicate). IDR was calculated using est.IDR function in calculate_idr.R. IDR was saved as an npy array.

Second, *_summary_mad.npy was obtained for each sample be rerunning summarize_all_sample_counts.py with an updated bootstrap_sam_file.py that added 'MAD' files for each sample's summary output.

## Obtaining replicate summary statistics for enrichment of ChIP over input ##

To obtain RE and RSE summary signals we used the samplers from the previous step
as input to the `bootstrapped_chip_enrichment.py` script.

Here
run_info is the csv provided by the UofM sequencing core with Sample_IDs and Descriptions
arg0 is the sample prefix
arg1 and arg2 are samplers for the WT extracted samples
arg3 and arg4 are samplers for the WT input samples

Note that arg1 and arg3 are paired, arg2 and arg4 are paired etc.

```bash
python3 ~/src/2018_Lrp_ChIP/ChIP_analysis/bootstrapped_chip_enrichment.py --sample_name_lut run_info --genome_size 4215607 --out_prefix arg0 --ChIP_samps arg1 arg2 --inp_samps arg3 arg4 --num_replicates 1 --identity -s 1234 -p 8 --save_summaries 0.05 --resolution 10 2> arg0.log
```

## Obtaining bootstrap replicate summary statistics ##
To obtain the bootstrap MAD stat (as well as additional statistics) from 1000
bootstrap replicates for each 10 bp location in the genome we also used the
`bootstrapped_chip_enrichment.py` script.

Here
arg0 is the sample prefix
arg1 is the sampler for the ChIP sample ('Sample_*.ob.npy')
arg2 is the sampler for the corresponding input sample

Note that arg1 and arg2 are paired

```bash
python3 bootstrapped_chip_enrichment.py --sample_name_lut run_info --genome_size 4215607 --ChIP_samps arg1 --inp_samps arg2 --num_replicates 1000 -s 1234 -p 8 --save_summaries 0.05 --resolution 10 --numba 2> arg0.log
```


### Determining the effect of ppGpp on ChIP signal ###
Overall goal: Determine enrichment for bootstrapped coverage for each WT and p0 sample (before and after RHX). Then use this enrichment for each biological replicate (rep1, rep2, rep3) to take the ratio (wt_post / wt_pre)/(p0_post / p0_pre) to obtain an array for the effect of ppGpp on the ChIP signal after RHX treatment. Then take the median of the array for each biological replicate. Finally, calculate the mean and standard deviation of the biological replicates.

All commands run with 'calculate_effect_of_ppGpp.py' script.

To determine the effect of ppGpp on ChIP signal, bootstrapped .npy arrays of ChIP coverage for taken for each sample (file = 'Sample_X_sampled.npy'). Effect of ppGpp was determined with (enrich_wt_post / enrich_wt_pre) / (enrich_p0_post / enrich_p0_pre). The median bootstrap value was determined for each 'effect_of_ppGpp' biological replicate. Means and standard deviations were calculated from the medians of each biological replicate. The output files were "mean_ppGpp_effect_linear.npy" and "sd_ppGpp_effect_linear.npy". These files were used in the "ppGpp_effect.R" script for plotting.
