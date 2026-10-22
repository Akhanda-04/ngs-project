# Pre-trim vs Post-trim Quality Comparison

**Dataset:** ERR008613  
**Purpose:** Evaluate sequencing-read quality before and after Trimmomatic trimming.

## 1. Adapter Content

The post-trim FastQC reports indicated zero detectable adapter content in both R1 and R2. This suggests that adapter removal was effective. The pre-trim and post-trim Adapter Content plots should be compared to confirm the change.

## 2. Per-base Sequence Quality

The per-base sequence quality profiles improved after trimming in both R1 and R2.

The overall base-quality statistics also improved:

- **R1:** Q20 increased from 89.86% to 98.64%, while Q30 increased from 82.24% to 92.42%.
- **R2:** Q20 increased from 87.79% to 98.00%, while Q30 increased from 79.37% to 90.83%.

Average base quality increased from Q13.06 to Q29.41 for R1 and from Q12.35 to Q27.80 for R2. These results indicate substantially improved quality among the retained reads.

## 3. Read Retention

The original dataset contained 14,214,324 read pairs.

After trimming, 12,137,629 reads remained in each paired FASTQ file.

- Original read pairs: 14,214,324
- Retained paired reads: 12,137,629
- Pairs not retained as paired reads: 2,076,695
- Paired-read retention rate: 85.39%
- Percentage not retained as paired reads: 14.61%

The R1 and R2 paired FASTQ files contain equal numbers of reads, confirming consistent paired-file counts.

## 4. Read Length

Before trimming, both R1 and R2 had an average read length of 100 bp.

After trimming:

- **R1:** Average read length decreased to 92.7 bp, a reduction of 7.3 bp (7.3%).
- **R2:** Average read length decreased to 91.4 bp, a reduction of 8.6 bp (8.6%).

The minimum retained read length was 50 bp for both R1 and R2, while the maximum remained 100 bp.

These results indicate that trimming shortened the reads moderately while retaining most of the original read pairs.

## 5. GC Distribution and Overrepresented Sequences

The average GC content changed as follows:

- **R1:** 50.76% before trimming and 50.02% after trimming.
- **R2:** 50.84% before trimming and 49.94% after trimming.

The overall GC percentages remained close to 50%. However, the FastQC Per Sequence GC Content plots should be examined to determine whether the distribution shapes are unusual.

No overrepresented sequences were reported in the post-trim FastQC results. This finding is encouraging but does not independently rule out contamination.

If unusual GC distributions persist, further investigation using an appropriate contamination-screening or taxonomic-classification tool may be necessary.

## Overall Conclusion

Trimming retained 85.39% of the original read pairs and improved the overall base-quality statistics in both read directions. The average read length decreased by 7.3 bp for R1 and 8.6 bp for R2. Both paired FASTQ files contain 12,137,629 reads.

The reported adapter-content and overrepresented-sequence results suggest effective removal of detectable adapter-related sequences. The GC-content distribution should be reviewed further to determine whether the observed patterns are biologically expected or indicate possible contamination.

The trimmed paired-end reads will be considered for downstream bacterial genome assembly after the remaining quality-control checks are reviewed.

## Data Source

Read statistics were calculated using `seqkit stats -a` on the raw and trimmed FASTQ files. Quality metrics and adapter-content observations were obtained from the corresponding FastQC reports.

## Tools

- FastQC: Read-quality assessment
- Trimmomatic: Adapter removal and read trimming
- SeqKit: FASTQ statistics and read-length analysis
