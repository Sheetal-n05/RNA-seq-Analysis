# RNA-seq-Analysis


**1. Experimental Setup**

The process involves breaking down DNA or cDNA samples into random fragments, which are then joined with adapter molecules at both ends. These fragments with adapters are either amplified or cleaned up for further use

Library is loaded into a flow cell; fragments are captured on a lawn of oligos complementary to the adapters. 
Each fragment is then amplified into distinct, clonal clusters

Illumina technology uses a reversible terminator–based method that detects single bases.
All four reversible terminator–bound dNTPs are present during each sequencing cycle, minimizing incorporation bias and error rates.

**2. Pipeline processing**

Post sequencing, raw data (.fq) processed through pipeline in AWS. Processing gets data analysis ready. 

Here, we can view MultiQC reports for all the samples and verify whether the pipeline has run properly through all of its components.


**3. Data Analysis and visualization**

During data analysis and alignment, the newly identified sequence reads are aligned to a reference genome following alignment, many variations of analysis are possible.

**With analysis-ready dataset, can perform:**
- Expression based analysis
- Sample similarity
- Differential expression
- Cell type identity
