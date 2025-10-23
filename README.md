# RNA-seq Analysis Pipeline with Snakemake
A comprehensive RNA-seq analysis pipeline for gene expresssion and alternative splicing analysis using a Snakemake workflow management.

**Overview** <br/>
This pipeline performs RNA-seq analysis from raw FASTQ files download to gene expression analysis with splicing information.
The work flow includes quality control, adapter trimming, alignment and quantifications <br/>

**Data used** <br/>
📊 Data from a published study PMID: 34632719, Bioproject PRJNA682755  <br/>
🏋🏾‍♀️ Data as a training set to develop the pipeline, which was tested to replicate results of this study <br/>

**Pipeline Features** <br/>
Quality Control: FastQC and MultiQC for QC check <br/>
Adapter Trimming: Adapter detection and removal <br/>
Genome Alignment: STAR aligner optimised for splice junction detection <br/>
Gene Quantification: FeatureCounts for gene level counting <br/>
Reproducibility: conda environment to ensure consistent software versions<br/>

**Requirements** <br/>
 🐍 Snakemake <br/>
 🐊 Conda/Mamba (for environment management)<br/>
 👩🏾‍💻 SLURM (cluster execution)<br/>

 **System requirements** <br/>
🧠 Memory: 86GB RAM (used, but could be more for human genome indexing)<br/>
🗄️ Storage: ∼ 100GB for genome indexing + analysis outputs <br/>
💾 CPU: Muliticore processor (used 16, but optimal is 8+)<br/>


**Setup** <br/>
1. Create directory, all the files would be here <br/>
```bash
###Create a project directory
mkdir rnaseq_analysis
cd rnaseq_analysis
```

2. Install snakemake
```bash
conda install -c conda-forge -c bioconda snakemake
```
3. Test the pipeline: Dry Run
   ```bash
   ##local machine
   snakemake -n 
   ```
