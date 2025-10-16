sample_id = ["SRR13202556",
"SRR13202557",
"SRR13202558",
"SRR13202559",
"SRR13202560",
"SRR13202561",
"SRR13202562",
"SRR13202563",
"SRR13202564",
"SRR13202565",
"SRR13202566",
"SRR13202567"
]

rule all: 
    input:
        expand("QC_results/fastqc/{accession}/{accession}_{pair}_fastqc.html", accession = sample_id, pair = (1, 2)),
        expand("QC_results/fastqc/{accession}/{accession}_{pair}_fastqc.zip",accession = sample_id, pair=(1,2)),
        expand("results/star/{accession}/Aligned.sortedByCoord.out.bam", accession = sample_id),
        expand("results/star/{accession}/Aligned.sortedByCoord.out.bam.bai", accession = sample_id),
        expand("results/star/{accession}/ReadsPerGene.out.tab", accession = sample_id),
        expand("results/star/{accession}/SJ.out.tab", accession = sample_id),
        expand("results/featurecounts/{accession}_counts.txt", accession = sample_id),
        "results/gene_expression/combined_counts_matrix.txt",
        "results/gene_expression/combined_counts_matrix_tpm.txt",
        "QC_results/multiqc/multiqc_report.html",
        "genome_index"

# New rule to download SRA data
rule download_sra:
    conda: "envs/qc.yml"
    output:
        fastq_1 = "rawreads/{accession}_1.fastq",
        fastq_2 = "rawreads/{accession}_2.fastq"
    params:
        outdir = "rawreads"
    log:
        out = "log/download_sra/{accession}.out",
        err = "log/download_sra/{accession}.err"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p log/download_sra
        fasterq-dump --split-files {wildcards.accession} -O {params.outdir} >{log.out} 2>{log.err}
        """

rule fastqc:
    conda: "envs/qc.yml"
    input:
        fastq_1 = "rawreads/{accession}_1.fastq",
        fastq_2 =  "rawreads/{accession}_2.fastq"
    output:
        folder = directory("QC_results/fastqc/{accession}"),
        fastq_1 = "QC_results/fastqc/{accession}/{accession}_1_fastqc.zip",
        fastq_2 = "QC_results/fastqc/{accession}/{accession}_2_fastqc.zip",
        report_1 = "QC_results/fastqc/{accession}/{accession}_1_fastqc.html",
        report_2 = "QC_results/fastqc/{accession}/{accession}_2_fastqc.html",
    log:
        out = "log/fastqc/{accession}.out",
        err = "log/fastqc/{accession}.err",
    shell:
        """
        mkdir -p {output.folder}
        mkdir -p log/fastqc
        fastqc {input.fastq_1} {input.fastq_2} -o {output.folder} >{log.out} 2>{log.err}
        """

rule multiqc:
    conda: "envs/qc.yml"
    input:
        expand("QC_results/fastqc/{accession}/{accession}_{pair}_fastqc.html", accession=sample_id, pair=(1, 2))
    output:
        "QC_results/multiqc/multiqc_report.html"
    log: "log/multiqc.log"
    shell:
        """
        mkdir -p QC_results/multiqc
        mkdir -p log
        multiqc QC_results/fastqc/ -o QC_results/multiqc >{log} 2>&1
        """

rule download_reference:
    output:
        fasta = "reference/human.fna.gz"
    shell: 
        """
        mkdir -p reference
        curl 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz' -o {output.fasta}
        """

rule download_gtf:
    output:
        gtf = "reference/human.gtf.gz"
    shell: 
        """
        mkdir -p reference
        curl 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz' -o {output.gtf}
        """

rule decompress_fna:
    input: "reference/human.fna.gz"
    output: "reference/human.fna"
    shell:
        "gzip -d -c {input} >{output}"

use rule decompress_fna as decompress_gtf with:
    input: "reference/human.gtf.gz"
    output: "reference/human.gtf"

rule index_genome:
    conda: "envs/star.yaml"
    threads: 8
    input:
        fasta = "reference/human.fna",
        gtf = "reference/human.gtf"
    output:
        directory("genome_index")
    params:
        sjdbOverhang = 149
    log: "log/index_genome.log"
    shell:
        """
        mkdir -p {output}
        mkdir -p log
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang} \
             --genomeSAindexNbases 14 \
             --limitGenomeGenerateRAM 31000000000 >{log} 2>&1
        """

rule star_align:
    conda: "envs/star.yaml"
    shadow: "shallow"
    threads: 8
    input:
        index = "genome_index",
        gtf = "reference/human.gtf",
        fastq_1 = "rawreads/{accession}_1.fastq",
        fastq_2 =  "rawreads/{accession}_2.fastq"
    output: 
        bam = "results/star/{accession}/Aligned.sortedByCoord.out.bam",
        splice_junctions = "results/star/{accession}/SJ.out.tab",
        counts = "results/star/{accession}/ReadsPerGene.out.tab",
        chimeric = "results/star/{accession}/Chimeric.out.junction"
    params:
        outdir = "results/star/{accession}",
        prefix = "results/star/{accession}/",
        sjdb_overhang=149,
        align_sjdb_overhang_min=3,
        align_sj_stitch_mismatch_nmax="5 -1 5 5",
        out_filter_mismatch_nmax=2,
        out_filter_multimap_nmax=10,
        align_intron_max=500000,
        align_mates_gap_max=500000
    log: "log/star_align/{accession}.log"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p log/star_align
        STAR --runMode alignReads \
            --genomeDir {input.index} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdb_overhang} \
            --readFilesIn {input.fastq_1} {input.fastq_2} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --runThreadN {threads} \
            --outSAMattributes All \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
            --alignSJDBoverhangMin {params.align_sjdb_overhang_min} \
            --alignSJstitchMismatchNmax {params.align_sj_stitch_mismatch_nmax} \
            --outFilterMismatchNmax {params.out_filter_mismatch_nmax} \
            --outFilterMultimapNmax {params.out_filter_multimap_nmax} \
            --alignIntronMax {params.align_intron_max} \
            --alignMatesGapMax {params.align_mates_gap_max} \
            --twopassMode Basic \
            --chimOutType Junctions SeparateSAMold \
            --chimSegmentMin 20 \
            --chimJunctionOverhangMin 20 \
            --chimOutJunctionFormat 1 \
            --outSAMunmapped Within \
            --outSAMmapqUnique 60 \
            --outMultimapperOrder Random \
            --runRNGseed 12345 >{log} 2>&1
        """

rule index_bam:
    conda: "envs/indexing.yaml"
    input:
        bam = "results/star/{accession}/Aligned.sortedByCoord.out.bam"
    output:
        bai = "results/star/{accession}/Aligned.sortedByCoord.out.bam.bai"
    log:
        "log/index_bam/{accession}.log"
    shell:
        """
        mkdir -p log/index_bam
        samtools index {input.bam} > {log} 2>&1
        """

rule featurecounts:
    conda: "envs/featurecounts.yaml"
    input:
        bam = "results/star/{accession}/Aligned.sortedByCoord.out.bam",
        gtf = "reference/human.gtf"
    output:
        counts = "results/featurecounts/{accession}_counts.txt",
        summary = "results/featurecounts/{accession}_counts.txt.summary"
    threads: 4
    log: "log/featurecounts/{accession}.log"
    shell:
        """
        featureCounts -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -g gene_id \
            -t exon \
            -s 0 \
            -p \
            --countReadPairs \
            -B \
            -C \
            {input.bam} >{log} 2>&1
        """

