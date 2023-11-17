configfile: "./config/config.yaml"
num_threads = workflow.cores

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule all:
    input:
        "results/mapping/done",
        "results/parsing_done",
        "results/mapping/kma_GMS_summary.tsv",
        "results/mapping/kma_abfhpv_summary.tsv",
        expand("results/readstats/{sample}.tsv", sample=config["samples"]),
        "results/plots/SeqRunStats.png",
        "results/plots/heatmap_background_species.png",
        "results/plots/Time_and_AMR.png",
        "results/plots/ARG_exp_byclass.png"

rule get_readstats:
    input:
        fastq=get_bwa_map_input_fastqs
    output:
        stats_tsv="results/readstats/{sample}.tsv"
    conda:
        "FARMED_onsite.yml"
    shell:
        """
        python GetFastqStats.py {input.fastq} > {output.stats_tsv}
        """


rule filter_reads:
    input:
        fastq=get_bwa_map_input_fastqs
    output:
        fastq="results/fastq/{sample}.HQ.fastq"
    threads:
        num_threads
    log: "logs/{sample}.log"
    shell:
        """
        gunzip -c {input.fastq} | NanoFilt -q 7 -l 300 | gzip > {output.fastq.gz} 2>{log}
        """

rule map_reads_GMSdb:
    input:
        "results/fastq/{sample}.HQ.fastq.gz"
    output:
        sam=temp("results/mapping/{sample}_GMSdb.sam"),
        bam_unsorted=temp("results/mapping/{sample}_GMSdb.unsorted.bam"),
        bam="results/mapping/{sample}_GMSdb.bam",
        kmares="results/mapping/{sample}_GMSdb",
        kmaresfile="results/mapping/{sample}_GMSdb.res",
        kmamapstat="results/mapping/{sample}_GMSdb.mapstat"
    threads: num_threads
    log: "logs/{sample}.log"
    params:
        db=config['databases']['GMSdb'],
        kma="-mem_mode -bc 0.7 -bcNano -ID 0.0 -ef -proxi 0.9 -na -nf -nc -1t1 -ca -verbose 2",
        samview="-u -h -F0x004"
    shell:
        """
        kma -i {input} -o {output.kmares} -t_db {params.db} -tmp $(pwd)/ \
        {params.kma} -t {threads} -sam > {output.sam} 2>>{log} 
        
        samtools view {params.samview} -@{threads} -o {output.bam_unsorted} {output.sam} 2>>{log}
        samtools sort -@{threads} -o {output.bam} {output.bam_unsorted} 2>>{log}
        touch {output.kmares}
        echo "Done mapping"
        """

rule map_reads_ResFdb:
    input:
        "results/fastq/{sample}.HQ.fastq.gz"
    output:
        sam=temp("results/mapping/{sample}_ResFdb.sam"),
        bam_unsorted=temp("results/mapping/{sample}_ResFdb.unsorted.bam"),
        bam="results/mapping/{sample}_ResFdb.bam",
        kmares="results/mapping/{sample}_ResFdb",
        kmaresfile="results/mapping/{sample}_ResFdb.res",
        kmamapstat="results/mapping/{sample}_ResFdb.mapstat"
    threads: num_threads
    log: "logs/{sample}.log"
    params:
        db=config['databases']['ResFinder'],
        kma="-mem_mode -bc 0.7 -bcNano -ID 0.0 -ef -proxi 0.9 -na -nf -nc -ca -verbose 2",
        samview="-u -h -F0x004"
    shell:
        """
        kma -i {input} -o {output.kmares} -t_db {params.db} -tmp $(pwd)/ \
        {params.kma} -t {threads} -sam > {output.sam} 2>>{log}

        samtools view {params.samview} -@{threads} -o {output.bam_unsorted} {output.sam} 2>>{log}
        samtools sort -@{threads} -o {output.bam} {output.bam_unsorted} 2>>{log}
        touch {output.kmares}
        echo "Done mapping"
        """

rule map_reads_RefSeq:
    input:
        "results/fastq/{sample}.HQ.fastq.gz"
    output:
        sam=temp("results/mapping/{sample}_abfhpv.sam"),
        bam_unsorted=temp("results/mapping/{sample}_abfhpv.unsorted.bam"),
        bam="results/mapping/{sample}_abfhpv.bam",
        kmares="results/mapping/{sample}_abfhpv",
        kmaresfile="results/mapping/{sample}_abfhpv.res",
        kmamapstat="results/mapping/{sample}_abfhpv.mapstat"
    threads: num_threads
    log: "logs/{sample}.log"
    params:
        db=config['databases']['abfhpv'],
        kma="-mem_mode -bc 0.7 -bcNano -ID 0.0 -ef -proxi 0.9 -na -nf -nc -1t1 -ca -verbose 2",
        samview="-u -h -F0x004"
    shell:
        """
        kma -i {input} -o {output.kmares} -t_db {params.db} -tmp $(pwd)/ \
        {params.kma} -t {threads} -sam > {output.sam} 2>>{log}

        samtools view {params.samview} -@{threads} -o {output.bam_unsorted} {output.sam} 2>>{log}
        samtools sort -@{threads} -o {output.bam} {output.bam_unsorted} 2>>{log}
        touch {output.kmares}
        echo "Done mapping"
        """

rule finish_mapping:
    input:
        bamGMS=expand("results/mapping/{sample}_GMSdb.bam", sample=config["samples"]),
        bamResF=expand("results/mapping/{sample}_ResFdb.bam", sample=config["samples"]),
        bamabfhpv=expand("results/mapping/{sample}_abfhpv.bam", sample=config["samples"])
    output:
        done="results/mapping/done"
    shell:
        """
        touch {output.done}
        """

rule coverage_by_time_DMC:
    input:
        "results/mapping/{sample}_GMSdb.bam"
    output:
        tsv="results/CovByTime/{sample}_GMSdbByTime.tsv",
        plot="results/CovByTime/{sample}_GMSdbByTime"
    conda:
        "FARMED_onsite.yml"
    threads: num_threads
    log: "logs/{sample}.log"
    shell:
        """
        python GetCoverageByTime.py -o {output.plot} -t {threads} -r zymo_GMS.csv {input} > {output.tsv} 2>>{log}
        touch {output.plot}
        """

rule coverage_by_time_RefSeq:
    input:
        "results/mapping/{sample}_abfhpv.bam"
    output:
        tsv="results/CovByTime/{sample}_abfhpvByTime.tsv",
        plot="results/CovByTime/{sample}_abfhpvByTime"
    conda:
        "FARMED_onsite.yml"
    threads: num_threads
    log: "logs/{sample}.log"
    shell:
        """
        python GetCoverageByTime.py -o {output.plot} -t {threads} {input} > {output.tsv} 2>>{log}
        touch {output.plot}
        """

rule AMRlink_GMSdb:
    input:
        taxabam="results/mapping/{sample}_GMSdb.bam",
        resfbam="results/mapping/{sample}_ResFdb.bam"
    output:
        tsv="results/AMRlinking/{sample}_GMSvResF.tsv",
        plot="results/AMRlinking/{sample}_GMSvResF"
    conda:
        "FARMED_onsite.yml"
    threads: num_threads
    log: "logs/{sample}.log"
    shell:
        """
        python GetAMRlinksByTime.py -o {output.plot} -t {input.taxabam} -a {input.resfbam} > {output.tsv} 2>>{log}
        touch {output.plot}
        """

rule AMRlink_abfhpvdb:
    input:
        taxabam="results/mapping/{sample}_abfhpv.bam",
        resfbam="results/mapping/{sample}_ResFdb.bam"
    output:
        tsv="results/AMRlinking/{sample}_abfhpvvResF.tsv",
        plot="results/AMRlinking/{sample}_abfhpvvResF"
    conda:
        "FARMED_onsite.yml"
    threads: num_threads
    log: "logs/{sample}.log"
    shell:
        """
        python GetAMRlinksByTime.py -o {output.plot} -t {input.taxabam} -a {input.resfbam} > {output.tsv} 2>>{log}
        touch {output.plot}
        """

rule finish_alingment_parsing:
    input:
        CovByTimeTSV_GMSdb=expand("results/CovByTime/{sample}_GMSdbByTime.tsv", sample=config["samples"]),
        CovByTimeTSV_abfhpv=expand("results/CovByTime/{sample}_abfhpvByTime.tsv", sample=config["samples"]),
        AMRlinkTSV_GMSdb=expand("results/AMRlinking/{sample}_GMSvResF.tsv", sample=config["samples"]),
        AMRlinkTSV_abfhpv=expand("results/AMRlinking/{sample}_abfhpvvResF.tsv", sample=config["samples"])
    output:
        "results/parsing_done"
    shell:
        """
        touch {output}
        """

rule kmasummary_tsv:
    input:
        parsing_done="results/parsing_done",
        GMS=expand("results/mapping/{sample}_GMSdb.res", sample=config["samples"]),
        abfhpv=expand("results/mapping/{sample}_abfhpv.res", sample=config["samples"])
    output:
        GMS="results/mapping/kma_GMS_summary.tsv",
        abfhpv="results/mapping/kma_abfhpv_summary.tsv"
    conda:
        "FARMED_onsite.yml"
    log: "logs/kma_summary.log"
    shell:
        """
        python KMA_taxa_summary.py {input.GMS} > {output.GMS} 2>{log}
        python KMA_taxa_summary.py {input.abfhpv} > {output.abfhpv} 2>{log}
        """

rule plotseqrunstats:
    input:
        GMS="results/mapping/kma_GMS_summary.tsv",
        abfhpv="results/mapping/kma_abfhpv_summary.tsv",
        readstats=expand("results/readstats/{sample}.tsv", sample=config["samples"])
    output:
        seqrunstats="results/plots/SeqRunStats.png",
        bgheatmap="results/plots/heatmap_background_species.png"
    conda:
        "FARMED_onsite.yml"
    script:
        "Plot_SeqRunStats.R"

rule plot_time_amr_stats:
    input:
        CovByTimeTSV_GMSdb=expand("results/CovByTime/{sample}_GMSdbByTime.tsv", sample=config["samples"]),
        CovByTimeTSV_abfhpv=expand("results/CovByTime/{sample}_abfhpvByTime.tsv", sample=config["samples"]),
        AMRlinkTSV_GMSdb=expand("results/AMRlinking/{sample}_GMSvResF.tsv", sample=config["samples"]),
        AMRlinkTSV_abfhpv=expand("results/AMRlinking/{sample}_abfhpvvResF.tsv", sample=config["samples"])
    conda:
        "FARMED_onsite.yml"
    output:
        TimeAMR="/results/plots/Time_and_AMR.png",
        ARGbyclass="/results/plots/ARG_exp_byclass.png"
    script:
        "Species_ARGs_ByTime.R"
