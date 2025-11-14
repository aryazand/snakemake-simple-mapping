rule freebayes:
    input:
        alns=get_bam_2,
        idxs=get_bai,
        ref=rules.get_genome.output.fasta,
        refidx=rules.get_genome.output.fai,
    output:
        bcf="results/freebayes/call/{sample}.bcf",
    params:
        extra=config["variant_calling"]["freebayes"]["extra"],
    log:
        "results/freebayes/call/{sample}.log",
    message:
        "call variants using freebayes"
    threads: 4
    wrapper:
        "v7.0.0/bio/freebayes"
