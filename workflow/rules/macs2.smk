rule macs2_callpeaks:
    input:
        treatment=get_bam_2,
    output:
        peaks="results/macs2/{sample}_peaks.narrowPeak",
        summits="results/macs2/{sample}_summits.bed",
        xls="results/macs2/{sample}_peaks.xls",
    log:
        "results/macs2/{sample}.log",
    message:
        "call peaks with MACS2",
    params:
        extra=config["peak_calling"]["macs2"]["extra"]
    threads: 2
    wrapper:
        "v2.9.1/bio/macs2/callpeak"