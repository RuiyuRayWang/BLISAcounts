rule umi_tools_whitelist:
    input:
        "workflow/data/{user}/{project}/{library}/fastqs/{sample}.fq.gz"
    output:
        "workflow/data/{user}/{project}/{library}/outs/{sample}_whitelist.txt"
    log:
        "workflow/data/{user}/{project}/{library}/logs/{sample}_whitelist.log"
    threads:1
    shell:
        """
        umi_tools whitelist --extract-method=regex \
                            --bc-pattern='(?P<cell_1>.{{6}})CTTGTGGAAAGGACGAAACA{{s<=2}}(?P<cell_2>.{{6}})(?P<umi_1>.{{15}}).*' \
                            -L {log} \
                            --stdin {input} \
                            --set-cell-number=800 \
                            --plot-prefix=workflow/data/{wildcards.user}/{wildcards.project}/{wildcards.library}/outs/{wildcards.sample} \
                            --log2stderr > {output}
        """

rule umi_tools_extract:
    input:
        fq="workflow/data/{user}/{project}/{library}/fastqs/{sample}.fq.gz",
        whitelist="workflow/data/{user}/{project}/{library}/outs/{sample}_whitelist_washed.txt"
    output:
        "workflow/data/{user}/{project}/{library}/outs/{sample}_extracted.fastq.gz"
    log:
        "workflow/data/{user}/{project}/{library}/logs/{sample}_extract.log"
    threads:1
    shell:
        """
        umi_tools extract --extract-method=regex \
                          --bc-pattern='(?P<cell_1>.{{6}})CTTGTGGAAAGGACGAAACA{{s<=2}}(?P<cell_2>.{{6}})(?P<umi_1>.{{15}}).*' \
                          -L {log} \
                          --stdin {input.fq} \
                          --error-correct-cell \
                          --filter-cell-barcode \
                          --stdout {output} \
                          --whitelist={input.whitelist}
        """