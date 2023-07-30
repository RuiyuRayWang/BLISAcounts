rule wash_whitelist:
    input:
        whitelist="workflow/data/{user}/{project}/{library}/outs/{sample}_whitelist.txt",
        well_bc_ground_truth=config['plate']['well_settings'],
        ab_bc_ground_truth=config['plate']['ab_settings']
    output:
        "workflow/data/{user}/{project}/{library}/outs/{sample}_whitelist_washed.txt"
    script:
        "../scripts/wash_whitelist.py"

rule parse_bc_umi:
    input:
        "workflow/data/{user}/{project}/{library}/outs/{sample}_extracted.fastq.gz"
    output:
        temp("workflow/data/{user}/{project}/{library}/outs/{sample}_multiplex")
    threads:1
    shell:
        """
        zcat {input} | sed -n '1~4p' | awk 'BEGIN{{FS=\"_\";OFS=\"\\t\"}}{{print $2,$3}}' | \
        awk 'BEGIN{{FS=\" \";OFS=\"\\t\"}}{{print $1,$2}}' | sed 's/./&\\t/6' - > {output}
        """

rule extract_staggered:
    input:
        read="workflow/data/{user}/{project}/{library}/fastqs/{sample}.fastq.gz",
        well_bc_ground_truth=config['plate']['well_settings'],
        ab_bc_ground_truth=config['plate']['ab_settings']
    output:
        temp("workflow/data/{user}/{project}/{library}/outs/{sample}_staggered_multiplex")
    script:
        "../scripts/extract_staggered_seq.py"

rule demultiplex:
    input:
        get_demultiplex_inputs()
    output:
        "workflow/data/{user}/{project}/{library}/outs/{sample}_final_count"
    threads:1
    shell:
        """
        sort {input} | uniq | cut -f1,2 | sort | uniq -c | awk 'BEGIN{{FS=" ";OFS=\"\\t\"}}{{print $2,$3,$1}}' > {output}
        """
