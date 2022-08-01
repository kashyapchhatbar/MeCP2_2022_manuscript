locis = ["1kb_subtract_CAn", "CAn_1kb", "CANn_1kb", "except_CAn_1kb"]

rule mC_average:
    input:
        expand("window_averages/wgbs/{sample}.mm9.{loci}.tsv", sample=["lagger_CAC", "lister_CAC"], loci=locis)

rule run_bedtools_map:
    input:
        b = "wgbsseq/{sample}.bedGraph.gz",
        w = "wgbsseq/windows/mm9_{loci}.bed"
    output: "window_averages/wgbs/{sample}.mm9.{loci}.tsv"
    shell: "bedtools map -a {input.w} -b {input.b} -c 4 -o sum > {output}"
 
rule mC_average_boxer:
    input:
        expand("window_averages/wgbs/boxer_CAC.mm10.{r}.tsv", r=locis)

rule boxer_bedtools_map:
    input:
        bed = "wgbsseq/boxer_CAC.bedGraph.gz",
        win = "wgbsseq/windows/mm10_{r}.bed"
    output: "window_averages/wgbs/boxer_CAC.mm10.{r}.tsv"
    shell: "bedtools map -a {input.win} -b {input.bed} -c 4 -o sum > {output}"
    
rule screenshot_bw:
    input:
        expand("window_averages/wgbs/{sample}.mm9.1kb.bw", sample=["lagger_CAC", "lister_CAC"]),
        expand("window_averages/wgbs/{sample}.mm10.1kb.bw", sample=["boxer_CAC"])

rule run_screenshot_bw:
    input:
        b = "wgbsseq/{sample}.bedGraph.gz",
        s = "../methods/data/mm9.chrom.sizes"
    output:
        final="window_averages/wgbs/{sample}.mm9.1kb.bw",
        intermediate=temp("window_averages/wgbs/{sample}.mm9.1kb.bedGraph")
    run:
        shell("bedtools makewindows -g {input.s} -w 1000 | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.b} -c 4 -o sum -null 0 > {output.intermediate}")
        shell("bedGraphToBigWig {output.intermediate} {input.s} {output.final}")

rule run_screenshot_bw_boxer:
    input:
        b = "wgbsseq/{sample}.bedGraph.gz",
        s = "../methods/data/mm10.chrom.sizes"
    params:
        
    output: 
        final="window_averages/wgbs/{sample}.mm10.1kb.bw",
        intermediate=temp("window_averages/wgbs/{sample}.mm10.1kb.bedGraph")
    run:
        shell("bedtools makewindows -g {input.s} -w 1000 | sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.b} -c 4 -o sum -null 0 > {output.intermediate}")
        shell("bedGraphToBigWig {output.intermediate} {input.s} {output.final}")

rule liftover:
    output:
        bw="window_averages/wgbs/boxer_CAC.mm9.1kb.bw",
        b=temp("window_averages/wgbs/boxer_CAC.mm9.1kb.unsorted.bedGraph"),
        s=temp("window_averages/wgbs/boxer_CAC.mm9.1kb.sorted.bedGraph"),
        m=temp("window_averages/wgbs/boxer_CAC.mm9.1kb.sorted.merged.bedGraph")
    input:
        b="window_averages/wgbs/boxer_CAC.mm10.1kb.bw"        
    params:
        c="../methods/windows/mm10ToMm9.over.chain",
        s="../methods/data/mm9.chrom.sizes"
    run:
        shell("liftOver <(bigWigToBedGraph {input.b} /dev/stdout) {params.c} {output.b} /tmp/unmapped")
        shell("sort -k1,1 -k2,2n {output.b} > {output.s}")
        shell("bedtools merge -i {output.s} -d -1 -c 4 -o mean > {output.m}")
        shell("bedGraphToBigWig {output.m} {params.s} {output.bw}")

rule merge_CAn:
    input: expand("wgbsseq/windows/{assembly}.CAn.bed", assembly=["mm9", "mm10"])

rule run_merge_CAn:
    input: "../methods/data/{assembly}.CAn.bed"
    output: "wgbsseq/windows/{assembly}.CAn.bed"
    shell: "bedtools merge -i {input} > {output}"
    