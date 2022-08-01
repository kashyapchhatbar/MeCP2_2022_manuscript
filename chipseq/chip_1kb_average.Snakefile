bws, = glob_wildcards("data/GSE67293_DNA_mapping/deepTools_ChIP/bamCompare/{bw}.bw")
locis = ["1kb_subtract_CAn", "CAn_1kb", "CANn_1kb", "except_CAn_1kb"]

rule gabel:
    input: expand("window_averages/GSE67293/{bw}.mm9.{loci}.tsv", bw=bws, loci=locis)

rule gabel_average:
    input:
        b="data/GSE67293_DNA_mapping/deepTools_ChIP/bamCompare/{bw}.bw",
        w="windows/mm9_{loci}.bed"    
    output: "window_averages/GSE67293/{bw}.mm9.{loci}.tsv"
    shell: "bigWigAverageOverBed {input.b} {input.w} {output}" 

boxer_bws, = glob_wildcards("data/GSE139509_DNA_mapping/deepTools_ChIP/bamCompare/{boxer_bw}.bw")

rule boxer:
    input: expand("window_averages/GSE139509/{boxer_bw}.mm10.{loci}.tsv", boxer_bw=boxer_bws, loci=locis)

rule boxer_average:
    input:
        b="data/GSE139509_DNA_mapping/deepTools_ChIP/bamCompare/{boxer_bw}.bw",
        w="windows/mm10_{loci}.bed"    
    output: "window_averages/GSE139509/{boxer_bw}.mm10.{loci}.tsv"
    shell: "bigWigAverageOverBed {input.b} {input.w} {output}" 

rule liftover:
    input: expand("data/GSE139509_DNA_mapping/deepTools_ChIP/bamCompare/{boxer_bw}.mm9.bw", boxer_bw=boxer_bws)

rule boxer_average_liftover:
    input:
        b="data/GSE139509_DNA_mapping/deepTools_ChIP/bamCompare/{boxer_bw}.bw",
    params:
        c="../methods/windows/mm10ToMm9.over.chain",
        s="../methods/data/mm9.chrom.sizes"    
    output:
        b=temp("data/GSE139509_DNA_mapping/deepTools_ChIP/bamCompare/{boxer_bw}.mm9.unsorted.bedGraph"),
        s=temp("data/GSE139509_DNA_mapping/deepTools_ChIP/bamCompare/{boxer_bw}.mm9.sorted.bedGraph"),
        m=temp("data/GSE139509_DNA_mapping/deepTools_ChIP/bamCompare/{boxer_bw}.mm9.sorted.merged.bedGraph"),
        bw="data/GSE139509_DNA_mapping/deepTools_ChIP/bamCompare/{boxer_bw}.mm9.bw" 
    run:
        shell("liftOver <(bigWigToBedGraph {input.b} /dev/stdout) {params.c} {output.b} /tmp/unmapped")
        shell("sort -k1,1 -k2,2n {output.b} > {output.s}")
        shell("bedtools merge -i {output.s} -d -1 -c 4 -o mean > {output.m}")
        shell("bedGraphToBigWig {output.m} {params.s} {output.bw}")


chen_bws, = glob_wildcards("data/GSE66868_DNA_mapping/deepTools_ChIP/bamCompare/{chen_bw}.bw")

rule chen:
    input: expand("window_averages/GSE66868/{chen_bw}.mm9.{loci}.tsv", chen_bw=chen_bws, loci=locis)

rule chen_average:
    input:
        b="data/GSE66868_DNA_mapping/deepTools_ChIP/bamCompare/{chen_bw}.bw",
        w="windows/mm9_{loci}.bed"    
    output: "window_averages/GSE66868/{chen_bw}.mm9.{loci}.tsv"
    shell: "bigWigAverageOverBed {input.b} {input.w} {output}"
    