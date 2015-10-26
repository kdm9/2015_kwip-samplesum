GENOMES = ["A", "B", "C", "D"]
SAMPLES = [str(i) for i in range(1, 5)]
GENOME_SEEDS = {g: str(i) for i, g in enumerate(GENOMES)}
SAMPLE_SEEDS = {g: {s: str(i ** j) for j, s in enumerate(SAMPLES)}
                for i, g in enumerate(GENOMES)}
READ_NUMBERS = ["1e3", "2e3"]
METRICS = ['wip', 'ip']


rule all:
    input:
        expand("data/kwip/{rn}-{metric}.dist", rn=READ_NUMBERS, metric=METRICS),
        expand("data/kwip/{rn}.stat", rn=READ_NUMBERS),
        "data/aligned_genomes.fasta",
        expand("data/samples/{genome}-{sample}-{rn}_il.fastq.gz", genome=GENOMES,
               sample=SAMPLES, rn=READ_NUMBERS)


rule root_genome:
    output:
        "data/genome.fa"
    params:
        length="1024",
        seed="23"
    log:
        "data/log/genomes/genome.log"
    shell:
        "mason_genome"
        " -l {params.length}"
        " -o {output}"
        " --seed {params.seed}"
        " >{log} 2>&1"


rule genomes:
    input:
        "data/genome.fa"
    output:
        fa="data/genomes/{genome}.fa",
        vcf="data/genomes/{genome}.vcf"
    params:
        seed=lambda w: GENOME_SEEDS[w.genome]
    log:
        "data/log/genomes/{genome}.log"
    shell:
        "mason_variator"
        " -ir {input}"
        " --snp-rate 0.005"
        " --small-indel-rate 0.0001"
        " -ov {output.vcf}"
        " -of {output.fa}"
        " --seed {params.seed}"
        " >{log} 2>&1"


rule alignment:
    input:
        expand("data/genomes/{genome}.fa", genome=GENOMES)
    output:
        "data/aligned_genomes.fasta"
    log:
        "data/log/alignment.log"
    shell:
        "cat {input}"
        " | mafft  --nuc -"
        " >{output}"
        " 2>{log}"


rule samples:
    input:
        "data/genomes/{genome}.fa",
    output:
        r1=temp("data/samples/{genome}-{sample}-{rn}_R1.fastq.gz"),
        r2=temp("data/samples/{genome}-{sample}-{rn}_R2.fastq.gz")
    params:
        seed=lambda w: SAMPLE_SEEDS[w.genome][w.sample],
        rn=lambda w: str(int(float(w.rn)))
    log:
        "data/log/samples/{genome}-{sample}-{rn}.log"
    shell:
        "mason_simulator"
        " -ir {input}"
        " --illumina-read-length 101"
        " -o {output.r1}"
        " -or {output.r2}"
        " --seed {params.seed}"
        " -n {params.rn}"
        " >{log} 2>&1"

rule ilfq:
    input:
        r1="data/samples/{genome}-{sample}-{rn}_R1.fastq.gz",
        r2="data/samples/{genome}-{sample}-{rn}_R2.fastq.gz"
    output:
        "data/samples/{genome}-{sample}-{rn}_il.fastq.gz"
    log:
        "data/log/join/{genome}-{sample}-{rn}.log"
    shell:
        "pairs join"
        " {input.r1}"
        " {input.r2}"
        " 2>>{log}"
        " | sickle se "
        " -f /dev/stdin"
        " -t sanger"
        " -o /dev/stdout"
        " -n"
        " 2>>{log}"
        " | gzip > {output}"
        " 2>>{log}"


rule hash:
    input:
        "data/samples/{genome}-{sample}-{rn}_il.fastq.gz"
    output:
        "data/hashes/{genome}-{sample}-{rn}.ct.gz"
    params:
        x='1e7',
        N='1',
        k='20'
    log:
        "data/log/hashes/{genome}-{sample}-{rn}.log"
    shell:
        "load-into-counting.py"
        " -N {params.N}"
        " -x {params.x}"
        " -k {params.k}"
        " -b"
        " -s tsv"
        " {output}"
        " {input}"
        " >{log} 2>&1"


rule kwip:
    input:
        expand("data/hashes/{genome}-{sample}-{{rn}}.ct.gz",
               genome=GENOMES, sample=SAMPLES)
    output:
        d="data/kwip/{rn}-{metric}.dist",
        k="data/kwip/{rn}-{metric}.kern"
    params:
        metric= lambda w: '-U' if w.metric == 'ip' else ''
    log:
        "data/log/kwip/{rn}-{metric}.log"
    threads:
        4
    shell:
        "kwip"
        " {params.metric}"
        " -d {output.d}"
        " -k {output.k}"
        " -t {threads}"
        " {input}"
        " >{log} 2>&1"


rule kwip_stats:
    input:
        expand("data/hashes/{genome}-{sample}-{{rn}}.ct.gz",
               genome=GENOMES, sample=SAMPLES)
    output:
        "data/kwip/{rn}.stat"
    log:
        "data/log/kwip-entvec/{rn}.log"
    threads:
        4
    shell:
        "kwip-entvec"
        " -o {output}"
        " -t {threads}"
        " {input}"
        " >{log} 2>&1"
