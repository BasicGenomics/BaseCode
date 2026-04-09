rule generate_polya:
    input: bam = "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"]),
           bai = "results/{name}.stitched.molecules.sorted.bam.bai".format(name=config["name"])
    output: bed = temp("results/downstream/{name}_polyA.bed".format(name=config["name"])), 
            merged_bed = temp("results/downstream/{name}_polyA_merged.bed".format(name=config["name"])),
            genome = temp("results/downstream/{name}_genome_polyA.txt".format(name=config["name"])),
            bedgraph = "results/downstream/{name}_polyA.bedgraph".format(name=config["name"]),
            bed_gz = "results/downstream/{name}_polyA.bed.gz".format(name=config["name"]), 
            merged_bed_gz = "results/downstream/{name}_polyA_merged.bed.gz".format(name=config["name"]),
            done = "results/dones/{name}_generate_polyA.done".format(name=config["name"])
    params: index = REFINDEX
    shell: """
        python workflow/scripts/generate_sites_bed.py --input {input.bam} --output results/downstream/{config[name]} --site-type polya
        bedtools merge -i {output.bed} -c 4 -o count > {output.merged_bed}
        cut -f1,2 {params.index} | sort -k1,1V > {output.genome}
        bedtools genomecov -i {output.bed} -g {output.genome} -bg > {output.bedgraph} 
        bgzip -c {output.bed} > {output.bed_gz}
        bgzip -c {output.merged_bed} > {output.merged_bed_gz}
        tabix -p bed {output.bed_gz}
        tabix -p bed {output.merged_bed_gz}
        touch {output.done}
        """

rule generate_tss:
    input: bam = "results/{name}.stitched.molecules.sorted.bam".format(name=config["name"]),
           bai = "results/{name}.stitched.molecules.sorted.bam.bai".format(name=config["name"])
    output: bed = temp("results/downstream/{name}_TSS.bed".format(name=config["name"])), 
            merged_bed = temp("results/downstream/{name}_TSS_merged.bed".format(name=config["name"])),
            genome = temp("results/downstream/{name}_TSS_genome.txt".format(name=config["name"])),
            bedgraph = "results/downstream/{name}_TSS.bedgraph".format(name=config["name"]),
            bed_gz = "results/downstream/{name}_TSS.bed.gz".format(name=config["name"]), 
            merged_bed_gz = "results/downstream/{name}_TSS_merged.bed.gz".format(name=config["name"]),
            done = "results/dones/{name}_generate_TSS.done".format(name=config["name"])
    params: index = REFINDEX
    shell: """
        python workflow/scripts/generate_sites_bed.py --input {input.bam} --output results/downstream/{config[name]} --site-type tss
        bedtools merge -i {output.bed} -c 4 -o count > {output.merged_bed}
        cut -f1,2 {params.index} | sort -k1,1V > {output.genome}
        bedtools genomecov -i {output.bed} -g {output.genome} -bg > {output.bedgraph} 
        bgzip -c {output.bed} > {output.bed_gz}
        bgzip -c {output.merged_bed} > {output.merged_bed_gz}
        tabix -p bed {output.bed_gz}
        tabix -p bed {output.merged_bed_gz}
        touch {output.done}
        """