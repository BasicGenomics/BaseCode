def get_final_output(wildcards):
    final_output = ["results/{name}.stitched.molecules.sorted.bam".format(name=config["name"])]

    return final_output