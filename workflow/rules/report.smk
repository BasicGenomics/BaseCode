rule generate_report:
    input: samplesheet = "results/metadata/{name}_samplesheet.csv".format(name=config["name"]),
           summary_stats = "results/QC_files/{name}_summary_stats.csv".format(name=config["name"]),
           long_form_reconstruction_stats = "results/QC_files/{name}_long_form_reconstruction_stats.csv".format(name=config["name"])
    output: run_report = "results/{name}_run_report.pdf".format(name=config["name"]),
            done = "results/dones/{name}_generate_report.done".format(name=config["name"])
    log: "results/logs/generate_report.log"        
    shell: """
    python3 workflow/scripts/run_report.py --run-id {config[name]} --samplesheet {input.samplesheet} --summary {input.summary_stats} --long-form-summary {input.long_form_reconstruction_stats} --outdir results --background-pdf workflow/resources/Background_v1.1.png  --readtype-schema workflow/resources/Read_Type_Schema.png --font-regular workflow/resources/MonaSans-Regular.ttf --font-bold workflow/resources/MonaSans-Bold.ttf > {log} 2>&1
    touch {output.done}   
    """
