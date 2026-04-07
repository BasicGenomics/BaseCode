rule generate_report:
    input: samplesheet = "results/metadata/{name}_samplesheet.csv".format(name=config["name"]),
           summary = "results/QC_files/{name}_summary_stats.csv".format(name=config["name"]),
           long_form= "results/QC_files/{name}_long_form_reconstruction_stats.csv".format(name=config["name"])
    output: "results/{name}_run_report.pdf".format(name=config["name"])
    shell: """
    python3 workflow/scripts/run_report.py --run_id {config[name]} --samplesheet {input.samplesheet} --summary {input.summary} --long-form-summary {input.long_form} --outdir results \
       --background-pdf workflow/resources/Background_v1.1.png  --readtype-schema workflow/resources/Read_Type_Schema.png --font-regular workflow/resources/MonaSans-Regular.ttf --font-bold workflow/resources/MonaSans-Bold.ttf
    """
