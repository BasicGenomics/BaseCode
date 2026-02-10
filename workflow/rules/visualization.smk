rule generate_report:
    input: samplesheet = "results/metadata/{name}_samplesheet.csv".format(name=config["name"]), summary = "results/QC_files/{name}_summary_stats.csv".format(name=config["name"]), long_form_summary = "results/QC_files/{name}_long_form_reconstruction_stats.csv".format(name=config["name"])
    output: report = "results/{name}_run_report.pdf".format(name=config["name"])
    shell: "echo Generating Run Report &&  python3 workflow/scripts/run_report.py --samplesheet {input.samplesheet} --summary {input.summary} --long-form-summary {input.long_form_summary} --outdir results"