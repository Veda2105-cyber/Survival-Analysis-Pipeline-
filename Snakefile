
configfile: "config/config.yaml"
print("CONFIG LOADED:", config)
rule all:
    input:
        config["checked_output"],
        config["cleaned_output"],
        config["merged_data_output"],
        config["normalized_expression"],
        config["filtered_mutations_output"],
        config["top_mutated_genes"],
        config["merged_omics_output"],
        config["expression_heatmap"],
        config["methylation_heatmap"],
        config["cna_heatmap"],
        config["labeled_survival_output"],
        config["final_model_data_output"],
        config["top_cox_features"],
        config["tuned_cox_model"]        
        
rule check_input_files:
    output:
        config["checked_output"]
    input:
        config["input_files"]
    shell:
        """
        mkdir -p {config[results_dir]}
        echo "All input files found:" > {output}
        for file in {input}; do
            echo $file >> {output}
        done
        """
rule clean_data:
    input:
        config["input_files"]
    output:
        config["cleaned_output"]
    script:
        "scripts/clean_data.py"

rule merge_data:
    input:
        clinical_sample=config["input_files"][1],  
        mutations=config["input_files"][4]         
    output:
        config["merged_data_output"]
    script:
        "scripts/merge_data.py"

rule normalize_expression:
    input:
        "data/data_mrna_illumina_microarray.txt"
    output:
        "results/normalized_expression.txt"
    script:
        "scripts/normalize_expression.py"

rule filter_mutations:
    input: 
        "data/data_mutations.txt"
    output: 
        "results/filtered_mutations.txt"
    params:
        top_n=config["top_n_mutated_genes"]
    script: 
        "scripts/filter_mutations.py"

rule get_top_mutated_genes:
    input:
        config["filtered_mutations_output"]
    output:
        "results/top_mutated_genes.txt"
    params:
        n=config["top_n_mutated_genes"]
    script:
        "scripts/get_top_mutated_genes.py"


rule merge_top_gene_multiomics:
    input:
        top_genes=config["top_mutated_genes"],
        expression=config["normalized_expression"],
        methylation=config["methylation_file"],
        cna=config["cna_file"]
    output:
        config["merged_omics_output"]
    script:
        "scripts/merge_top_gene_multiomics.py"


rule visualize_expression_heatmap:
    input: 
        config["merged_omics_output"]
    output: 
	config["expression_heatmap"]
    params: 
	suffix="_expr", title="Expression Heatmap"
    conda: 
	"envs/heatmap.yaml"
    script: 
	"scripts/plot_heatmap.py"


rule visualize_methylation_heatmap:
    input: 
	config["merged_omics_output"]
    output: 
	config["methylation_heatmap"]
    params: 
	suffix="_meth", title="Methylation Heatmap"
    conda: 
	"envs/heatmap.yaml"
    script: 
	"scripts/plot_heatmap.py"


rule visualize_cna_heatmap:
    input: 
	config["merged_omics_output"]
    output: 
	config["cna_heatmap"]
    params: 
	suffix="_x", title="CNA Heatmap"
    conda: 
	"envs/heatmap.yaml"
    script: 
	"scripts/plot_heatmap.py"


rule label_survival_outcomes:
    input:
        config["clinical_patient_file"]
    output:
        config["labeled_survival_output"]
    conda:
        "envs/heatmap.yaml"
    script:
        "scripts/label_survival.py"


rule create_model_ready_matrix:
    input:
        omics = config["merged_omics_output"],
        labels = config["labeled_survival_output"]
    output:
        config["final_model_data_output"]
    conda:
        "envs/heatmap.yaml"
    script:
        "scripts/merge_with_labels.py"


rule select_top_cox_features:
    input:
        data=config["final_model_data_output"]
    output:
        selected=config["top_cox_features"]
    script:
        "scripts/select_top_cox_features.py"


rule tune_cox_model:
    input:
        data=config["final_model_data_output"]
    output:
        config["tuned_cox_model"]
    script:
        "scripts/tune_cox_model.py"
