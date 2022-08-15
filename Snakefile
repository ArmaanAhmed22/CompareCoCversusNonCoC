rule convert_seer_output_to_csvs:
    input:
        "Input/{file}/seer.csv"
    output:
        directory("Input/{file}/data")
    script:
        "Pipeline/convert_data_to_csvs.py"

rule compare:
    input:
        directory("Input/{file}/data")
    output:
        directory("Output/{file}/kmfs"),
        directory("Input/{file}/kmf_univariate_models")
    script:
        "Pipeline/compare.py"
    
rule generate_diff:
    input:
        directory("Input/{file}/kmf_univariate_models")
    output:
        directory("Output/{file}/diffs")
    script:
        "Pipeline/generate_diff.py"