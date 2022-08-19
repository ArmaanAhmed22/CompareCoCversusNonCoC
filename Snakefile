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
        directory("Input/{file}/kmf_univariate_models"),
        directory("Input/{file}/lifelines")
    script:
        "Pipeline/compare.py"
    
rule generate_diff:
    input:
        directory("Input/{file}/kmf_univariate_models")
    output:
        directory("Output/{file}/diffs")
    script:
        "Pipeline/generate_diff.py"

rule division_diff:
    input:
        directory("Input/{file}/kmf_univariate_models"),
        directory("Input/{file}/lifelines")
    output:
        "Output/{file}/division_diff.png",
        "Input/{file}/p_values_difference_end.csv"
    script:
        "Pipeline/generate_division_diff.py"