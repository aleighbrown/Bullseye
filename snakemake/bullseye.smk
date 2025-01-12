configfile: "config.yaml"
include: "helpers.py"

import sys
import os
import glob
import re
# Function to get files that match a specific pattern
def get_files(pattern,SAMPLE_NAMES):
    regex = re.compile(pattern)
    files = [f for f in SAMPLE_NAMES if regex.match(f)]
    return files


# Function to collect all files for a given contrast
def get_contrast_files(wildcards):
    files = []
    contrast = wildcards.contrast
    for control in CONTROL_NAMES:
        print(control)
        pattern = out_dir_FindEdit + f"{contrast}.{control}.bed"
        files.append(pattern)
    print(files)
    return files


in_bam_dir = config["input_dir"]
bam_suffix = config["bam_suffix"]
out_dir_parseBAM = config["out_dir_parseBAM"]
path_to_perl_code = config["code_dir"]
out_dir_FindEdit = config["out_dir_FindEdit"]

contrast_pattern = "full"
control_pattern = "control_truncated"

SAMPLE_NAMES = [f.replace(bam_suffix, "") for f in os.listdir(in_bam_dir) if f.endswith(bam_suffix)]
CONTRAST_NAMES = get_files(f".*{contrast_pattern}.*",SAMPLE_NAMES)
CONTROL_NAMES = get_files(f".*{control_pattern}.*",SAMPLE_NAMES)




wildcard_constraints:
    sample = "|".join(SAMPLE_NAMES),
    contrast = "|".join(CONTRAST_NAMES),
    baseline = "|".join(CONTROL_NAMES)

# Adjust the all rule to include the output of parseBAM
rule all:
    input:
        expand(os.path.join(out_dir_parseBAM, "{sample}.parse.output.gz"), sample=SAMPLE_NAMES),
        expand(os.path.join(out_dir_parseBAM, "{sample}.FINISHED"), sample=SAMPLE_NAMES),
        expand(os.path.join(out_dir_FindEdit, "{contrast}.merged.bed"), 
                contrast=[os.path.splitext(f)[0] for f in get_files(f".*{contrast_pattern}.*",SAMPLE_NAMES)]),
        expand(os.path.join(out_dir_FindEdit, "{contrast}.merged.RAC.bed"), 
                contrast=[os.path.splitext(f)[0] for f in get_files(f".*{contrast_pattern}.*",SAMPLE_NAMES)]),
        expand(os.path.join(out_dir_FindEdit,"{contrast}.{baseline}.bed"), 
               contrast=[os.path.splitext(f)[0] for f in get_files(f".*{contrast_pattern}.*",SAMPLE_NAMES)],
               baseline=[os.path.splitext(f)[0] for f in get_files( f".*{control_pattern}.*",SAMPLE_NAMES)])



# Rule for parsing BAM files
rule parseBAM:
    input:
        bam = os.path.join(in_bam_dir, "{sample}" + bam_suffix),
    output:
        bam_out = os.path.join(out_dir_parseBAM, "{sample}.parse.output.gz"),
        test_out = os.path.join(out_dir_parseBAM, "{sample}.FINISHED")
    params:
        parseBamParams = return_parsed_extra_params(config['ParseBam_params']),
        cpu = 4,
        perl_script = os.path.join(path_to_perl_code, "parseBAM.pl"),
        output_base = lambda wildcards: os.path.join(out_dir_parseBAM, wildcards.sample + ".parse.output")
    shell:
        """
        set +u;
        source activate Bullseye
        perl {params.perl_script} \
        --input {input.bam} \
        --output {params.output_base} \
        --cpu {params.cpu} \
        {params.parseBamParams}
        touch {output.test_out}
        """

rule process_FindEdit:
    input:
        contrast=out_dir_parseBAM + "{contrast}.parse.output.gz",
        baseline=out_dir_parseBAM + "{baseline}.parse.output.gz"
    output:
        out=out_dir_FindEdit + "{contrast}.{baseline}.bed"
    params:
        editSiteParams = return_parsed_extra_params(config['FindEditSite_params']),
        cpu = 4,
        perl_script = os.path.join(path_to_perl_code, "Find_edit_site.pl"),
        ref_flat = config['refFlat']

    shell:
        """
        set +u;
        source activate Bullseye
        perl {params.perl_script} \
        --annotationFile {params.ref_flat} \
        --EditedMatrix {input.contrast} \
        --controlMatrix {input.baseline} \
        --cpu {params.cpu} \
        --outfile {output.out} \
        --verbose \
        {params.editSiteParams}
        """

# Rule to merge files for each contrast
rule merge_contrast_files:
    input:
        # THESE FILES ARE GENERATED BY THE RULE ABOVE process_FindEdit
        contrast_files = get_contrast_files
    output:
        merged_file = out_dir_FindEdit + "{contrast}.merged.bed"
    params:
        perl_script = os.path.join(path_to_perl_code, "summarize_sites.pl"),
        the_value = out_dir_FindEdit + "{contrast}.*.bed"
    shell:
        """
        set +u;
        source activate Bullseye
        perl {params.perl_script} --MinRep 3 \
        --mut 2 \
        --repOnly \
        {params.the_value} > {output.merged_file}
        """

# Rule to merge files for each contrast
rule rac_filter:
    input:
        # This input is not directly used in the shell command, but ensures the rule waits for the creation of these files.
        merged_file = out_dir_FindEdit + "{contrast}.merged.bed"
    output:
        merged_file = out_dir_FindEdit + "{contrast}.merged.RAC.bed"
    params:
        bash_script = os.path.join(path_to_perl_code, "RACfilter.sh"),
        out_dir = out_dir_FindEdit,
        pattern = "{contrast}.*.bed"
    shell:
        """
        set +u;
        cd {params.out_dir}
        source activate Bullseye
        bash {params.bash_script} -r {input.merged_file} 
        """
