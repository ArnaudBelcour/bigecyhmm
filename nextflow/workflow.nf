params.str = 'EsMeCaTa + bigecyhmm + visualisation'

params.infile = "$launchDir/buchnera_workflow.tsv"
params.inAbundfile = "$launchDir/buchnera_workflow_abund.tsv"
params.precomputedDB = "$launchDir/esmecata_database.zip"
params.visualisationScript = "$launchDir/create_bigecyhmm_plot.py"

params.outputFolder = "$launchDir/output_folder"

input_file_path = Channel.fromPath(params.infile)
input_abundance_file_path = Channel.fromPath(params.inAbundfile)

precomputedDB_path = Channel.fromPath(params.precomputedDB)

outputFolder_path = Channel.fromPath(params.outputFolder)

visualisation_script_path = Channel.fromPath(params.visualisationScript)

process esmecata {
    input:
        path input_esmecata
        path esmecata_precomputed_db
    output:
        path 'output_1_esmecata', emit: output_1_esmecata, type: "dir"

    publishDir "${params.outputFolder}", mode: 'copy'

    script:
    """
    esmecata precomputed -i ${input_esmecata} -d ${esmecata_precomputed_db} -o output_1_esmecata
    """
}

process bigecyhmm {
    input:
        path esmecata_output_folder

    output:
        path 'output_2_bigecyhmm', emit: output_2_bigecyhmm, type: "dir"

    publishDir "${params.outputFolder}", mode: 'copy'

    script:
    """
    bigecyhmm  -i ${esmecata_output_folder}/1_clustering/reference_proteins_consensus_fasta -o output_2_bigecyhmm
    """
}

process visualisation{
    input:
        path visualisation_script_path
        path input_abundance_file_path
        path esmecata_output_folder
        path bigecyhmm_output_folder

    output:
        path 'output_3_visualisation', emit: output_3_visualisation, type: "dir"

    publishDir "${params.outputFolder}", mode: 'copy'

    script:
    """
    python3 ${visualisation_script_path} ${input_abundance_file_path} ${esmecata_output_folder} ${bigecyhmm_output_folder} output_3_visualisation
    """
}

workflow {
    esmecata(input_file_path, precomputedDB_path)
    bigecyhmm(esmecata.out.output_1_esmecata)
    visualisation(visualisation_script_path, input_abundance_file_path, esmecata.out.output_1_esmecata, bigecyhmm.out.output_2_bigecyhmm)
}