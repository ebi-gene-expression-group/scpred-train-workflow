#!/usr/bin/env nextflow 

// extract training data
TRAINING_DATA = Channel.fromPath(params.training_10x_dir)
TRAINING_METADATA = Channel.fromPath(params.metadata_file)

// if necessary, down-sample cells to avoid memory issues 
process downsample_cells {
    conda "${baseDir}/envs/label_analysis.yaml"

    memory { 32.GB * task.attempt }
    maxRetries 5
    errorStrategy { task.attempt<=5 ? 'retry' : 'ignore' }

    input:
        file(expression_data) from TRAINING_DATA
        file(training_metadata) from TRAINING_METADATA
        
    output:
        file("expr_data_downsampled") into TRAINING_DATA_DOWNSAMPLED
        file("metadata_filtered.tsv") into TRAINING_METADATA_DOWNSAMPLED

    """
    set +e
    downsample_cells.R\
        --expression-data ${expression_data}\
        --metadata ${training_metadata}\
        --exclusions ${params.exclusions}\
        --cell-id-field ${params.cell_id_col_name}\
        --cell-type-field ${params.cell_types_col_name}\
        --output-dir expr_data_downsampled\
        --metadata-upd metadata_filtered.tsv

    if [ \$? -eq  2 ];
    then
        cp -P ${expression_data} expr_data_downsampled
        cp -P ${training_metadata} metadata_filtered.tsv
        exit 0
    fi
    """
}

// parse data into Seurat object 
process read_training_data{
  conda "${baseDir}/envs/seurat.yaml"

  memory { 32.GB * task.attempt }
  maxRetries 5
  errorStrategy { task.attempt<=5 ? 'retry' : 'ignore' }

  input: 
    file(training_10X_data) from TRAINING_DATA_DOWNSAMPLED
    file(training_metadata) from  TRAINING_METADATA_DOWNSAMPLED

  output:
    file("training_seurat.rds") into TRAINING_SEURAT

  """
  seurat-read-10x.R\
            --data-dir ${training_10X_data}\
            --output-format 'seurat'\
            --metadata ${training_metadata}\
            --output-object-file training_seurat.rds\
  """
}

process get_features{
  conda "${baseDir}/envs/scpred.yaml"

  memory { 32.GB * task.attempt }
  maxRetries 5
  errorStrategy { task.attempt<=5 ? 'retry' : 'ignore' }

  input:
    file(scpred_training_object) from TRAINING_SEURAT

  output:
    file("scpred_training_features.rds") into TRAINING_FEATURES

  """
  scpred_get_feature_space.R\
          --input-object ${scpred_training_object}\
          --prediction-column ${params.cell_types_col_name}\
          --output-path scpred_training_features.rds
  """
}

process train_model{
  publishDir "${params.results_dir}", mode: 'copy'
  conda "${baseDir}/envs/scpred.yaml"

  memory { 32.GB * task.attempt }
  maxRetries 5
  errorStrategy { task.attempt<=5 ? 'retry' : 'ignore' }

  input:
    file(scpred_training_features) from TRAINING_FEATURES

  output:
    file("scpred_classifier.rds") into TRAINED_MODEL

  """
  scpred_train_model.R\
          --input-object ${scpred_training_features}\
          --train-id ${params.training_dataset_id}\
          --model ${params.model}\
          --allow-parallel ${params.allow_parallel}\
          --num-cores ${params.num_cores}\
          --get-scpred ${params.get_scpred}\
          --output-path scpred_classifier.rds\
          --train-probs-plot ${params.train_probs_plot_path}
  """
}

