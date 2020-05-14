#!/usr/bin/env nextflow 

// extract training data
TRAINING_DATA = Channel.fromPath(params.training_10x_dir)
TRAINING_METADATA = Channel.fromPath(params.metadata_file)

process read_training_data{
  conda "${baseDir}/envs/dropletutils.yaml"

  errorStrategy { task.exitStatus == 130 || task.exitStatus == 137  ? 'retry' : 'finish' }   
  maxRetries 10
  memory { 16.GB * task.attempt }

  maxForks 1

  input: 
    file(training_10X_data) from TRAINING_DATA
    file(training_metadata) from  TRAINING_METADATA

  output:
    file("training_sce.rds") into TRAINING_SCE

  """
  dropletutils-read-10x-counts.R\
            --samples ${training_10X_data}\
            --col-names ${params.col_names}\
            --metadata-files ${training_metadata}\
            --cell-id-column ${params.cell_id_col_name}\
            --output-object-file training_sce.rds
  """
}

process process_training_sce{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  maxForks 1

  input:
    file(training_sce) from TRAINING_SCE

  output:
    file("training_matrix.rds") into TRAINING_MATRIX
    file("training_labels.txt") into TRAINING_LABELS

  """
  scpred_preprocess_data.R\
            --input-sce-object ${training_sce}\
            --normalised-counts-slot ${params.normalised_counts_slot}\
            --output-matrix-object training_matrix.rds\
            --output-labels training_labels.txt 
  """
}

process eigen_decompose{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  maxForks 1

  input:
    file(training_matrix) from TRAINING_MATRIX
    file(training_labels) from TRAINING_LABELS

  output:
    file("scpred_training_object.rds") into TRAINING_OBJECT

  """
  scpred_eigen_decomp.R\
            --training-matrix ${training_matrix}\
            --log-transform ${params.log_transform}\
            --training-labels ${training_labels}\
            --output-path scpred_training_object.rds 
  """
}

process get_features{
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  maxForks 1

  input:
    file(scpred_training_object) from TRAINING_OBJECT

  output:
    file("scpred_training_features.rds") into TRAINING_FEATURES

  """
  scpred_get_feature_space.R\
          --input-object ${scpred_training_object}\
          --prediction-column ${params.cell_types_col_name}\
          --output-path scpred_training_features.rds\
          --eigenvalue-plot-path ${params.eigenvalue_plot_path}
  """
}

process train_model{
  publishDir "${params.results_dir}", mode: 'copy'
  conda "${baseDir}/envs/scpred.yaml"

  errorStrategy { task.attempt < 5  ? 'retry' : 'finish' }   
  maxRetries 5
  memory { 16.GB * task.attempt }

  maxForks 1

  input:
    file(scpred_training_features) from TRAINING_FEATURES

  output:
    file("scpred_trained_model.rds") into TRAINED_MODEL

  """
  scpred_train_model.R\
          --input-object ${scpred_training_features}\
          --train-id ${params.training_dataset_id}\
          --model ${params.model}\
          --output-path ${params.trained_model}\
          --train-probs-plot ${params.train_probs_plot_path}
  """
}

