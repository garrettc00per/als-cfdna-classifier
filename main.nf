#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bam_dir = "/home/ec2-user/celfie/bam"
params.outdir = "results"
params.metadata = "/home/ec2-user/celfie/celfie_cfDNA_ss.csv"

// Create channel from BAM files
Channel
    .fromPath("${params.bam_dir}/*.bam")
    .map { file -> 
        def sample_id = file.baseName.replaceAll(/\.deduplicated\.sorted_ds10mill_chr21/, '')
        tuple(sample_id, file)
    }
    .set { bam_ch }

process EXTRACT_INSERT_SIZES {
    tag "$sample_id"
    publishDir "${params.outdir}/insert_sizes", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_insert_sizes.txt")
    
    script:
    """
    samtools view -f 0x2 ${bam} | \\
        awk '{if (\$9 > 0) print \$9}' > ${sample_id}_insert_sizes.txt
    """
}

process EXTRACT_END_MOTIFS {
    tag "$sample_id"
    publishDir "${params.outdir}/end_motifs", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_end_motifs.txt")
    
    script:
    """
    ${projectDir}/bin/extract_end_motifs.sh ${bam} ${sample_id}
    """
}

process EXTRACT_METHYLATION {
    tag "$sample_id"
    publishDir "${params.outdir}/methylation", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    path("${sample_id}_methylation.txt")
    
    script:
    """
    ${projectDir}/bin/extract_methylation.sh ${bam} ${sample_id}
    """
}

process EXTRACT_POSITIONS {
    tag "$sample_id"
    publishDir "${params.outdir}/positions", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple path("${sample_id}_position_stats.txt"), path("${sample_id}_position_bins.txt")
    
    script:
    """
    ${projectDir}/bin/extract_positions.sh ${bam} ${sample_id}
    """
}

process COMPUTE_INSERT_STATS {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(insert_sizes)
    
    output:
    path("${sample_id}_stats.txt")
    
    script:
    """
    ${projectDir}/bin/compute_insert_stats.sh ${insert_sizes} ${sample_id} > ${sample_id}_stats.txt
    """
}

process COMBINE_INSERT_STATS {
    
    input:
    path(stats_files)
    
    output:
    path("insert_size_summary.txt")
    
    script:
    """
    echo "sample_id mean median stddev min max q25 q75 n_fragments" > insert_size_summary.txt
    cat ${stats_files} >> insert_size_summary.txt
    """
}

process COMBINE_METHYLATION {
    
    input:
    path(meth_files)
    
    output:
    path("methylation_summary.txt")
    
    script:
    """
    echo "sample_id cpg_meth_rate chg_meth_rate chh_meth_rate overall_meth_rate total_cpg total_chg total_chh n_reads" > methylation_summary.txt
    cat ${meth_files} >> methylation_summary.txt
    """
}

process COMBINE_POSITIONS {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(pos_stats)
    
    output:
    path("position_summary.txt")
    
    script:
    """
    echo "sample_id region1_count region2_count region3_count total_fragments" > position_summary.txt
    cat ${pos_stats} >> position_summary.txt
    """
}

process CREATE_MOTIF_FEATURES {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(motif_files)
    
    output:
    path("motif_features.csv")
    
    script:
    """
    mkdir -p motif_dir
    cp ${motif_files} motif_dir/
    ${projectDir}/bin/create_motif_features.py motif_dir motif_features.csv
    """
}

process COMBINE_ALL_FEATURES {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(insert_summary)
    path(meth_summary)
    path(motif_features)
    path(metadata)
    
    output:
    path("all_features_combined.csv")
    
    script:
    """
    ${projectDir}/bin/combine_all_features.py ${insert_summary} ${meth_summary} ${motif_features} ${metadata} all_features_combined.csv
    """
}

process CLASSIFY {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(features)
    
    output:
    path("classification_results.csv")
    path("classification_output.txt")
    
    script:
    """
    ${projectDir}/bin/classify_combined.py ${features} classification_results.csv > classification_output.txt
    """
}

process VISUALIZE {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path(features)
    
    output:
    path("visualizations_*.png")
    path("visualizations_feature_importances.csv")
    
    script:
    """
    ${projectDir}/bin/create_visualizations.py ${features} visualizations
    """
}

workflow {
    // Extract all features from BAM files
    insert_sizes_ch = EXTRACT_INSERT_SIZES(bam_ch)
    end_motifs_ch = EXTRACT_END_MOTIFS(bam_ch)
    methylation_ch = EXTRACT_METHYLATION(bam_ch)
    positions_ch = EXTRACT_POSITIONS(bam_ch)
    
    // Compute insert size stats
    stats_ch = COMPUTE_INSERT_STATS(insert_sizes_ch)
    insert_summary_ch = COMBINE_INSERT_STATS(stats_ch.collect())
    
    // Combine methylation stats
    meth_summary_ch = COMBINE_METHYLATION(methylation_ch.collect())
    
    // Combine position stats
    pos_summary_ch = COMBINE_POSITIONS(positions_ch.map { it[0] }.collect())
    
    // Create motif features
    motif_features_ch = CREATE_MOTIF_FEATURES(
        end_motifs_ch.map { it[1] }.collect()
    )
    
    // Load metadata
    metadata_ch = Channel.fromPath(params.metadata)
    
    // Combine ALL features with disease labels (now including positions!)
    all_features_ch = COMBINE_ALL_FEATURES(
        insert_summary_ch,
        meth_summary_ch,
        motif_features_ch,
        metadata_ch
    )
    
    // Run classification
    classification_ch = CLASSIFY(all_features_ch)
    
    // Create visualizations
    VISUALIZE(all_features_ch)
}
