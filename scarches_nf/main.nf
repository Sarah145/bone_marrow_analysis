#!/usr/bin/env nextflow

ch_input = Channel.fromPath("data/*.h5ad").map{ it -> [it.baseName, it]}

// stage, split parameters for each
hvg_stage = params.hvg? params.hvg.split(',').collect{it.trim()} : []
knn_stage = params.knn? params.knn.split(',').collect{it.trim()} : []
latent_stage = params.latent ? params.latent.split(',').collect{it.trim()} : []

// place in channel
hvg_ch = Channel.value(hvg_stage).flatten()
knn_ch = Channel.value(knn_stage).flatten()
latent_ch = Channel.value(latent_stage).flatten()

// Cartesian product
scvi_params = hvg_ch.combine(latent_ch)

process RUN_SCVI{
        publishDir "./data/latent"
        tag "scVI $database $hvg hvgs $latent latent dims"
        
        input:
        tuple val(base), file(database), val(hvg), val(latent) from ch_input.combine(scvi_params)

        output:
        set val(base), val(hvg), val(latent), file('*_latent.h5ad') into UMAP

        script:
        """
        run_scvi.py\
                --input_object ${database}\
		--output_prefix ${base}_${hvg}_${latent}\
		--batch_key ${params.batch_key}\
                --hvg ${hvg}\
                --latent ${latent}
        """
}

process umap{
        publishDir "./data/umap"
        tag "umap $datain"
        
        input:
        set val(knn), val(base), val(hvg), val(latent), file(datain) from knn_ch.combine(UMAP)

        output:
        set val(base), val(hvg), val(latent), val(knn), file('*.csv') into EVAL

        script:
        """
        run_umap.py\
                --input_object ${datain}\
                --output_prefix ${base}_${hvg}_${latent}_${knn}\
                --knn $knn\
		--batch_key ${params.batch_key}\
		--celltype_key ${params.celltype_key}
        """
}

process eval{
        publishDir "./data/eval"
        tag "eval $datain"
        
        input:
        set val(base), val(hvg), val(latent), val(knn), file(datain) from EVAL

        output:
        stdout to out
        file('*_eval.csv')

        script:
        """
        eval.R\
                --input_object ${datain}\
                --output_prefix ${base}_${hvg}_${latent}_${knn}\
		--batch_key ${params.batch_key}\
		--celltype_key ${params.celltype_key}
        """
}

