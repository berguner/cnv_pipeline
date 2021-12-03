version 1.0

workflow cnv_pipeline {
    input{
        String project_name
        String project_folder
        String pipeline_folder
        String assembly
        File bed_file
        File sample_annotation_sheet
        File gtf_file

        File all_gene_symbol_file
        File canonical_file
        File pli_score_file
        File hi_score_file
    }
    call count_reads {
        input:
            project_folder = project_folder,
            pipeline_folder = pipeline_folder,
            bed_file = bed_file
    }

    call cluster_samples {
        input:
            project_folder = project_folder,
            pipeline_folder = pipeline_folder,
            read_counts_list = count_reads.read_counts_list
    }

    Array[Array[String]] clusters = read_tsv(cluster_samples.sample_clusters)

    scatter(cluster in clusters) {
        call codex {
            input:
                project_folder = project_folder,
                pipeline_folder = pipeline_folder,
                gtf_file = gtf_file,
                bed_file = bed_file,
                assembly = assembly,
                cluster_label = cluster[0],
                sample_names = cluster[1]
        }

        call exomedepth {
            input:
                project_folder = project_folder,
                pipeline_folder = pipeline_folder,
                bed_file = bed_file,
                assembly = assembly,
                cluster_label = cluster[0],
                sample_names = cluster[1]
        }
    }

    call extract_codex_sample_results {
        input:
            project_folder = project_folder,
            pipeline_folder = pipeline_folder,
            codex_qc = codex.codex_qc,
            assembly = assembly
    }

    call annotate_cnv {
        input:
            project_folder = project_folder,
            pipeline_folder = pipeline_folder,
            assembly = assembly,
            codex_log = extract_codex_sample_results.log_file,
            exomedepth_logs = exomedepth.exomedepth_log
    }

    call aggregate_cnv_deletions {
        input:
            project_name = project_name,
            project_folder = project_folder,
            pipeline_folder = pipeline_folder,
            sample_annotation_sheet = sample_annotation_sheet,
            gtf_file = gtf_file,
            all_gene_symbol_file = all_gene_symbol_file,
            canonical_file = canonical_file,
            pli_score_file = pli_score_file,
            hi_score_file = hi_score_file,
            codex_log = extract_codex_sample_results.log_file,
            exomedepth_logs = exomedepth.exomedepth_log
    }

    output {
        File read_counts_list = count_reads.read_counts_list
        File sample_clusters = cluster_samples.sample_clusters
    }
}

task count_reads {
    input {
        String project_folder
        String pipeline_folder
        File bed_file
        String bam_folder
        String bam_prefix
        String bam_suffix

        # runtime parameters
        Int cpus = 8
        Int memory = 16000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    String output_folder = "~{project_folder}/read_counts"
    String count_script= "~{pipeline_folder}/codex_exome_coverage.R"
    command <<<
        [ ! -d "~{output_folder}" ] && mkdir -p ~{output_folder};

        for i in ~{bam_folder}/*~{bam_suffix}; do \
            BN=$(basename $i);
            SN=${BN/~{bam_prefix}/};
            SN=${SN/~{bam_suffix}/};
            [ ! -f "~{output_folder}/${SN}_exome_RC.csv.gz" ] && echo -e "Rscript --vanilla  ~{count_script} --bed ~{bed_file} --project_folder ~{project_folder} --bam ~{bam_folder}/${BN} --sample_name ${SN} > ~{output_folder}/${SN}_exome_RC.csv.gz.log 2>&1";
        done | parallel -j ~{cpus} :::

        for i in ~{output_folder}/*exome_RC.csv.gz; do \
            BN=$(basename $i);
            SN=${BN/_exome_RC.csv.gz/};
            echo -e "${SN}\t${i}"; done > "~{output_folder}/read_counts_list.tsv"

    >>>

    output {
        File read_counts_list = "~{output_folder}/read_counts_list.tsv"
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }

}

task cluster_samples {
    input {
        String project_folder
        String pipeline_folder
        File read_counts_list
        Int min_mean_rc = 20
        Int min_cluster_size = 12

        # runtime parameters
        Int cpus = 2
        Int memory = 8000
        String partition = "tinyq"
        String time = "2:00:00"
        String? rt_additional_parameters
        String? rt_image
    }
    String cluster_script = "~{pipeline_folder}/cluster_samples.py"
    command <<<
        cat ~{read_counts_list};
        python3 ~{cluster_script} \
            -p ~{project_folder} \
            --minimum_mean_read_count ~{min_mean_rc} \
            --min_cluster_size ~{min_cluster_size} > ~{project_folder}/sample_clusters.log 2>&1;
    >>>

    output {
        File sample_clusters = "~{project_folder}/sample_clusters.tsv"
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }
}

task exomedepth {
    input {
        String project_folder
        String pipeline_folder
        File bed_file
        String cluster_label
        String sample_names
        String assembly

        # runtime parameters
        Int cpus = 2
        Int memory = 8000
        String partition = "shortq"
        String time = "12:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    String output_folder = "~{project_folder}/exomedepth_results"
    String exomedepth_script = "~{pipeline_folder}/exomedepth_CNV.R"
    command <<<
        [ ! -d "~{output_folder}" ] && mkdir -p ~{output_folder};
        Rscript --vanilla ~{exomedepth_script} \
            --bed ~{bed_file} \
            --cluster ~{cluster_label} \
            --project_folder ~{project_folder} \
            --sample_names ~{sample_names} \
            --force n \
            --assembly ~{assembly} > ~{output_folder}/exomedepth_cls~{cluster_label}.log 2>&1
    >>>

    output {
        File exomedepth_log = "~{output_folder}/exomedepth_cls~{cluster_label}.log"
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }
}

task codex {
    input {
        String project_folder
        String pipeline_folder
        File bed_file
        String cluster_label
        String sample_names
        String assembly
        File gtf_file

        # runtime parameters
        Int cpus = 8
        Int memory = 32000
        String partition = "shortq"
        String time = "12:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    String output_folder = "~{project_folder}/codex_results"
    String codex_script = "~{pipeline_folder}/codex_chromosome_CNV.R"
    command <<<
        [ ! -d "~{output_folder}" ] && mkdir -p ~{output_folder};

        chromosomes=''
        for i in {1..22}; do if [ "~{assembly}" = "b37" ] ; then chromosomes="$chromosomes $i"; else chromosomes="$chromosomes chr$i"; fi done

        for chromosome in ${chromosomes}; do \
            echo -e "Rscript --vanilla ~{codex_script} --chromosome ${chromosome} --bed ~{bed_file} --gtf ~{gtf_file} --cluster ~{cluster_label} --project_folder ~{project_folder} --sample_names ~{sample_names} --force y --assembly ~{assembly} --write_sample_output n > ~{output_folder}/codex_cls~{cluster_label}_chr${chromosome}.log 2>&1";
        done | parallel -j ~{cpus} :::
    >>>

    output {
        File codex_qc = "~{output_folder}/codex_qc_cls_~{cluster_label}.tsv"
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }
}

task extract_codex_sample_results {
    input {
        String project_folder
        String pipeline_folder
        Array[File] codex_qc
        String assembly

        # runtime parameters
        Int cpus = 1
        Int memory = 1000
        String partition = "tinyq"
        String time = "2:00:00"
        String? rt_additional_parameters
        String? rt_image
    }
    String extraction_script = "~{pipeline_folder}/extract_codex_sample_results.py"
    command <<<
        CHROMOSOME_PREFIX=''
        if [ "~{assembly}" = "b37" ] ; then CHROMOSOME_PREFIX=""; else CHROMOSOME_PREFIX="chr"; fi
        python3 ~{extraction_script} \
            -c ${CHROMOSOME_PREFIX} \
            -p ~{project_folder}  > ~{project_folder}/codex_results/extract_codex_sample_results.log 2>&1;
    >>>

    output {
        File log_file = "~{project_folder}/codex_results/extract_codex_sample_results.log"
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }
}

task annotate_cnv {
    input {
        String project_folder
        String pipeline_folder
        File codex_log
        String assembly
        Array[File] exomedepth_logs

        # runtime parameters
        Int cpus = 1
        Int memory = 4000
        String partition = "shortq"
        String time = "12:00:00"
        String? rt_additional_parameters
        String? rt_image
    }
    String annotation_script = "~{pipeline_folder}/annotate_cnv.sh"
    command <<<
        [ ! -d "~{project_folder}/annotated_results" ] && mkdir ~{project_folder}/annotated_results ;
        ASSEMBLY=''
        if [ "~{assembly}" = "b37" ] ; then ASSEMBLY="GRCh37"; else ASSEMBLY="GRCh38"; fi

        for i in ~{project_folder}/read_counts/*RC.csv.gz;
        do
            BN=$(basename $i);
            SN=${BN/_exome_RC.csv.gz/};
            [ ! -f "~{project_folder}/annotated_results/${SN}_AnnotSV.tsv" ] && bash ~{annotation_script} ~{project_folder} ${SN} ${ASSEMBLY} >> ~{project_folder}/annotated_results/AnnotSV.log 2>&1;
        done
    >>>

    output {
        File log_file = "~{project_folder}/annotated_results/AnnotSV.log"
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }
}

task aggregate_cnv_deletions {
    input {
        String project_name
        String project_folder
        String pipeline_folder
        String sample_annotation_sheet
        File gtf_file
        File all_gene_symbol_file
        File canonical_file
        File pli_score_file
        File hi_score_file
        File codex_log
        Array[File] exomedepth_logs

        # runtime parameters
        Int cpus = 1
        Int memory = 8000
        String partition = "shortq"
        String time = "12:00:00"
        String? rt_additional_parameters
        String? rt_image
    }
    String aggregation_script = "~{pipeline_folder}/aggregate_cnv_deletions.R"
    command <<<
        [ ! -d "~{project_folder}/aggregated" ] && mkdir ~{project_folder}/aggregated ;

        Rscript --vanilla ~{aggregation_script} \
            --project_id ~{project_name} \
            --sample_sheet ~{sample_annotation_sheet} \
            --project_folder ~{project_folder} \
            --gtf ~{gtf_file} \
            --all_gene_symbol_file ~{all_gene_symbol_file} \
            --canonical_file ~{canonical_file} \
            --pli_score_file ~{pli_score_file} \
            --hi_score_file ~{hi_score_file} > ~{project_folder}/aggregated/aggregate_cnv_deletions.log 2>&1
    >>>

    output {
        File log_file = "~{project_folder}/aggregated/aggregate_cnv_deletions.log"
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }
}