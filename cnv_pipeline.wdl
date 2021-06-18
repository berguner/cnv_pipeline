version 1.0


workflow cnv_pipeline {
    input{
        String project_name
        String project_folder
        String pipeline_folder
        String assembly
        File bed_file
        File sample_annotation_sheet

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
            --sample_names ~{sample_names}
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
        Int cpus = 16
        Int memory = 32000
        String partition = "covid"
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

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }
}