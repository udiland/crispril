import server_run

server_run.run_server(
    fasta_file="/groups/itay_mayrose/udiland/crispys_off_target/test_crunch/fam_rice_omer/HOM05D001899/HOM05D001899.fa",
    output_path="/groups/itay_mayrose/udiland/crispys_off_target/test_crunch/fam_rice_omer/HOM05D001899",
    alg='default',
    where_in_gene=0.8,
    use_thr=1,
    omega=0.6,
    scoring_function="cfd_funct",
    min_length=20,
    max_length=20,
    start_with_g=False,
    internal_node_candidates=200,
    max_target_polymorphic_sites=12,
    pams=0,
    off_targets='Not_selected')
