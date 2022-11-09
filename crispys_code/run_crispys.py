import Stage0

Stage0.CRISPys_main("/groups/itay_mayrose/udiland/crispys_off_target/test_crunch/fam_rice_omer/HOM05D001900/HOM05D001900.fa",
                    "/groups/itay_mayrose/udiland/crispys_test/test_files_git/for_debug/out1",
                    alg='gene_homology', singletons=0, max_target_polymorphic_sites=12, internal_node_candidates=200,
                    omega=0.45, set_cover=False, pams=0,
                    get_n_candidates=2)



# 3 genes out1 "/groups/itay_mayrose/udiland/crispys_off_target/test_crunch/fam_rice_omer/HOM05D001900/HOM05D001900.fa"
# 2 genes out2 /groups/itay_mayrose/udiland/crispys_off_target/test_crunch/fam_rice_omer/HOM05D001890/HOM05D001890.fa
#
# out 8 genes /groups/itay_mayrose/udiland/crispys_test/test_files_server/HOM04D000632/HOM04D000632.txt
# alg='default'