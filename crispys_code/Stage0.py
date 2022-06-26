__author__ = 'Gal Hyams'
import sys
import CasSites
import Stage1
import Stage1_h
import Distance_matrix_and_UPGMA #from UPGMA.py, more informative name
import timeit
import pickle
import Metric
import argparse
import os
import random
import make_tree_display_CSV #from server


random.seed(1234)

# get the path of this script file
PATH = os.path.dirname(os.path.realpath(__file__))

def sort_expectation(candidates_DS, homology):
    def sort_subgroup(candidates_DS):
        candidates_DS.sort(key = lambda item: (item.cut_expectation, item.total_num_of_mm()), reverse=True)
    if not homology:
        sort_subgroup(candidates_DS)
    else:
        for i in range(len(candidates_DS)):
            sort_subgroup(candidates_DS[i].candidate_lst)

def sort_thr(candidates_DS, Omega, homology):
    '''sort the candidates DS by num of genes with cut prob> Omega and then by the probobility to cleave all of these genes'''
    def sort_subgroup(candidates_DS, Omega):
        for candidate in candidates_DS:
            num_of_genes_above_thr = 0
            cleave_all = 1
            for gene, score in candidate.genes_score_dict.items():
                if score >= Omega:
                    cleave_all *= score
                    num_of_genes_above_thr += 1
            candidate.cleve_all_above_thr = cleave_all
            candidate.num_of_genes_above_thr = num_of_genes_above_thr
        candidates_DS.sort(key = lambda item: (item.num_of_genes_above_thr, item.cleave_all_above_thr), reverse = True)
    if not homology:
        sort_subgroup(candidates_DS, Omega)
    else:
        for i in range(len(candidates_DS)):
            sort_subgroup(candidates_DS[i].candidate_lst, Omega)

def leave_only_relevant_sgRNA(res):
    if len(res) < 1:
        return
    candidates_to_del = []
    for i in range(len(res) -1, -1,-1):
        if res[i].cut_expectation < 1:
            del res[i]
        elif i < len(res) - 1:
            for j in range(i+1,len(res)):
                if j >= len(res):
                    continue
                if res[i].seq == res[j].seq: # there is no need in both of them. Who is sutable for more genes?
                    if res[i].cut_expectation <= res[j].cut_expectation:
                        del res[i]
                    else:
                        del res[j]


def CRISPys_main(fasta_file, path, alg = 'A', where_in_gene = 1, use_thr = 0, Omega = 1, df_targets = Metric.cfd_funct, protdist_outfile = "outfile", min_length= 20, max_length = 20,start_with_G = False, internal_node_candidates = 10, PS_number = 12, PAMs=0):
    start = timeit.default_timer()
    random.seed(1234) #Omer 16/06/22
    cfd_dict = None
    if isinstance(where_in_gene, str):
        where_in_gene = float(where_in_gene.strip())
    if isinstance(Omega, str):
        Omega = float(Omega.strip())
    if isinstance(use_thr, str):
        use_thr = int(use_thr.strip())
    #choosing the distance function
    if df_targets == "gold_off" or df_targets == "goldoff":
       df_targets = Distance_matrix_and_UPGMA.gold_off_func
    if df_targets == "MITScore" or df_targets == "CrisprMIT":
        df_targets = Distance_matrix_and_UPGMA.MITScore
    if df_targets == "cfd_funct" or df_targets == "cfd_func" or df_targets == "cfd"\
            or df_targets == Metric.cfd_funct:
        df_targets = Metric.cfd_funct
        cfd_dict = pickle.load(open(PATH + "/cfd_dict.p",'rb'))
        
    if df_targets == "CCTop" or df_targets == "ccTop" :
        df_targets = Distance_matrix_and_UPGMA.ccTop
    # add an option to different pam (taken from server version by Udi)
    if PAMs == 0:
        PAMs = ['GG']
    elif PAMs == 1:
        PAMs = ['GG', 'AG']
    protdist_outfile = path + "/" + protdist_outfile

    #print(df_targets)
    original_range_in_gene = [0, where_in_gene]
    genes_sg_dict = {}
    sg_genes_dict = {}
    genesNames = []
    genesList = []
    f = open(fasta_file,'r')
    gene_name = ""
    gene_seq = ""
    lines = f.readlines()
    i = 0
    genes_exons_dict = {}  #key: gene name. value: list of exons
    while i <= len(lines):
    #stage 1: make  gene: sequence dictionary
        if i == len(lines) or lines[i][0] == '>':
            if len(gene_seq) > 0 and gene_name != "": #add the gene
                if gene_name not in genes_exons_dict:
                    genes_exons_dict[gene_name] = [gene_seq]
                else:
                    genes_exons_dict[gene_name] = genes_exons_dict[gene_name] + [gene_seq]
                gene_seq = ""
            if i != len(lines): # lines[i-1][0] == '>':
                gene_name = lines[i][1:].strip() #without the '>' and the '\n'
        elif lines[i] != "\n":
            gene_seq += lines[i].strip()
        i+=1
    #stage 2: find the target sites
    for gene_name in genes_exons_dict.keys():

        genes_sg_dict[gene_name] = CasSites.get_targets_sites_from_exons_lst(genes_exons_dict[gene_name],df_targets, original_range_in_gene, min_length, max_length, start_with_G, PAMs)
        genesNames.append(gene_name)
        genesList.append("".join(genes_exons_dict[gene_name]))
        #filling up the sg_genes_dict
        for sg in genes_sg_dict[gene_name]:
            if sg in sg_genes_dict:
                sg_genes_dict[sg] = sg_genes_dict[sg] + [gene_name]
            else:
                sg_genes_dict[sg] = [gene_name]
    if alg == 'E':

        res = Stage1_h.call_it_all(genesList, genesNames, sg_genes_dict, genes_sg_dict, Omega, protdist_outfile, path, df_targets, internal_node_candidates, cfd_dict, PS_number)

    else:
        res = Stage1.call_it_all(genesList, genesNames, sg_genes_dict, genes_sg_dict, Omega, protdist_outfile, path, df_targets, cfd_dict, PS_number) #thies line have been change to be sutable for wrapper
    if use_thr:
        sort_thr(res, Omega, alg == 'E')
    else:
        sort_expectation(res, alg == 'E')

    #remove the folowing two lines when using CRISPysCover
    # if len(res)>200: #commented by Udi 03032022
    #     res = res[:200]

    # old output function. commented by Udi 13/04/2022
    # Stage1.print_res_to_csvV2(res, sg_genes_dict, genesList, genesNames, path, alg == 'E')
    # Stage1.print_res_to_csvV3(res, sg_genes_dict, genesList, genesNames, path, alg =='E')

    pickle.dump(res, open(path + "/res_in_lst.p", "wb"))
    pickle.dump(genesNames, open(path + "/genesNames.p", "wb"))
    # add saving the geneList in pickle in order to produce the results like in the server version - Udi 28/02/22
    pickle.dump(genesList, open(path + '/genesList.p', 'wb'))
    pickle.dump(sg_genes_dict, open(path + "/sg_genes.p", "wb"))

    # new output function taken from the crispys server code. Udi 13/04/2022
    make_tree_display_CSV.tree_display(path, alg == 'E')

    # # make a removed repetition results. taken from server
    # removed_rep = remove_repetitions_in_targets_sites(res, alg, use_thr, Omega)
    # pickle.dump( removed_rep, open(path + '/res_in_lst_removed_rep.p', 'wb'))

    stop = timeit.default_timer()
    "time: ", stop - start
    with open(f"{path}/time.txt", 'w') as f: # Omer 15/06/22
        f.write(str(stop - start))
    return res




def parse_arguments(parser):
    #def CRISPys_main(fasta_file, path , alg = 'A', where_in_gene = 1, use_thr = 0,  Omega = 1, df_targets = Metric.cfd_funct, protdist_outfile = "outfile", min_length= 20, max_length = 20,start_with_G = False, internal_node_candidates = 10, PS_number = 12):
    parser.add_argument('fasta_file', type=str, metavar='<fasta_file>', help='The path to the input fasta file')
    parser.add_argument('path', type=str, metavar='<path>', help='THe path to the directory in which the output files will be written')
    parser.add_argument('--alg', type=str, default='A', help='Choose E for considering homology')
    parser.add_argument('--where_in_gene', type=float, default=1, help='input a number between 0 to 1 in order to ignore targets sites downstream to the corresponding gene prefix')
    parser.add_argument('--t', type=bool, default=0, help='for using sgRNA to gain maximal gaining score among all of the input genes or 1 for the maximal cleavage likelihood only among genes with score higher than the average. Default: 0.')
    parser.add_argument('--v', type=float, default=0.43, help='the value of the threshold. A number between 0 to 1 (included). Default: 0.43')
    parser.add_argument('--s', type=str, default='cfd_funct', help='the scoring function of the targets. Optional scoring systems are: cfd_funct (default), gold_off, CrisprMIT and CCtop. Additinal scoring function may be added by the user or by request.')
    parser.add_argument('--p', type=str, default='outfile', help='protDist output file name. Default: "outfile"')
    parser.add_argument('--l', type=int, default=20, help='minimal length of the target site. Default:20')
    parser.add_argument('--m', type=bool, default=20, help = 'maximal length of the target site, Default:20')
    parser.add_argument('--g', type=bool, default=0, help='1 if the target sites are obligated to start with a G codon or 0 otherwise. Default: 0.')
    parser.add_argument('--i', type=int, default=10, help='when choosing the consider homology option, this is the number of sgRNAs designed for each homology sub-group. Default: 10')
    parser.add_argument('--ps', type=int, default=12, help='the maximal number of possible polymorphic sites in a target. Default: 12')
    parser.add_argument('--PAMs', type=int, default=0, help='0 to search NGG pam or 1 to search for NGG and NAG. Default: 0')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    CRISPys_main(fasta_file = args.fasta_file,
                 path = args.path,
                 alg = args.alg,
                 where_in_gene = args.where_in_gene,
                 use_thr = args.t,
                 Omega=args.v,
                 df_targets = args.s,
                 protdist_outfile= args.p,
                 min_length=args.l,
                 max_length=args.m,
                 start_with_G = args.g,
                 internal_node_candidates=args.i,
                 PS_number=args.ps,
                 PAMs=args.PAMs)
    #CRISPys_main(*sys.argv[1:])

