#!/usr/bin/env python3
import sys
import os
from os.path import isfile, isdir, join
import tempfile

import logging

import subprocess
import resource

import scripts.load_mappings
import scripts.input_sanity_check

import argparse
import yaml

execution_path = os.path.dirname(os.path.abspath(__file__))

def align_hmms(hmms_file, graph_file, k, evalue, maxsize, threads, out_dir):
    com = execution_path + "/aligners/pathracer " + hmms_file + " " + graph_file + " " + str(k) \
        + " --output " + out_dir + " --rescore --annotate-graph --threads " + str(threads) \
        + " -E " + evalue + " --domE " + evalue + " --max-size " + maxsize + " > " + out_dir + ".log"
    logging.info( u'Running HMM alignment. See log in ' + out_dir + u'.log')
    logging.debug( u'Running: ' + com)
    return_code = subprocess.call([com], shell=True)
    return out_dir, return_code

def align_sequences(graph_file, k, protein_file, threads, out_prefix):
    com = execution_path + "/aligners/spaligner " + execution_path + "/aligners/spaligner_cfg.yaml -g " + graph_file + \
                    " -k " + str(k) + " -s " + protein_file + " -t " + str(threads) + \
                     " -d protein -o " + out_prefix + " > " + out_prefix + ".log"
    logging.info( u'Running Sequence Alignment. See log in ' + out_prefix + u'.log')
    logging.debug( u'Running: ' + com)
    return_code = subprocess.call([com], shell=True)
    return out_prefix + "/alignment.fasta", return_code

def find_true_hmm_alignments(hmmer_path, proteins_file, hmms_file, evalue, threads, out_file):
    com = hmmer_path + "hmmsearch --domtblout " + out_file + \
                    ".dtbl -E " + evalue + " --cpu " + str(threads)  + " " + hmms_file + " " + proteins_file \
                    + " > " + out_file + "_true.log"
    logging.info( u'Running Start Codon Distance Estimator. See log in ' + out_file + u'_true.log')
    logging.debug( u'Running: ' + com)
    return_code = subprocess.call([com], shell=True)
    return out_file + ".dtbl", return_code

def extract_ORFs_from_graph(hmms_alignments, proteins_alignments, graph_file, k, proteins_file, \
                            hmms_true_alignments, threads, out_dir):
    com = execution_path + "/scripts/identify_gene_ends.py "
    if os.path.exists(proteins_alignments):
        com += "-s " + proteins_alignments
    if os.path.exists(hmms_alignments):
        com += " -m " + hmms_alignments
    if proteins_file != None:
        com += " -p " + proteins_file \
               +" -d " + hmms_true_alignments
    com += " -g " + graph_file + " -k " + str(k) + " -t " + str(threads) +" -o " + out_dir
    logging.info( u'Extracting ORFs from assembly graph. See log in ' + out_dir + u'/orfs_raw.log')
    logging.debug( u'Running: ' + com)
    return_code = subprocess.call([com], shell=True)
    return out_dir + "/orfs_raw.fasta", return_code

def filter_orfs(orfs_sequences, graph, proteins_file, contigs_file, threads, print_all, out_dir):

    com = execution_path + "/scripts/filter_ORFs.py -s " + orfs_sequences + \
                        " -g " + graph + " -t " + str(threads) + " -o " + out_dir
    if not print_all:
        if contigs_file != None:
            com += " -c " + contigs_file
        if proteins_file != None:
            com += " -p " + proteins_file
    logging.info( u'Filtering ORFs. See log in ' + out_dir + u'/orfs_final.log')
    logging.debug( u'Running: ' + com)
    return_code = subprocess.call([com], shell=True)
    return out_dir + "/orfs_final.fasta", return_code

def form_yaml(args, outputdir):
    config_name = os.path.join(outputdir, "config.yaml")
    config = {}
    config["pathracer"] = {}
    config["pathracer"]["evalue"] = args.pr_evalue
    config["pathracer"]["min_length"] = args.pr_minlen
    config["pathracer"]["max_size"] = args.pr_maxsize

    config["orfs_search"] = {}
    config["orfs_search"]["longest"] = args.os_longest
    config["orfs_search"]["max_restorable_length"] = args.os_maxrestorelen
    config["orfs_search"]["max_path_length"] = args.os_maxpathlen
    config["orfs_search"]["min_path_length"] = args.os_minpathlen
    config["orfs_search"]["max_paths_num"] = args.os_maxpathnum
    config["orfs_search"]["min_unique_edge_length"] = args.os_minuniqueedgelen

    config["orfs_filtering"] = {}
    config["orfs_filtering"]["min_reliable_length"] = args.of_minreliablelen
    config["orfs_filtering"]["reliable_length"] = args.of_reliablelen
    config["orfs_filtering"]["identity"] = args.of_identity
    config["orfs_filtering"]["ed_threshold"] = args.of_ed
    config["orfs_filtering"]["all"] = args.of_all

    config_file = open(config_name, "w")
    yaml.dump(config, config_file)
    config_file.close()
    return config_name

def identify_graph_ksize(gfafile):
    k = -1
    with open(gfafile, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith("L"):
                k = int(ln.strip().split("\t")[-1][:-1])
    if k == -1:
        logging.error( u'K-mer size cannot be identified from gfa-file')
        exit(-1)
    return k

def main(args):
    parser = argparse.ArgumentParser(description='Searches for potential genes in assembly graphs')
    parser.add_argument('-m', '--hmms', help='file with hmms in HMMer format', required=True)
    parser.add_argument('-s', '--sequences', help='fasta-file with proteins(optional)', required=False)
    parser.add_argument('-r', '--runspaligner', help='if set, run SPAligner on proteins', action='store_true')
    parser.add_argument('-g', '--graph', help='assembly graph in gfa-format', required=True)
    parser.add_argument('-c', '--contigs', help='fasta-file with assembly contigs', required=False)
    parser.add_argument('-t', '--threads', help='number of threads', required=False)
    parser.add_argument('-o', '--out', help='output directory', required=True)
    parser.add_argument('--pr-evalue', help='PathRacer: e-value (DEFAULT: 1.0e-09)', default=0.000000001, type=float, required=False)
    parser.add_argument('--pr-minlen', help='PathRacer: minimum length of alignment (DEFAULT: 0.9)', default=0.9, type=float, required=False)
    parser.add_argument('--pr-maxsize', help='PathRacer: maximum component size to process (DEFAULT: 2500000)', default=2500000, type=int, required=False)
    parser.add_argument('--os-longest', help='ORF search: generate ORFs that have stop codon before start codon (DEFAULT: False)', action='store_false')
    parser.add_argument('--os-maxrestorelen', help='ORF search: maximum distance in the graph from alignment where start/stodon can be found (DEFAULT: 3000)', default=3000, type=int, required=False)
    parser.add_argument('--os-maxpathlen', help='ORF search: maximum length of path to start/stop codon (DEFAULT: 3000)', default=1500, type=int, required=False)
    parser.add_argument('--os-minpathlen', help='ORF search: minimum length of path to start/stop codon (DEFAULT: 0)', default=0, type=int, required=False)
    parser.add_argument('--os-maxpathnum', help='ORF search: maximum number of paths can be generated during reconstruction of paths to start/stop codons (DEFAULT: 2000)', default=2000, type=int, required=False)
    parser.add_argument('--os-minuniqueedgelen', help='ORF search: length of edges that considered us unique during filtering based on contigs (DEFAULT: 0)', default=0, type=int, required=False)
    parser.add_argument('--of-minreliablelen', help='ORF filtering: second minimum length of edge for filtering (used only if there is no edges with length > reliable length) (DEFAULT: 100)', default=100, type=int, required=False)
    parser.add_argument('--of-reliablelen', help='ORF filtering: minimum length of edge for filtering (DEFAULT: 3000)', default=300, type=int, required=False)
    parser.add_argument('--of-identity', help='ORF filtering: single-linkage clustering cluster two sequence into the same cluster if they are similar > indetity (DEFAULT: 90)', default=90, type=int, required=False)
    parser.add_argument('--of-ed', help='ORF filtering: threshold for edit distance calculation (DEFAULT: 0.2)', default=0.2, type=float, required=False)
    parser.add_argument('--of-all', help='ORF filtering: do not perform filtering based on contigs or known IPGs (DEFAULT: False)', required=False, action='store_false')

    is_test = False
    if len(args) == 2 and args[1] == "--test":
        is_test = True
        test_dir = "tiny_dataset_test"
        p = os.path.abspath(__file__)
        args = [args[0], "-m", p[:-len("orf_search.py")]  + "/tiny_dataset/ricinb_lectin2.hmm", "-s", p[:-len("orf_search.py")] + "/tiny_dataset/toxin.fasta", \
                         "-r", "-g", p[:-len("orf_search.py")]  + "/tiny_dataset/graph.gfa", "-o", test_dir]
    args = parser.parse_args(args[1:])

    if not os.path.exists(args.out):
        os.makedirs(args.out)

    config_file = form_yaml(args, args.out)

    logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level = logging.DEBUG)
    #logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level = logging.DEBUG, filename = args.out + u'/orf_search.log')

    if is_test:
        logging.info(u'Start test on a small dataset...')

    hmmer_path = ""
    evalue, maxsize = args.pr_evalue, args.pr_maxsize
    kmer = identify_graph_ksize(args.graph)

    t = 1
    if args.threads != None:
        t = int(args.threads)

    if args.sequences != None:
        check_return_code, sequences_filename = scripts.input_sanity_check.check_sequences(args.sequences, args.out)
        if not check_return_code:
            logging.warning( u'Proteins file contains non-unique sequences, equal sequences were collapsed')
    else:
        sequences_filename = None

    hmms_name = ".".join(args.hmms.split("/")[-1].split(".")[0:-1])
    hmm_return_str, hmm_return_code = align_hmms(args.hmms, args.graph, kmer, str(evalue), str(maxsize), min(t, 16), join(args.out, hmms_name))
    if hmm_return_code != 0:
        logging.warning( u'HMM alignment failed')
    if hmm_return_code == 0 and sequences_filename != None and not os.path.exists(join(args.out, hmms_name + ".dtbl")):
        find_true_hmm_alignments(hmmer_path, sequences_filename, args.hmms, str(evalue), t, join(args.out, hmms_name))

    if sequences_filename != None and args.runspaligner:
        seq_name = "seq_aln"
        seq_return_str, seq_return_code = align_sequences(args.graph, kmer, sequences_filename, t, join(args.out, seq_name))
        if seq_return_code != 0:
            logging.warning( u'Sequence alignment failed')
    else:
        seq_return_str = ""
        seq_return_code = 0

    if seq_return_code != 0 and hmm_return_code != 0:
        logging.error( u'No data was generated to find orfs')
        exit(-1)

    orfs_fasta, return_code = extract_ORFs_from_graph(hmm_return_str, seq_return_str, args.graph, kmer, \
                                                        sequences_filename, join(args.out, hmms_name + ".dtbl"), \
                                                        t, args.out)

    if return_code != 0:
        logging.error( u'Orfs generation failed')
        exit(-1)

    final_orfs_str, return_code = filter_orfs(orfs_fasta, args.graph, sequences_filename, args.contigs, t, args.of_all, args.out)
    if return_code != 0:
        logging.error( u'Filtering failed')
        exit(-1)

    logging.info( u'ORFs search finished. Please find results in: ' + args.out + "/")
    logging.info( u'Potential novel ORFs can be found in: ' + args.out + "/orfs_final_most_reliable.fasta")

    if is_test:
        logging.info(u'The test finished successfully!')



if __name__ == "__main__":
    main(sys.argv)
