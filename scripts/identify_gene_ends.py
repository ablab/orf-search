#!/usr/bin/env python3
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from joblib import Parallel, delayed

import sys
import argparse
import logging

import os
from os import listdir
from os.path import isfile, isdir, join

import multiprocessing

import load_mappings
from assembly_graph import Graph
from gene_ends_finder import GeneEndsFinder
from ends_path_constructor import PathConstructor, EndsPathConstructor, translate_str

def overlap(p1, p2, g):
    e1 = set(p1)
    e2 = set(p2)
    sm = 0
    for e in e1 & e2:
        sm += len(g.edges[e])
    return sm

def supported_by_contig(p1, p2, g, overlap_len):
    p1_str = ",".join(p1)
    p2_str = ",".join(p2)
    if p1_str in p2_str or p2_str in p1_str:
        return True
    has_good_overlap = False
    for i1 in range(len(p1)):
        for i2 in range(len(p2)):
            p1_str = ",".join(p1[:i1 + 1])
            p2_str = ",".join(p2[i2:])
            if p1_str == p2_str:
                if overlap(p1[:i1 + 1], p2[i2:], g) >= overlap_len:
                    has_good_overlap = True
                    break
            p1_str = ",".join(p1[i1:])
            p2_str = ",".join(p2[:i2 + 1])
            if p1_str == p2_str:
                if overlap(p1[i1:], p2[:i2 + 1], g) >= overlap_len:
                    has_good_overlap = True
                    break
        if has_good_overlap:
            break
    return has_good_overlap

def compare_with_contig_paths(name, paths, g, uniqueedge_len):
    res_path = []
    for p in paths:
        contigs = set()
        for e in p["Edges"]:
            if e in g.edge_paths:
                for c in g.edge_paths[e]:
                    contigs.add(c)
        supported_edges = set()
        supported = True
        for c in contigs:
            overlap_len = overlap(p["Edges"], g.paths[c], g)
            has_good_overlap = supported_by_contig(p["Edges"], g.paths[c], g, overlap_len)
            if not has_good_overlap:
                for e in g.paths[c]:
                    if e in p["Edges"] and len(g.edges[e]) > uniqueedge_len:
                        if len(g.graph[e].keys()) > 0 and len(g.graph[g.revert_edge(list(g.graph[e].keys())[0])]) > 1 and len(g.graph[g.revert_edge(e)]) > 1:
                            supported = False
                            break
            if not supported:
                break
        if supported:
            res_path.append(p)
    return res_path

def is_inside(aln1, aln2, path_constructor):
    path1_s = ",".join(aln1["path"])
    path2_s = ",".join(aln2["path"])
    if path1_s in path2_s and \
      (path1_s != path2_s or \
      (path1_s == path2_s and (aln1["start"] > aln2["start"] and aln1["end"] < aln2["end"]))):
        index = 0
        for i in range(len(aln2["path"])):
            j = 0
            while j < len(aln1["path"]) and i + j < len(aln2["path"]) and aln2["path"][i + j] == aln1["path"][j]:
                j += 1
            if j == len(aln1["path"]):
                index = i
                break
        if i == 0 and aln2["start"] > aln1["start"] - 1:
            return False
        ln = path_constructor.restore_path_len(aln2["start"], aln1["start"] - 1, aln2["path"][:i + 1])
        if ln % 3 == 0:
            return True
        else:
            return False
    else:
        return False

def remove_covered_orfs(alns, startcodon_dists, g):
    res_alns1 = []
    res_startcodon_dists1 = []
    path_constructor = PathConstructor(g)
    for i in range(len(alns)):
        is_covered = False
        path_str = path_constructor.restore_path(alns[i]["start"], alns[i]["end"], alns[i]["path"])
        if len(path_str) % 3 == 0:
            t_s = translate_str(path_str)
            if t_s.find("*") == len(t_s) - 1 or t_s.find("*") == -1:
                res_alns1.append(alns[i])
                res_startcodon_dists1.append(startcodon_dists[i])
    res_alns = []
    res_startcodon_dists = []
    for i in range(len(res_alns1)):
        is_covered = False
        for j in range(len(res_alns1)):
            if is_inside(res_alns1[i], res_alns1[j], path_constructor):
                is_covered = True
                break
        if not is_covered:
            res_alns.append(res_alns1[i])
            res_startcodon_dists.append(res_startcodon_dists1[i])
    logging.info( u'Reduced number of anchors from ' + str(len(alns)) + ' to ' + str(len(res_alns)))
    return res_alns, res_startcodon_dists

def generate_orf(args):
    aln, g, startcodon_dist, only_longest, uniqueedge = args[0], args[1], args[2], args[3], args[4]
    logging.debug(aln["name"])
    genes_finder = GeneEndsFinder(g)
    start_codon_pos, stop_codon_pos = genes_finder.find_ends_positions(aln["start"], aln["end"], aln["path"], startcodon_dist, only_longest)
    logging.debug(u'Start codon num=' + str(len(start_codon_pos)) + ' Stop codons num=' + str(len(stop_codon_pos)))
    ends_path_contructor = EndsPathConstructor(g)
    all_paths = ends_path_contructor.restore_full_paths(aln, start_codon_pos, stop_codon_pos, startcodon_dist)
    all_paths = compare_with_contig_paths(aln["name"], all_paths, g, uniqueedge)
    logging.debug(u'Paths num=' + str(len(all_paths)))
    return {"name": aln["name"], "all_paths": all_paths}

def generate_orfs(output, output_shortest, alns, g, startcodon_dists, only_longest, uniqueedge, t):
    logging.debug( u'Threads ' + str(t) + u' alns ' + str(len(alns)))
    alns, startcodon_dists = remove_covered_orfs(alns, startcodon_dists, g)
    all_orfs = Parallel(n_jobs=t, require='sharedmem')(delayed(generate_orf)([alns[i], g, startcodon_dists[i], only_longest, uniqueedge]) for i in range(len(alns)))
    with open(output, "a+") as fout:    
        for orf in all_orfs:
            name = orf["name"]
            for path in orf["all_paths"]:
                fout.write(">" + path["prefix"] + "_" + name + "|Edges=" + "_".join(path["Edges"]) + "|" + \
                           "|".join([k + "=" + str(path[k]) for k in path.keys() if k not in {"Edges", "seq", "prefix"}]) + "\n")
                fout.write(path["seq"] + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print potential ORFs for all given graph alignments in fasta-format')
    parser.add_argument('-s', '--sequences', help='fasta-file with SPAligner alignments', required=False)
    parser.add_argument('-m', '--hmms', help='path to PathRacer results folder', required=False)
    parser.add_argument('-g', '--graph', help='path to assembly graph (in gfa-format)', required=True)
    parser.add_argument('-k', '--kmer', help='k-mer size in graph', required=True)
    parser.add_argument('-p', '--proteins',  help='list of genes to estimate hmms position on genes (domtbl has to be set)', required=False)
    parser.add_argument('-d', '--domtbl',  help='HMMer alignment of hmms to genes (genes has to be set)', required=False)
    parser.add_argument('-e', '--evalue',  help='minimum e-value for HMM alignment', default=0.000000001)
    parser.add_argument('-l', '--minlen',  help='minimum length', default=0.9)
    parser.add_argument('-f', '--longestorf', help='generate ORFs that have stop codon before start codon', action='store_true')
    parser.add_argument('-u', '--uniqueedge',  help='length of unique edges in nucs', default=500, type=int)
    parser.add_argument('-o', '--out',  help='output prefix', required=True)
    parser.add_argument('-t', '--threads', help='threads number', required=False)
    args = parser.parse_args()
    logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level = logging.DEBUG, filename = args.out + u'.log')
    if args.hmms == None and args.sequences == None:
        logging.error( u'Please provide sequences or hmm alignments')
        exit(-1)
    K = int(args.kmer)
    g = Graph(args.graph)
    output = args.out + ".fasta"
    output_shortest = args.out + "_shortest.fasta"
    if os.path.exists(output):
        os.remove(output)
    if os.path.exists(output_shortest):
        os.remove(output_shortest)

    alns = []
    startcodon_dists = []
    if args.hmms != None:
        if (args.proteins == None or args.domtbl == None):
            logging.warning( u'Information about HMMs positions is not provided')
            hmm_hits = {}
        else:
            hmm_hits = load_mappings.load_true_crygenes(args.proteins, args.domtbl, {}, K)

        filenames = [f for f in listdir(args.hmms) \
                     if isfile(join(args.hmms, f)) and \
                        isfile(join(args.hmms, f[:-9] + "edges.fa")) and \
                        f.endswith("domtblout") ]
        for f in filenames:
            f_path = join(args.hmms, f)
            alns.extend(load_mappings.load_pathracer_mapping(f_path, g.edges, float(args.minlen), float(args.evalue), K))
        for aln in alns:
            startcodon_dist = []
            if aln["name"] in hmm_hits:
                startcodon_dist = [x["start"] for x in hmm_hits[aln["name"]] ]
            startcodon_dists.append(startcodon_dist)
        logging.info( u'HMMs: ' + str(len(alns)))

    if args.sequences != None:
        seq_alns = load_mappings.load_spaligner_mapping(args.sequences)
        alns.extend(seq_alns)
        for a in seq_alns:
            startcodon_dists.append([a["d"]])
        logging.info( u'Seqs: ' + str(len(alns)))

    generate_orfs(output, output_shortest, alns, g, startcodon_dists, args.longestorf, args.uniqueedge, int(args.threads))



