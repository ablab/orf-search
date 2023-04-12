#!/usr/bin/env python3
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from joblib import Parallel, delayed

import subprocess
import logging
import yaml

import sys
import argparse

import os
from os import listdir
from os.path import isfile, isdir, join

import re
import edlib

MIN_RELIABLE_LENGTH = 100
RELIABLE_LENGTH = 300
IDENTITY = 90
ED_THRESHOLD = 0.2

def edist(lst):
    if len(str(lst[0])) == 0:
        return -1, ""
    if len(str(lst[1])) == 0:
        return -1, ""
    ed_er = int(ED_THRESHOLD*max(len(str(lst[0])), len(str(lst[1]))))
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path", k = ed_er)
    return result["editDistance"], result["cigar"]

def aai(ar):
    p1, p2 = str(ar[0]), str(ar[1])
    if p1.endswith("*"):
        p1 = p1[:-1]
    if p2.endswith("*"):
        p2 = p2[:-1]
    ed, cigar = edist([str(p1), str(p2)])
    if ed == -1:
        return 0
    matches = re.findall(r'\d+=', cigar)
    aai = 0.0
    for m in matches:
        aai += int(m[:-1])
    aai /= max(len(p1), len(p2))
    return aai*100

def make_rc(s):
    seq = Seq(s)
    return str(seq.reverse_complement())

def load_gfa_edges(gfa_filename):
    res = {}
    coverage = {}
    with open(gfa_filename, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith("S"):
                lst = ln.strip().split("\t")[1:]
                node_id, seq, kc = lst[0], lst[1], lst[-1]
                res[node_id + "+"] = seq
                res[node_id + "-"] = make_rc(seq)
                coverage[node_id + "+"] = float(kc[len("KC:i:"):])/len(seq)
                coverage[node_id + "-"] = float(kc[len("KC:i:"):])/len(seq)
    return res, coverage
    
def load_fasta(filename):
    record_lst = list(SeqIO.parse(filename, "fasta"))
    return record_lst

def save_fasta(filename, orfs):
    with open(filename + ".fasta", "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

def align_with_nucmer(orfs, orfs_filename, contigs_filename, outdir, threads):
    res = []
    com = "nucmer -t " + threads + " --sam-short " + outdir + "/orfs_final.sam " + contigs_filename + " " + orfs_filename
    logging.info( u'Running nucmer: ' + com)
    subprocess.call([com], shell=True)
    aligned_set = set()
    pattern = re.compile("^[0-9]+M$")
    with open(outfile + ".sam", "r") as fin:
        for ln in fin.readlines():
            if not ln.startswith("@"):
                lst = ln.split("\t")
                name = lst[0]
                cigar = lst[5]
                if pattern.match(cigar):
                    aligned_set.add(name)
    #os.remove(outfile+".sam")
    logging.info( u'Contained in contigs: ' + str(len(aligned_set)))
    for orf in orfs:
        if orf.name not in aligned_set:
            res.append(orf) 
    return res


def filter_in_contigs(orfs, contigs):
    res = []
    for orf in orfs:
        found = False
        for c in contigs:
            if orf.seq in c.seq:
                found = True
                break
        if not found:
            res.append(orf)
    return res

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)

def translate_orfs(orfs):
    res = []
    for orf in orfs:
        res.append(make_record(orf.seq.translate(), orf.id, orf.name))
    return res

def leave_unique(orfs):
    res = []
    seq_set = set()
    for orf in orfs:
        seq_set.add(str(orf.seq))
    for s in seq_set:
        sid = []
        name = []
        desc = []
        for orf in orfs:
            if s == str(orf.seq):
                name.append(orf.name)
                sid.append(orf.id)
                desc.append(orf.description)
        res.append(make_record(Seq(s), ";".join(name), ";".join(sid), ";".join(desc)))
    return res

def leave_unknown(orfs, known_proteins):
    res = []
    for orf in orfs:
        unknown = True
        for p in known_proteins:
            if p.seq == orf.seq or (orf.seq[-1] == "*" and p.seq == orf.seq[:-1]):
                unknown = False
                break
        if unknown:
            res.append(orf)
    return res

def get_metainfo(s):
    keys = {"Edges":[], "apriori_startd_prob": 0, "coverage": 0}
    ss = s.split(";")[0]
    for c in ss.split("|"):
        for k in keys:
            if c.startswith(k):
                if k == "Edges":
                    keys[k] = c[len(k)+1:].split("_")
                if k in {"apriori_startd_prob", "coverage"}:
                    if c[len(k)+1:] != "None":
                        keys[k] = float(c[len(k)+1:])

    return keys

def find_set(ind, parent):
    if ind == parent[ind]:
        return ind
    parent[ind] = find_set(parent[ind], parent)
    return parent[ind]

def union_sets(a, b, parent, rank):
    a = find_set(a, parent)
    b = find_set(b, parent)
    if a != b:
        if rank[a] < rank[b]:
            a, b = b, a
        parent[b] = a
        if rank[a] == rank[b]:
            rank[a] += 1

def in_one_cluster(args):
    s1, s2, comparison_res, ind = args[0], args[1], args[2], args[3]
    in_one = str(s1) in str(s2) if len(s1) < len(s2) else str(s2) in str(s1)
    comparison_res[ind] = in_one or aai([s2, s1]) > IDENTITY

def divide_into_clusters(orfs, t):
    clusters = []
    parent = [i for i in range(len(orfs))]
    rank = [0 for _ in range(len(orfs))]
    for i in range(len(orfs)):
        comparison_res = [False for _ in range(i + 1, len(orfs))]
        Parallel(n_jobs=t, require='sharedmem')(delayed(in_one_cluster)([orfs[i].seq, orfs[j].seq, comparison_res, j - i - 1]) for j in range(i + 1, len(orfs)))
        for j in range(i + 1, len(orfs)):
            if comparison_res[j - i - 1]:
                union_sets(i, j, parent, rank)

    clusters_id = {}
    for i in range(len(orfs)):
        ind = find_set(i, parent)
        if ind not in clusters_id:
            clusters_id[ind] = {}
        clusters_id[ind][orfs[i].name] = orfs[i]

    for cl in clusters_id:
        clusters.append(clusters_id[cl])

    return clusters

def pick_representatives(clusters, orfs, graph):
    edges, coverage = graph[0], graph[1]
    res = []
    cl_ind = 0
    for cl in clusters:
        cl_edges = set()
        cl_meta = []
        for c in cl:
            o = cl[c]
            keys = get_metainfo(o.name)
            path = keys["Edges"]
            cl_meta.append({"name": o.name, "len": len(o.seq), \
                            "start_prob": keys["apriori_startd_prob"] if keys["apriori_startd_prob"] != "None" else 0,\
                            "coverage": keys["coverage"],
                            "edges": set(path)})
            for e in path:
                if len(edges[e]) > MIN_RELIABLE_LENGTH:
                    cl_edges.add(e)
        sorted_edges = sorted([{"name": e, "len": len(edges[e])} for e in cl_edges], key=lambda x: -x["len"])
        cl_edges = set([sorted_edges[i]["name"] for i in range(len(sorted_edges)) if i < 5 or sorted_edges[i]["len"] > RELIABLE_LENGTH])
        sorted_cl = sorted(cl_meta, key=lambda x: (-x["start_prob"], -x["coverage"], -x["len"]))
        num = 0
        for orf in sorted_cl:
            if len(cl_edges & orf["edges"]) > 0 or len(cl_edges) == 0:
                res.append(make_record(cl[orf["name"]].seq, str(cl_ind) + "|" + cl[orf["name"]].id.split(";")[0], str(cl_ind) + "|" + orf["name"].split(";")[0]))
                num += 1
            cl_edges = cl_edges - orf["edges"]
            if len(cl_edges) == 0:
                break
        logging.info( u'Number of representatives ' + str(num) + u' Total number ' + str(len(sorted_cl)))
        cl_ind += 1
    return res

def cluster_orfs(orfs, graph, t):
    logging.info( u'Number of orfs: ' + str(len(orfs)))
    clusters = divide_into_clusters(orfs, t)
    logging.info( u'Number of clusters: ' + str(len(clusters)))
    reprentatives_cl = pick_representatives(clusters, orfs, graph)
    res_cl = []
    for i in range(len(clusters)):
        for o in clusters[i]:
            orf = clusters[i][o]
            res_cl.append(make_record(orf.seq, str(i) + "|" + orf.id.split(";")[0], str(i) + "|" + orf.name.split(";")[0]))
    return reprentatives_cl, res_cl

def load_yaml(outdir):
    with open(outdir + "/config.yaml", 'r') as stream:
        try:
            res = yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            logging.error(exc)
            exit(-1)
    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter ORFs')
    parser.add_argument('-s', '--orfs', help='fasta-file with ORFs', required=True)
    parser.add_argument('-g', '--graph', help='gfa-file with assembly graph', required=True)
    parser.add_argument('-c', '--contigs', help='fasta-file with contigs', required=False)
    parser.add_argument('-p', '--proteins',  help='list of known genes', required=False)
    parser.add_argument('-o', '--out',  help='output directory', required=True)
    parser.add_argument('-t', '--threads', help='threads number', required=False)
    args = parser.parse_args()
    logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level = logging.DEBUG, filename = args.out + u'/filtering.log')
    t = args.threads
    if t == None:
        t = "1"
    config = load_yaml(args.out)

    MIN_RELIABLE_LENGTH = config["orfs_filtering"]["min_reliable_length"]
    RELIABLE_LENGTH = config["orfs_filtering"]["reliable_length"]
    IDENTITY = config["orfs_filtering"]["identity"]
    ED_THRESHOLD = config["orfs_filtering"]["ed_threshold"]

    orfs = load_fasta(args.orfs)
    orfs_new = translate_orfs(orfs)
    orfs_new = leave_unique(orfs_new)
    logging.info( u'Total orfs number: ' + str(len(orfs_new)))
    save_fasta(args.out + "/orfs_final_total", orfs_new)
    if args.contigs != None:
        contigs = load_fasta(args.contigs)
        orfs = align_with_nucmer(orfs, args.orfs, args.contigs, args.out, t)
        orfs = translate_orfs(orfs)
        orfs = leave_unique(orfs)
        logging.info( u'Graph only orfs number: ' + str(len(orfs)))
        save_fasta(args.out + "/orfs_final_graphonly", orfs)
    else:
        orfs = translate_orfs(orfs)
        orfs = leave_unique(orfs)

    if args.proteins != None:
        known_proteins = load_fasta(args.proteins)
        orfs = leave_unknown(orfs, known_proteins)
        logging.info( u'Novel ORFs: ' + str(len(orfs)))
        save_fasta(args.out + "/orfs_final_novel", orfs)

    repres, clustered_orfs = cluster_orfs(orfs, load_gfa_edges(args.graph), int(t))
    logging.info( u'Most reliable ORFs number: ' + str(len(repres)))
    save_fasta(args.out + "/orfs_final_clustered", clustered_orfs)
    save_fasta(args.out + "/orfs_final_most_reliable", repres)

