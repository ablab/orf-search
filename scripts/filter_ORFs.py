#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import subprocess
import logging

import sys
import argparse

import os
from os import listdir
from os.path import isfile, isdir, join

import re
import edlib

def edist(lst):
    if len(str(lst[0])) == 0:
        return len(str(lst[1]))
    if len(str(lst[1])) == 0:
        return len(str(lst[0]))
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    return result["editDistance"], result["cigar"]

def aai(ar):
    p1, p2 = str(ar[0]), str(ar[1])
    if p1.endswith("*"):
        p1 = p1[:-1]
    if p2.endswith("*"):
        p2 = p2[:-1]
    ed, cigar = edist([str(p1), str(p2)])
    matches = re.findall(r'\d+=', cigar)
    aai = 0.0
    for m in matches:
        aai += int(m[:-1])
    aai /= max(len(p1), len(p2))
    return aai*100
    
def load_fasta(filename):
    record_lst = list(SeqIO.parse(filename, "fasta"))
    return record_lst

def save_fasta(filename, orfs):
    with open(filename + ".fasta", "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

def align_with_nucmer(orfs, orfs_filename, contigs_filename, outfile, threads):
    res = []
    com = "nucmer -t " + threads + " --sam-short " + outfile + ".sam " + contigs_filename + " " + orfs_filename
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
        seq_set.add(orf.seq)

    for s in seq_set:
        sid = []
        name = []
        desc = []
        for orf in orfs:
            if s == orf.seq:
                name.append(orf.name)
                sid.append(orf.id)
                desc.append(orf.description)
        res.append(make_record(s, ";".join(name), ";".join(sid), ";".join(desc)))
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

def cluster_orfs_new(orfs):
    clusters = []
    for i in range(len(orfs)):
        new_clusters = []
        first_cluster = -1
        for cl in clusters:
            in_cluster = False
            for o in cl:
                if aai([o.seq, orfs[i].seq]) > 90:
                    #print o.id, orfs[i].id, aai([o.seq, orfs[i].seq])
                    in_cluster = True
                    break
            if in_cluster:
                if first_cluster == -1:
                    new_clusters.append(cl)
                    new_clusters[-1].add(orfs[i])
                    first_cluster = len(new_clusters) - 1
                else:
                    new_clusters[first_cluster] |= cl
            else:
                new_clusters.append(cl)
        if first_cluster == -1:
            new_clusters.append(set({orfs[i]}))
        clusters = []
        for cl in new_clusters:
            clusters.append(cl)
    print len(clusters)
    new_clusters = []
    for i in range(len(clusters)):
        is_covered = False
        for j in range(len(clusters)):
            if i != j and not is_covered:
                for it_i in clusters[i]:
                    for it_j in clusters[j]:
                        if it_i.seq in it_j.seq:
                            is_covered = True
                            break
                    if is_covered:
                        break
            else:
                cur_cluster = set()
                for it_i in clusters[i]:
                    is_covered2 = False
                    for it_j in clusters[j]:
                        if len(it_i.seq) < len(it_j.seq) and it_i.seq in it_j.seq:
                            is_covered2 = True
                            break
                    if not is_covered2:
                        cur_cluster.add(it_i)
        if not is_covered:
            new_clusters.append(cur_cluster)
    res = []
    for i in range(len(new_clusters)):
        for orf in new_clusters[i]:
            res.append(make_record(orf.seq, str(i) + "|" + orf.id.split(";")[0], str(i) + "|" + orf.name.split(";")[0]))
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter ORFs')
    parser.add_argument('-s', '--orfs', help='fasta-file with ORFs', required=True)
    parser.add_argument('-c', '--contigs', help='fasta-file with contigs', required=False)
    parser.add_argument('-p', '--proteins',  help='list of known genes', required=False)
    parser.add_argument('-o', '--out',  help='output prefix', required=True)
    parser.add_argument('-t', '--threads', help='threads number', required=False)
    args = parser.parse_args()
    logging.basicConfig(format = u'%(levelname)-8s [%(asctime)s] %(message)s', level = logging.DEBUG, filename = args.out + u'.log')
    t = args.threads
    if t == None:
        t = "1"

    orfs = load_fasta(args.orfs)
    orfs_new = translate_orfs(orfs)
    orfs_new = leave_unique(orfs_new)
    logging.info( u'Total orfs number: ' + str(len(orfs_new)))
    save_fasta(args.out + "_total", orfs_new)
    if args.contigs != None:
        contigs = load_fasta(args.contigs)
        logging.info( u'ORFs num: ' + str(len(orfs)) + u' Contigs num: ' + str(len(contigs)))
        orfs = align_with_nucmer(orfs, args.orfs, args.contigs, args.out, t)
        orfs = translate_orfs(orfs)
        orfs = leave_unique(orfs)
        logging.info( u'Graph only orfs number: ' + str(len(orfs)))
        save_fasta(args.out + "_graphonly", orfs)
    else:
        orfs = translate_orfs(orfs)
        orfs = leave_unique(orfs)

    if args.proteins != None:
        known_proteins = load_fasta(args.proteins)
        orfs = leave_unknown(orfs, known_proteins)

    logging.info( u'Resulting ORFs: ' + str(len(orfs)))
    save_fasta(args.out + "_novel", orfs)
    clustered_orfs = cluster_orfs_new(orfs)
    logging.info( u'Resulting clusters: ' + str(len(clustered_orfs)))
    save_fasta(args.out + "_novel_clustered", clustered_orfs)

