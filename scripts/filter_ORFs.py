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
            if p.seq == orf.seq:
                unknown = False
                break
        if unknown:
            res.append(orf)
    return res

def cluster_orfs(orfs):
    clusters = []
    for o in orfs:
        c_id = -1
        paths = set()
        for n in o.name.split(";"):
            if len(n) > 0:
                lst = n.split("|")
            paths.add(lst[1])
        for i in range(len(clusters)):
            if len(paths & clusters[i][0]) > 0:
                c_id = i
                break
        if c_id == -1:
            clusters.append([paths, set({o})])
        else:
            clusters[c_id][0] = paths & clusters[c_id][0]
            clusters[c_id][1].add(o)

    res = []
    for i in range(len(clusters)):
        for orf in clusters[i][1]:
            res.append(make_record(orf.seq, str(i) + "|" + orf.id, str(i) + "|" + orf.name))
    best_orfs = []
    i = 0
    for c in clusters:
        suffix = "_best"
        after_stopcodon = set()
        for o in c[1]:
            for n in o.name.split(";"):
                if len(n) > 0:
                    lst = n.split("|")
                if lst[-1].endswith("True"):
                    after_stopcodon.add(o)
        if len(after_stopcodon) == 0:
            after_stopcodon = c[1]
        best_prob = 0
        best_orf = None
        has_best = False
        for o in after_stopcodon:
            for n in o.name.split(";"):
                if len(n) > 0:
                    lst = n.split("|")
                prob = 0
                for it in lst:
                    if it.startswith("apriori_startd_prob=") and not it.endswith("None"):
                        c_prob = float(it[len("apriori_startd_prob="):])
                        if c_prob > prob:
                            prob = c_prob
            if prob > best_prob:
                best_orf = o
                has_best = True
        if not has_best:
            suffix = "_longest"
            for o in after_stopcodon:
                if not has_best or len(o.seq) > len(best_orf.seq):
                    best_orf = o
                    has_best = True
        best_orfs.append(make_record(best_orf.seq, str(i) + suffix + "|" + best_orf.id, str(i) + suffix + "|"  + best_orf.name))
        i += 1
    return best_orfs, res


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
    if args.contigs != None:
        contigs = load_fasta(args.contigs)
        logging.info( u'ORFs num: ' + str(len(orfs)) + u' Contigs num: ' + str(len(contigs)))
        orfs = align_with_nucmer(orfs, args.orfs, args.contigs, args.out, t)
    orfs = translate_orfs(orfs)
    orfs = leave_unique(orfs)
    # if args.proteins != None:
    #     known_proteins = load_fasta(args.proteins)
    #     orfs = leave_unknown(orfs, known_proteins)
    logging.info( u'Resulting ORFs: ' + str(len(orfs)))
    best_orfs, orfs = cluster_orfs(orfs)
    logging.info( u'Resulting clusters: ' + str(len(best_orfs)))
    save_fasta(args.out, orfs)
    save_fasta(args.out + "_best_in_cluster", best_orfs)

