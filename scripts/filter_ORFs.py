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

def merge(orfs):
    seq_set1 = set()
    seq_set2 = set()
    for orf in orfs:
        if "pathracer" in orf.name:
            seq_set1.add(orf.name)
        if "spaligner" in orf.name:
            seq_set2.add(orf.name)


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
    if args.proteins != None:
        known_proteins = load_fasta(args.proteins)
        orfs = leave_unknown(orfs, known_proteins)
    logging.info( u'Resulting ORFs: ' + str(len(orfs)))
    save_fasta(args.out, orfs)
    merge(orfs)

