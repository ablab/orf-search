#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import sys
import argparse
from collections import deque

import os
from os import listdir
from os.path import isfile, isdir, join

import multiprocessing

import load_mappings

sys.setrecursionlimit(1000000)

MET = set({"ATG"})
cMet = set({"CAT"})

sd_seq = "AGGAGG"
asd_seq = "CCTCCT"

Stop_codons = set({"TAA", "TAG", "TGA"})
cStop_codons = set({"TTA", "CTA", "TCA"})
K = 55

def make_rc(s):
    seq = Seq(s, generic_dna)
    return str(seq.reverse_complement())

def revert(edge):
    if edge.endswith("-"):
        c_edge = edge[:-1] + '+'
    else:
        c_edge = edge[:-1] + '-'
    return c_edge

def translate_str(s):
    seq = Seq(s)
    return str(seq.translate())

def load_gfa_edges(gfa_filename):
    res = {}
    graph = {}
    rev = {"+": "-", "-": "+"}
    with open(gfa_filename, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith("S"):
                node_id, seq = ln.strip().split("\t")[1:3]
                res[node_id + "+"] = seq
                res[node_id + "-"] = make_rc(seq)
                graph[node_id + "+"] = {}
                graph[node_id + "-"] = {}
            elif ln.startswith("L"):
                _, node_id1, pos1, node_id2, pos2, match  = ln.strip().split("\t")
                graph[node_id1+pos1][node_id2+pos2] = 1
                graph[node_id2+rev[pos2]][node_id1+rev[pos1]] = 1
    return res, graph


def load_subpaths(filename):
    alns = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            name, start_pos, end_pos, path = ln.strip().split("\t")
            alns.append({"name": name, "start": int(start_pos), "end": int(end_pos), "path": paths.split(",")})
    return alns

def check_for_sdseq(graph, edges, p, edge, max_length = 20):
    queue = deque()
    color = {}
    queue.append({"edge": edge, "prev": "", "pos": p, "dist": 0})
    if p == len(edges[edge]):
        queue.append({"edge": edge, "prev": "", "pos": p - 1, "dist": 0})  
        start_state = edge + "_" + str(p - 1) + "_"
    else:
        queue.append({"edge": edge, "prev": "", "pos": p, "dist": 0})
        start_state = edge + "_" + str(p) + "_"

    found_SD = False
    while len(queue) > 0:
        cur_state = queue.popleft()
        cur_str = cur_state["prev"] + edges[cur_state["edge"]][cur_state["pos"]]
        if cur_str == asd_seq:
            found_SD = True
            break
        if len(cur_str) == 6:
            prev = cur_str[1:]
        else:
            prev = cur_str

        if cur_state["dist"] < max_length:
            if cur_state["pos"] + 1 > len(edges[cur_state["edge"]]) - 1:
                for a_edge in graph[cur_state["edge"]]:
                    if a_edge + "_" + str(K)  + "_" + prev not in color: 
                        queue.append({"edge": a_edge, "prev": prev, "pos": K, "dist": cur_state["dist"] + 1})
                        color[a_edge + "_" + str(K) + "_" + prev] = cur_state["edge"] + "_" + str(cur_state["pos"])  + "_" + cur_state["prev"] 
            else:
                if cur_state["edge"] + "_" + str(cur_state["pos"] + 1)  + "_" + prev not in color:
                    queue.append({"edge": cur_state["edge"], "prev": prev, "pos": cur_state["pos"] + 1, "dist": cur_state["dist"] + 1})
                    color[cur_state["edge"] + "_" + str(cur_state["pos"] + 1)  + "_" + prev] = \
                                    cur_state["edge"] + "_" + str(cur_state["pos"])  + "_" + cur_state["prev"]
                                    
    return found_SD

def restore_subpath(color, s, f):
    path = []
    pos = int(f.split("_")[1])
    c = f
    while c != s:
        if c.split("_")[0] != color[c].split("_")[0] or \
           (c.split("_")[0] ==  color[c].split("_")[0] and int(c.split("_")[1]) < int(color[c].split("_")[1])) :
            path.append(c.split("_")[0])
        c = color[c]
    path = path[::-1]
    return [path, pos]

def find_ends(graph, edges, p, edge, find_set, stop_set, max_length):
    queue = deque()
    color = {}
    if p == len(edges[edge]):
        queue.append({"edge": edge, "prev": "", "pos": p - 3, "dist": 0})  
        start_state = edge + "_" + str(p - 3) + "_" 
    else:
        queue.append({"edge": edge, "prev": "", "pos": p, "dist": 0})
        start_state = edge + "_" + str(p) + "_"
    final_state = None
    res = []
    cnt = 0
    while len(queue) > 0:
        cur_state = queue.popleft()
        #print([cur_state["edge"], len(edges[cur_state["edge"]]), cur_state["pos"], cur_state])
        cur_str = cur_state["prev"] + edges[cur_state["edge"]][cur_state["pos"]]
        if len(cur_str) == 3 and cur_str in find_set:
            final_state = cur_state["edge"] + "_" + str(cur_state["pos"]) + "_" + cur_state["prev"] 
            res.append({"path": restore_subpath(color, start_state, final_state), "d": cur_state["dist"]/3})

        if len(cur_str) == 3:
            prev = ""
        else:
            prev = cur_str
        if (len(cur_str) != 3 or cur_str not in stop_set) and cur_state["dist"] < max_length:
            if cur_state["pos"] + 1 > len(edges[cur_state["edge"]]) - 1:
                for a_edge in graph[cur_state["edge"]]:
                    if a_edge + "_" + str(K) + "_" + prev not in color: 
                        queue.append({"edge": a_edge, "prev": prev, "pos": K, "dist": cur_state["dist"] + 1})
                        color[a_edge + "_" + str(K) + "_" + prev] = cur_state["edge"] + "_" + str(cur_state["pos"]) + "_" + cur_state["prev"]
            else:
                if cur_state["edge"] + "_" + str(cur_state["pos"] + 1) + "_" + prev not in color:
                    queue.append({"edge": cur_state["edge"], "prev": prev, "pos": cur_state["pos"] + 1, "dist": cur_state["dist"] + 1})
                    color[cur_state["edge"] + "_" + str(cur_state["pos"] + 1) + "_" + prev] = \
                            cur_state["edge"] + "_" + str(cur_state["pos"]) + "_" + cur_state["prev"]
    return res

def restore_path(edges, s_p, f_p, path):
    if len(path) == 1:
        return edges[path[0]][s_p:f_p + 1]

    cur_pos = s_p
    res = edges[path[0]][cur_pos:len(edges[path[0]]) - K]
    prev = path[0]
    for p in path[1:-1]:
        cur_pos -= (len(edges[prev]) - K)
        res += edges[p][max(cur_pos, 0):len(edges[p]) - K]
        prev = p
    cur_pos -= (len(edges[prev]) - K)
    res += edges[path[-1]][max(cur_pos, 0):f_p + 1]
    return res

def restore_path_len(edges, s_p, f_p, path):
    if len(path) == 1:
        return f_p - s_p + 1
    cur_pos = s_p
    res = max(0, len(edges[path[0]]) - K - cur_pos)
    prev = path[0]
    for p in path[1:-1]:
        cur_pos -= (len(edges[prev]) - K)
        res += max(0, len(edges[p]) - K - max(cur_pos, 0))
        prev = p
    cur_pos -= (len(edges[prev]) - K)
    res += f_p + 1 - max(cur_pos, 0)
    return res

def find_start_codons(graph, edges, s_p, s_edge, max_length):
    cs_edge = revert(s_edge)
    if len(edges[cs_edge]) - (s_p + 3) > 0:
        cs_p = len(edges[cs_edge]) - (s_p + 3)
    else:
        cs_p = len(edges[cs_edge]) - s_p
    cstart_codon_pos = find_ends(graph, edges, cs_p, cs_edge, cMet, cStop_codons, max_length)
    start_codon_pos = []
    for potential_start_pos in cstart_codon_pos:
        cpos = potential_start_pos["path"][1]
        cpath = potential_start_pos["path"][0] 
        if len(cpath) > 0:
            prefix = []
            for e in cpath:
                if e.endswith("-"):
                    prefix.append(e[:-1] + '+')
                else:
                    prefix.append(e[:-1] + '-')
            prefix = prefix[::-1]
            start_codon_pos.append({"path": [prefix, len(edges[prefix[0]]) - cpos - 1], "d": potential_start_pos["d"]})
        else:
            pos = len(edges[cs_edge]) - cpos - 1
            start_codon_pos.append({"path":[[], pos], "d": potential_start_pos["d"]})
    for p in start_codon_pos:
        if len(p["path"][0]) == 0:
            e = revert(s_edge)
        else:
            e = revert(p["path"][0][0])

        if len(edges[e]) - p["path"][1] > 0:
            cs_p = len(edges[e]) - (p["path"][1] + 1)
        else:
            cs_p = len(edges[e]) - p["path"][1]

        if check_for_sdseq(graph, edges, cs_p, e):
            p["has_sd"] = True
        else:
            p["has_sd"] = False

    return start_codon_pos

def find_stop_codons(graph, edges, f_p, f_edge, max_length):
    stop_codon_pos = find_ends(graph, edges, f_p + 1, f_edge, Stop_codons, Stop_codons, max_length)
    return stop_codon_pos


def find_paths(graph, edges, s_p, f_p, inner_path, startcodon_dist):
    s_edge = inner_path[0]
    f_edge = inner_path[-1]

    if len(startcodon_dist) == 0:
        max_length = 3000
    else:
        max_length = 3*max(startcodon_dist) + 300

    start_codon_pos = []
    for i in range(3):
        start_codon_pos.extend(find_start_codons(graph, edges, s_p + i, s_edge, max_length))

    stop_codon_pos = []
    for i in range(3):
        stop_codon_pos.extend(find_stop_codons(graph, edges, f_p - i, f_edge, 3000))
    return start_codon_pos, stop_codon_pos


def find_connected_edges(cur_edge, graph, color):
    color.add(cur_edge)
    q = [cur_edge]
    l = 0
    while len(q) - l > 0:
        c_edge = q[l]
        l += 1
        for a_edge in graph[c_edge]:
            if a_edge not in color:
                color.add(a_edge)
                q.append(a_edge)

def search_all_path(cur_edge, s_pos, e_pos, cur_len, final_edge, path, paths, all_paths, cur_edges, max_path_num, max_length, min_length, edges_intersection, edges, graph):
    if cur_edges[cur_edge] == 2:
        return
    cur_edges[cur_edge] += 1

    for a_edge in graph[cur_edge]:
        if a_edge in edges_intersection and len(all_paths) < max_path_num:
            new_path = path[:]
            new_path.append(a_edge)
            if a_edge != final_edge and cur_len + len(edges[a_edge]) - K < max_length:
                search_all_path(a_edge, s_pos, e_pos, cur_len + len(edges[a_edge]) - K, final_edge,\
                                new_path, paths, all_paths, cur_edges, max_path_num, max_length, min_length, edges_intersection, edges, graph)
            else:
                path_len = restore_path_len(edges, s_pos, e_pos, new_path)
                all_paths.append(path_len)
                if path_len < max_length and path_len%3 == 0 and path_len > min_length:
                    paths.append([new_path, path_len/3])
    cur_edges[cur_edge] -= 1
    return

def generate_all_paths(graph, edges, s_edge, f_edge, s_pos, e_pos, max_path_num, max_length = 1500, min_length = 0):
    edges_intersection = set()
    color = set()
    find_connected_edges(s_edge, graph, color)
    connected_to_start = color.copy()
    color = set()
    cf_edge = revert(f_edge)
    find_connected_edges(cf_edge, graph, color)
    connected_to_finish = color
    for e in connected_to_start:
        c_e = revert(e)
        if c_e in connected_to_finish:
            edges_intersection.add(e)
    paths = []
    cur_edges = {}
    for e in edges:
        cur_edges[e] = 0
    all_paths = []
    search_all_path(s_edge, s_pos, e_pos, 0, f_edge, [s_edge], paths, all_paths, cur_edges, max_path_num, max_length, min_length, edges_intersection, edges, graph)
    print([max_length, min_length, "paths", len(paths), len(all_paths), max_path_num])
    return paths, len(all_paths) 

def find_all_paths(graph, edges, aln, start_codons, stop_codons, startcodon_dist):
    s_edge = aln["path"][0]
    f_edge = aln["path"][-1]
    max_path_num = 5000
    paths = []
    start_codons_paths = {}
    #print("Find prefix")
    for s in start_codons:
        if len(s["path"][0]) > 0:
            k = s["path"][0][0] + "_" + str(s["path"][1])
            start_codons_paths[k] = {}
            start_codons_paths[k]["paths"], start_codons_paths[k]["is_all_paths"] = \
                generate_all_paths(graph, edges, s["path"][0][0], s_edge, s["path"][1], aln["start"]-1, max_path_num, s["d"]*3 + 1000, s["d"]*3 - 10)
            start_codons_paths[k]["has_sd"] = s["has_sd"]
        else:
            start_codons_paths[s_edge + "_" + str(s["path"][1])] = {"paths": [[[s_edge],s["d"]]], "is_all_paths": 1, "has_sd": s["has_sd"]}
    #print("Find suffix")
    stop_codons_paths = {}
    for s in stop_codons:
        if len(s["path"][0]) > 0:
            k = s["path"][0][-1] + "_" + str(s["path"][1])
            stop_codons_paths[k] = {}
            stop_codons_paths[k]["paths"], stop_codons_paths[k]["is_all_paths"] = \
                generate_all_paths(graph, edges, f_edge, s["path"][0][-1], aln["end"] + 1, s["path"][1], max_path_num, s["d"]*3 + 1000, s["d"]*3 - 10)
        else:
            stop_codons_paths[f_edge + "_" + str(s["path"][1])] = {"paths": [[[f_edge],s["d"]]], "is_all_paths": 1}

    total_num_start, total_num_end = 0, 0
    for sp in start_codons_paths.keys():
        for ep in stop_codons_paths.keys():
            start_num, stop_num = len(start_codons_paths[sp]["paths"]), len(stop_codons_paths[ep]["paths"]) 
            total_num_start, total_num_end = total_num_start + start_num, total_num_end + stop_num

    for sp in start_codons_paths.keys():
        for ep in stop_codons_paths.keys():
            has_sd = start_codons_paths[sp]["has_sd"]
            start_all, stop_all = start_codons_paths[sp]["is_all_paths"], stop_codons_paths[ep]["is_all_paths"]
            not_all = start_all >= max_path_num or stop_all >= max_path_num
            start_pos, end_pos = int(sp.split("_")[1]), int(ep.split("_")[1])
            start_num, stop_num = len(start_codons_paths[sp]["paths"]), len(stop_codons_paths[ep]["paths"]) 
            for start_path in start_codons_paths[sp]["paths"]:
                for end_path in stop_codons_paths[ep]["paths"]:
                    path = start_path[0][:-1] + aln["path"] + end_path[0][1:]
                    score = 0.5
                    if len(startcodon_dist) > 0:
                        for d in startcodon_dist:
                            if abs(start_path[1] - d) < 50:
                                score += 1.0
                        score /= len(startcodon_dist)
                    path_str = restore_path(edges, start_pos, end_pos, path)
                    if len(path_str) % 3 == 0 and len(path) > 1:
                        t_s = translate_str(path_str)
                        if t_s.find("*") == len(t_s) - 1:
                            paths.append({"Edges": path, "seq": path_str,\
                                          "apriori_startd_prob": score, "starts_cnt": start_num, "stops_cnt": stop_num,\
                                           "start_cnt": total_num_start, "stop_cnt": total_num_end, \
                                           "generated_all": not not_all, "cur_paths_cnt": (start_all*stop_all), \
                                           "has_sd": has_sd})
    return paths


def generate_orf(args):
    aln, graph, edges, startcodon_dist = args[0], args[1], args[2], args[3]
    if aln["name"] != "Cry22_MR":
        print(aln)
        start_codon_pos, stop_codon_pos = find_paths(graph, edges, aln["start"], aln["end"], aln["path"], startcodon_dist)
        print("Start codon num=" + str(len(start_codon_pos)) + " Stop codons num=" + str(len(stop_codon_pos)))
        all_paths = find_all_paths(graph, edges, aln, start_codon_pos, stop_codon_pos, startcodon_dist)
        print("Paths num=" + str(len(all_paths)))
        return {"name": aln["name"], "all_paths": all_paths}
    else:
        return {"name": aln["name"], "all_paths": []}


def generate_orfs(output, output_shortest, alns, edges, graph, startcodon_dists, t, prefix):
    print("Threads " + str(t))
    pool = multiprocessing.Pool(processes = t)
    args = zip(alns, [graph for _ in range(len(alns))], [edges for _ in range(len(alns))], startcodon_dists)
    all_orfs = pool.map(generate_orf, args)

    with open(output, "a+") as fout:    
        for orf in all_orfs:
            name = orf["name"]
            for path in orf["all_paths"]:
                fout.write(">" + prefix + "_" + name + "|Edges=" + "_".join(path["Edges"]) + "|" + \
                           "|".join([k + "=" + str(path[k]) for k in path.keys() if k not in {"Edges", "seq"}]) + "\n")
                fout.write(path["seq"] + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print potential ORFs for all given graph alignments in fasta-format')
    parser.add_argument('-s', '--sequences', help='fasta-file with SPAligner alignments', required=False)
    parser.add_argument('-m', '--hmms', help='path to PathRacer results folder', required=False)
    parser.add_argument('-g', '--graph', help='path to assembly graph (in gfa-format)', required=True)
    parser.add_argument('-k', '--kmer', help='k-mer size in graph', required=True)
    parser.add_argument('-p', '--proteins',  help='list of genes to estimate hmms position on genes (domtbl has to be set)', required=False)
    parser.add_argument('-d', '--domtbl',  help='HMMer alignment of hmms to genes (genes has to be set)', required=False)
    parser.add_argument('-o', '--out',  help='output prefix', required=True)
    parser.add_argument('-t', '--threads', help='threads number', required=False)
    args = parser.parse_args()
    
    if args.hmms == None and args.sequences == None:
        print("Error: Please provide sequences or hmm alignments")
        exit(-1)
    K = int(args.kmer)
    edges, graph = load_gfa_edges(args.graph)
    output = args.out + ".fasta"
    output_shortest = args.out + "_shortest.fasta"
    if os.path.exists(output):
        os.remove(output)
    if os.path.exists(output_shortest):
        os.remove(output_shortest)

    if args.hmms != None:
        if (args.proteins == None or args.domtbl == None):
            print("Warning: information about HMMs positions is not provided")
            hmm_hits = {}
        else:
            hmm_hits = load_mappings.load_true_crygenes(args.proteins, args.domtbl, {}, K)

        filenames = [f for f in listdir(args.hmms) \
                     if isfile(join(args.hmms, f)) and \
                        isfile(join(args.hmms, f[:-9] + "edges.fa")) and \
                        f.endswith("domtblout") ]
        alns = []
        for f in filenames:
            f_path = join(args.hmms, f)
            alns.extend(load_mappings.load_pathracer_mapping(f_path, edges, K))
        startcodon_dists = []
        for aln in alns:
            startcodon_dist = []
            if aln["name"] in hmm_hits:
                startcodon_dist = [x["start"] for x in hmm_hits[aln["name"]] ]
            startcodon_dists.append(startcodon_dist)
        generate_orfs(output, output_shortest, alns, edges, graph, startcodon_dists, int(args.threads), "pathracer")

    if args.sequences != None:
        alns = load_mappings.load_spaligner_mapping(args.sequences)
        startcodon_dists = []
        for a in alns:
            startcodon_dists.append([a["d"]])
        generate_orfs(output, output_shortest, alns, edges, graph, startcodon_dists, int(args.threads), "spaligner")


