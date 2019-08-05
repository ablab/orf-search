#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from joblib import Parallel, delayed

import sys
import argparse
from collections import deque
import logging

import os
from os import listdir
from os.path import isfile, isdir, join

import multiprocessing

import load_mappings


sys.setrecursionlimit(1000000)

MET = set({"ATG", "GTG", "TTG"})
cMet = set({"CAT", "CAC", "CAA"})

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

def revert2(edge, pos, edges):
    if edge.endswith("-"):
        c_edge = edge[:-1] + '+'
    else:
        c_edge = edge[:-1] + '-'
    c_pos = len(edges[edge]) - pos
    return c_edge, c_pos

def translate_str(s):
    seq = Seq(s)
    return str(seq.translate())

class Graph:
    def __init__(self, edges, graph, coverage, paths, edge_paths):
        self.edges = edges
        self.graph = graph
        self.coverage = coverage
        self.paths = paths
        self.edge_paths = edge_paths


def load_gfa_edges(gfa_filename):
    res = {}
    graph = {}
    coverage = {}
    paths = {}
    edge_paths = {}
    rev = {"+": "-", "-": "+"}
    with open(gfa_filename, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith("S"):
                lst = ln.strip().split("\t")[1:]
                node_id, seq, kc = lst[0], lst[1], lst[-1] 
                res[node_id + "+"] = seq
                res[node_id + "-"] = make_rc(seq)
                graph[node_id + "+"] = {}
                graph[node_id + "-"] = {}
                coverage[node_id + "+"] = float(kc[len("KC:i:"):])/len(seq)
                coverage[node_id + "-"] = float(kc[len("KC:i:"):])/len(seq)
            elif ln.startswith("L"):
                _, node_id1, pos1, node_id2, pos2, match  = ln.strip().split("\t")
                K = int(match[:-1])
                graph[node_id1+pos1][node_id2+pos2] = 1
                graph[node_id2+rev[pos2]][node_id1+rev[pos1]] = 1
            elif ln.startswith("P"):
                _, name, nodes_lst, _ = ln.strip().split("\t")
                paths[name] = nodes_lst.split(",")
                for n in paths[name]:
                    if n not in edge_paths:
                        edge_paths[n] = []
                    edge_paths[n].append(name)
                path_rc = []
                for n in paths[name]:
                    n_rc = n[:-1] + ("+" if n[-1] == "-" else "-")
                    if n_rc not in edge_paths:
                        edge_paths[n_rc] = []
                    edge_paths[n_rc].append(name + "_rc")
                    path_rc.append(n_rc)
                paths[name + "_rc"] = path_rc[::-1]
    return Graph(res, graph, coverage, paths, edge_paths)


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

def check_for_stopcodon(graph, edges, p, edge, find_set, stop_set, max_length):
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
    found_StopCodon = False
    while len(queue) > 0:
        cur_state = queue.popleft()
        cur_str = cur_state["prev"] + edges[cur_state["edge"]][cur_state["pos"]]
        if len(cur_str) == 3 and cur_str in stop_set:
            found_StopCodon = True

        if len(cur_str) == 3:
            prev = ""
        else:
            prev = cur_str
        if (len(cur_str) != 3 or cur_str not in find_set) and cur_state["dist"] < max_length:
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
    return found_StopCodon

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

def find_start_codons(graph, edges, s_p, s_edge, max_length, only_longest):
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
    start_codon_pos_filtered = []
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

        if check_for_stopcodon(graph, edges, cs_p, e, cStop_codons, cMet, max_length):
            p["after_stopcodon"] = True
            start_codon_pos_filtered.append(p)
        else:
            p["after_stopcodon"] = False
            if not only_longest:
                start_codon_pos_filtered.append(p)


    return start_codon_pos_filtered

def find_stop_codons(graph, edges, f_p, f_edge, max_length):
    stop_codon_pos = find_ends(graph, edges, f_p - 3 + 1, f_edge, Stop_codons, Stop_codons, max_length)
    return stop_codon_pos


def find_paths(graph, edges, s_p, f_p, inner_path, startcodon_dist, only_longest):
    s_edge = inner_path[0]
    f_edge = inner_path[-1]

    if len(startcodon_dist) == 0:
        max_length = 3000
    else:
        max_length = 3*max(startcodon_dist) + 300

    start_codon_pos = find_start_codons(graph, edges, s_p, s_edge, max_length, only_longest)

    stop_codon_pos = find_stop_codons(graph, edges, f_p, f_edge, 3000)
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

def find_connected_edges_tagged(cur_edge, cur_frame, graph, edges, color):
    color.add(cur_edge + "_" + str(cur_frame))
    q = [[cur_edge, cur_frame]]
    l = 0
    while len(q) - l > 0:
        c_edge, c_frame = q[l][0], q[l][1]
        l += 1
        for a_edge in graph[c_edge]:
            if a_edge + "_" + str((c_frame + len(edges[a_edge]) - K) % 3) not in color:
                color.add(a_edge + "_" + str((c_frame + len(edges[a_edge]) - K) % 3))
                q.append([a_edge, (c_frame + len(edges[a_edge]) - K) % 3])

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

def search_all_path_tagged(cur_edge, s_pos, e_pos, cur_len, final_edge, path, paths, all_paths, cur_edges, max_path_num, max_length, min_length, edges_intersection, edges, graph):
    if cur_edges[cur_edge] == 2:
        return
    cur_edges[cur_edge] += 1
    for a_edge in graph[cur_edge]:
        if len(all_paths) < max_path_num:
            if a_edge + "_" + str((cur_len + len(edges[a_edge]) - K )%3) in edges_intersection or\
                a_edge == final_edge and (cur_len + e_pos + 1) % 3 == 0: 
                new_path = path[:]
                new_path.append(a_edge)
                if a_edge != final_edge:
                    if cur_len + len(edges[a_edge]) - K < max_length:
                        search_all_path_tagged(a_edge, s_pos, e_pos, cur_len + len(edges[a_edge]) - K, final_edge,\
                                        new_path, paths, all_paths, cur_edges, max_path_num, max_length, min_length, edges_intersection, edges, graph)
                else:
                    path_len = restore_path_len(edges, s_pos, e_pos, new_path)
                    all_paths.append(path_len)
                    if path_len < max_length and path_len%3 == 0 and path_len > min_length:
                        paths.append([new_path, path_len/3])
    cur_edges[cur_edge] -= 1
    return

def find_subgraph(graph, edges, s_edge, f_edge):
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
    return edges_intersection

def find_subgraph_tagged(graph, edges, s_edge, s_pos, f_edge, f_pos):
    edges_intersection = set()
    color = set()
    if s_pos > len(edges[s_edge]) - K:
        start = 3 - (s_pos - (len(edges[s_edge]) - K)) % 3
    else:
        start = ((len(edges[s_edge]) - K) - s_pos) % 3
    find_connected_edges_tagged(s_edge, start, graph, edges, color)
    connected_to_start = color.copy()
    color = set()
    cf_edge, cf_pos = revert2(f_edge, f_pos, edges)
    if cf_pos > len(edges[cf_edge]) - K:
        start = 3 - (cf_pos - (len(edges[cf_edge]) - K)) % 3
    else:
        start = ((len(edges[cf_edge]) - K) - cf_pos) % 3
    find_connected_edges_tagged(cf_edge, start, graph, edges, color)
    connected_to_finish = color
    for e in connected_to_start:
        c_e = revert(e.split("_")[0])
        if c_e == cf_edge:
            if (f_pos + 1 < len(edges[c_e]) - K and int(e.split("_")[1]) == (len(edges[c_e]) - K - f_pos - 1)%3) \
                or (f_pos + 1 >= len(edges[c_e]) - K and (int(e.split("_")[1]) + f_pos + 1 - (len(edges[c_e]) - K))%3 == 0):
                edges_intersection.add(e.split("_")[0] + "_0")
        elif c_e + "_" + str( (3 + int(e.split("_")[1]) - (len(edges[c_e]) - 2*K)%3 ) % 3 ) in connected_to_finish:
            edges_intersection.add(e)
    return edges_intersection

def generate_all_paths(graph, edges, s_edge, f_edge, s_pos, e_pos, max_path_num, max_length = 1500, min_length = 0):
    edges_intersection = find_subgraph(graph, edges, s_edge, f_edge)
    # edges_intersection_tagged = find_subgraph_tagged(graph, edges, s_edge, s_pos, f_edge, e_pos)
    paths = []
    cur_edges = {}
    for e in edges:
        cur_edges[e] = 0
    all_paths = []
    if s_pos > len(edges[s_edge]) - K:
        start = 3 - (s_pos - (len(edges[s_edge]) - K)) % 3
    else:
        start = ((len(edges[s_edge]) - K) - s_pos)
    # search_all_path_tagged(s_edge, s_pos, e_pos, start, f_edge, [s_edge], paths, all_paths, cur_edges, max_path_num, max_length, min_length, edges_intersection_tagged, edges, graph)
    search_all_path(s_edge, s_pos, e_pos, 0, f_edge, [s_edge], paths, all_paths, cur_edges, max_path_num, max_length, min_length, edges_intersection, edges, graph)
    # print([max_length, min_length, "paths", len(paths), len(all_paths), max_path_num])
    return paths, len(all_paths) 

def get_coverage(path, cov):
    res = []
    for e in path:
        res.append(cov[e])
    return sorted(res)[len(res)/2]

def find_all_paths(graph, edges, coverage, aln, start_codons, stop_codons, startcodon_dist):
    s_edge = aln["path"][0]
    f_edge = aln["path"][-1]
    max_path_num = 5000
    paths = []
    start_codons_paths = {}
    for s in start_codons:
        if len(s["path"][0]) > 0:
            k = s["path"][0][0] + "_" + str(s["path"][1])
            start_codons_paths[k] = {}
            start_codons_paths[k]["paths"], start_codons_paths[k]["is_all_paths"] = \
                generate_all_paths(graph, edges, s["path"][0][0], s_edge, s["path"][1], aln["start"]-1, max_path_num, s["d"]*3 + 1000, s["d"]*3 - 10)
            start_codons_paths[k]["has_sd"] = s["has_sd"]
            start_codons_paths[k]["after_stopcodon"] = s["after_stopcodon"]
        else:
            start_codons_paths[s_edge + "_" + str(s["path"][1])] = {"paths": [[[s_edge],s["d"]]], "is_all_paths": 1, \
                                                                   "has_sd": s["has_sd"], "after_stopcodon": s["after_stopcodon"]}
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
                    score = None
                    if len(startcodon_dist) > 0:
                        score = 0
                        for d in startcodon_dist:
                            if abs(start_path[1] - d) < 50:
                                score += 1.0
                        score /= len(startcodon_dist)
                    path_str = restore_path(edges, start_pos, end_pos, path)
                    if len(path_str) % 3 == 0:
                        med_cov = get_coverage(path, coverage)
                        t_s = translate_str(path_str)
                        if t_s.find("*") == len(t_s) - 1:
                            paths.append({"Edges": path, "seq": path_str, "prefix": aln["prefix"],\
                                          "apriori_startd_prob": score, "starts_cnt": start_num, "stops_cnt": stop_num,\
                                           "start_cnt": total_num_start, "stop_cnt": total_num_end, \
                                           "generated_all": not not_all, "cur_paths_cnt": (start_all*stop_all), \
                                           "has_sd": has_sd, "coverage": med_cov,"after_stopcodon": start_codons_paths[sp]["after_stopcodon"]})
    return paths

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


def compare_with_contig_paths(name, paths, g):
    res_path = []
    for p in paths:
        contigs = set()
        for e in p["Edges"]:
            if e in g.edge_paths:
                for c in g.edge_paths[e]:
                    contigs.add(c)
        supported_edges = set()
        for c in contigs:
            overlap_len = overlap(p["Edges"], g.paths[c], g)
            has_good_overlap = supported_by_contig(p["Edges"], g.paths[c], g, overlap_len)
            if has_good_overlap:
                supported_edges |= set(g.paths[c])
        if len(set(p["Edges"]) - supported_edges) == 0:
            logging.debug(u'Supported by contigs')
            res_path.append(p)
        else:
            unsupported = set(p["Edges"]) - supported_edges
            still_good_path = True
            for i in range(len(p["Edges"])):
                if p["Edges"][i] in unsupported:
                    if i == 0 or i == len(p["Edges"]) - 1:
                        still_good_path = False
                        break
                    e_p, e, e_n = p["Edges"][i - 1], p["Edges"][i], p["Edges"][i + 1]
                    if len(set(g.edge_paths[e_p]) & set(g.edge_paths[e_n])) != 0 \
                        or len(set(g.edge_paths[e_p]) & set(g.edge_paths[e])) != 0 \
                        or len(set(g.edge_paths[e]) & set(g.edge_paths[e_n])) != 0:
                        still_good_path = False
                        break
            if still_good_path:
                logging.debug(u'Supported by contigs')
                res_path.append(p)
            else:
                logging.debug(u'Not supported by contigs')
    return res_path

def compare_with_contig_paths2(name, paths, g, uniqueedge_len):
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
                        if len(g.graph[e].keys()) > 0 and len(g.graph[revert(g.graph[e].keys()[0])]) > 1 and len(g.graph[revert(e)]) > 1:
                            supported = False
                            break
            if not supported:
                break
        if supported:
            res_path.append(p)
    return res_path


def generate_orf(args):
    aln, g, startcodon_dist, only_longest, uniqueedge = args[0], args[1], args[2], args[3], args[4]
    if aln["name"] != "Cry22_MR":
        logging.debug(aln["name"])
        start_codon_pos, stop_codon_pos = find_paths(g.graph, g.edges, aln["start"], aln["end"], aln["path"], startcodon_dist, only_longest)
        logging.debug(u'Start codon num=' + str(len(start_codon_pos)) + ' Stop codons num=' + str(len(stop_codon_pos)))
        all_paths = find_all_paths(g.graph, g.edges, g.coverage, aln, start_codon_pos, stop_codon_pos, startcodon_dist)
        all_paths = compare_with_contig_paths2(aln["name"], all_paths, g, uniqueedge)
        logging.debug(u'Paths num=' + str(len(all_paths)))
        return {"name": aln["name"], "all_paths": all_paths}
    else:
        return {"name": aln["name"], "all_paths": []}

def is_inside(edges, aln1, aln2):
    path1_s = ",".join(aln1["path"])
    path2_s = ",".join(aln2["path"])
    if path1_s in path2_s and \
      (path1_s != path2_s or \
      (path1_s == path2_s and (aln1["start"] > aln2["start"] or aln1["end"] < aln2["end"]))):
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
        ln = restore_path_len(edges, aln2["start"], aln1["start"] - 1, aln2["path"][:i + 1])
        if ln % 3 == 0:
            return True
        else:
            return False
    else:
        return False

def remove_covered_orfs(edges, alns, startcodon_dists):
    res_alns1 = []
    res_startcodon_dists1 = []
    for i in range(len(alns)):
        is_covered = False
        path_str = restore_path(edges, alns[i]["start"], alns[i]["end"], alns[i]["path"])
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
            if is_inside(edges, res_alns1[i], res_alns1[j]):
                is_covered = True
                break
        if not is_covered:
            res_alns.append(res_alns1[i])
            res_startcodon_dists.append(res_startcodon_dists1[i])
    logging.info( u'Reduced number of anchors from ' + str(len(alns)) + ' to ' + str(len(res_alns)))
    return res_alns, res_startcodon_dists


def generate_orfs(output, output_shortest, alns, g, startcodon_dists, only_longest, uniqueedge, t):
    logging.debug( u'Threads ' + str(t) + u' alns ' + str(len(alns)))
    alns, startcodon_dists = remove_covered_orfs(g.edges, alns, startcodon_dists)
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
    g = load_gfa_edges(args.graph)
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



