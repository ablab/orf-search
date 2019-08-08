import sys

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

sys.setrecursionlimit(1000000)

def make_rc(s):
    seq = Seq(s, generic_dna)
    return str(seq.reverse_complement())

class Graph:
    def __init__(self, gfa_filename):
        K = None
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
        self.edges = res
        self.graph = graph
        self.coverage = coverage
        self.paths = paths
        self.edge_paths = edge_paths
        self.K = K

    def revert_edge(self, edge):
        if edge.endswith("-"):
            c_edge = edge[:-1] + '+'
        else:
            c_edge = edge[:-1] + '-'
        return c_edge

    def revert_position(self, edge, pos):
        if edge.endswith("-"):
            c_edge = edge[:-1] + '+'
        else:
            c_edge = edge[:-1] + '-'
        c_pos = len(self.edges[edge]) - pos
        return c_edge, c_pos

    def get_coverage(self, path):
        res = []
        for e in path:
            res.append(self.coverage[e])
        return sorted(res)[len(res)//2]

    def find_connected_edges(self, cur_edge, color):
        color.add(cur_edge)
        q = [cur_edge]
        l = 0
        while len(q) - l > 0:
            c_edge = q[l]
            l += 1
            for a_edge in self.graph[c_edge]:
                if a_edge not in color:
                    color.add(a_edge)
                    q.append(a_edge)

    def find_subgraph(self, s_edge, f_edge):
        edges_intersection = set()
        color = set()
        self.find_connected_edges(s_edge, color)
        connected_to_start = color.copy()
        color = set()
        cf_edge = self.revert_edge(f_edge)
        self.find_connected_edges(cf_edge, color)
        connected_to_finish = color
        for e in connected_to_start:
            c_e = self.revert_edge(e)
            if c_e in connected_to_finish:
                edges_intersection.add(e)
        return edges_intersection