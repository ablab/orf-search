from Bio.Seq import Seq

from assembly_graph import Graph

def translate_str(s):
        seq = Seq(s)
        return str(seq.translate())

class PathConstructor():
    def __init__(self, g):
        self.g = g

    def restore_path(self, s_p, f_p, path):
        if len(path) == 1:
            return self.g.edges[path[0]][s_p:f_p + 1]

        cur_pos = s_p
        res = self.g.edges[path[0]][cur_pos:len(self.g.edges[path[0]]) - self.g.K]
        prev = path[0]
        for p in path[1:-1]:
            cur_pos -= (len(self.g.edges[prev]) - self.g.K)
            res += self.g.edges[p][max(cur_pos, 0):len(self.g.edges[p]) - self.g.K]
            prev = p
        cur_pos -= (len(self.g.edges[prev]) - self.g.K)
        res += self.g.edges[path[-1]][max(cur_pos, 0):f_p + 1]
        return res

    def restore_path_len(self, s_p, f_p, path):
        if len(path) == 1:
            return f_p - s_p + 1
        cur_pos = s_p
        res = max(0, len(self.g.edges[path[0]]) - self.g.K - cur_pos)
        prev = path[0]
        for p in path[1:-1]:
            cur_pos -= (len(self.g.edges[prev]) - self.g.K)
            res += max(0, len(self.g.edges[p]) - self.g.K - max(cur_pos, 0))
            prev = p
        cur_pos -= (len(self.g.edges[prev]) - self.g.K)
        res += f_p + 1 - max(cur_pos, 0)
        return res

class EndsPathConstructor():
    
    def __init__(self, g, config):
        self.g = g
        self.max_path_num = config["orfs_search"]["max_paths_num"]
        self.max_length = config["orfs_search"]["max_path_length"]
        self.min_length = config["orfs_search"]["min_path_length"]
        self.path_constructor = PathConstructor(g)

    def search_all_path(self, cur_edge, s_pos, e_pos, cur_len, final_edge, path, paths, all_paths, cur_edges, edges_intersection):
        if cur_edges[cur_edge] == 2:
            return
        cur_edges[cur_edge] += 1

        for a_edge in self.g.graph[cur_edge]:
            if a_edge in edges_intersection and len(all_paths) < self.max_path_num:
                new_path = path[:]
                new_path.append(a_edge)
                if a_edge != final_edge and cur_len + len(self.g.edges[a_edge]) - self.g.K < self.max_length:
                    self.search_all_path(a_edge, s_pos, e_pos, cur_len + len(self.g.edges[a_edge]) - self.g.K, final_edge,\
                                    new_path, paths, all_paths, cur_edges, edges_intersection)
                else:
                    path_len = self.path_constructor.restore_path_len(s_pos, e_pos, new_path)
                    all_paths.append(path_len)
                    if path_len < self.max_length and path_len%3 == 0 and path_len > self.min_length:
                        paths.append([new_path, path_len//3])
        cur_edges[cur_edge] -= 1
        return

    def generate_possible_paths(self, s_edge, f_edge, s_pos, e_pos):
        edges_intersection = self.g.find_subgraph(s_edge, f_edge)
        paths = []
        cur_edges = {}
        for e in self.g.edges:
            cur_edges[e] = 0
        all_paths = []
        if s_pos > len(self.g.edges[s_edge]) - self.g.K:
            start = 3 - (s_pos - (len(self.g.edges[s_edge]) - self.g.K)) % 3
        else:
            start = ((len(self.g.edges[s_edge]) - self.g.K) - s_pos)
        self.search_all_path(s_edge, s_pos, e_pos, 0, f_edge, [s_edge], paths, all_paths, cur_edges, edges_intersection)
        return paths, len(all_paths) 

    def generate_full_paths_endings(self, aln, start_codons, stop_codons):
        s_edge = aln["path"][0]
        f_edge = aln["path"][-1]
        start_codons_paths = {}
        for s in start_codons:
            if len(s["path"][0]) > 0:
                k = s["path"][0][0] + "_" + str(s["path"][1])
                self.max_length = s["d"]*3 + 1000
                self.min_length = s["d"]*3 - 10
                start_codons_paths[k] = {}
                start_codons_paths[k]["paths"], start_codons_paths[k]["is_all_paths"] = \
                    self.generate_possible_paths(s["path"][0][0], s_edge, s["path"][1], aln["start"]-1)
                start_codons_paths[k]["after_stopcodon"] = s["after_stopcodon"]
            else:
                start_codons_paths[s_edge + "_" + str(s["path"][1])] = {"paths": [[[s_edge],s["d"]]], "is_all_paths": 1, \
                                                                       "after_stopcodon": s["after_stopcodon"]}
        stop_codons_paths = {}
        for s in stop_codons:
            if len(s["path"][0]) > 0:
                k = s["path"][0][-1] + "_" + str(s["path"][1])
                self.max_length = s["d"]*3 + 1000
                self.min_length = s["d"]*3 - 10
                stop_codons_paths[k] = {}
                stop_codons_paths[k]["paths"], stop_codons_paths[k]["is_all_paths"] = \
                    self.generate_possible_paths(f_edge, s["path"][0][-1], aln["end"] + 1, s["path"][1])
            else:
                stop_codons_paths[f_edge + "_" + str(s["path"][1])] = {"paths": [[[f_edge],s["d"]]], "is_all_paths": 1}

        return start_codons_paths, stop_codons_paths

    def construct_full_paths(self, aln, start_codons_paths, stop_codons_paths, startcodon_dist):
        paths = []

        total_num_start, total_num_end = 0, 0
        for sp in start_codons_paths.keys():
            for ep in stop_codons_paths.keys():
                start_num, stop_num = len(start_codons_paths[sp]["paths"]), len(stop_codons_paths[ep]["paths"]) 
                total_num_start, total_num_end = total_num_start + start_num, total_num_end + stop_num

        for sp in start_codons_paths.keys():
            for ep in stop_codons_paths.keys():
                start_all, stop_all = start_codons_paths[sp]["is_all_paths"], stop_codons_paths[ep]["is_all_paths"]
                not_all = start_all >= self.max_path_num or stop_all >= self.max_path_num
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
                        path_str = self.path_constructor.restore_path(start_pos, end_pos, path)
                        if len(path_str) % 3 == 0:
                            med_cov = self.g.get_coverage(path)
                            t_s = translate_str(path_str)
                            if t_s.find("*") == len(t_s) - 1:
                                paths.append({"Edges": path, "seq": path_str, "prefix": aln["prefix"],\
                                              "apriori_startd_prob": score, "starts_cnt": start_num, "stops_cnt": stop_num,\
                                               "start_cnt": total_num_start, "stop_cnt": total_num_end, \
                                               "generated_all": not not_all, "cur_paths_cnt": (start_all*stop_all), \
                                               "coverage": med_cov,"after_stopcodon": start_codons_paths[sp]["after_stopcodon"]})
        return paths


    def restore_full_paths(self, aln, start_codons, stop_codons, startcodon_dist):
        start_codon_paths, stop_codons_paths = self.generate_full_paths_endings(aln, start_codons, stop_codons)     
        paths = self.construct_full_paths(aln, start_codon_paths, stop_codons_paths, startcodon_dist)
        return paths
