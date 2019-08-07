from collections import deque
from assembly_graph import Graph

MET = set({"ATG", "GTG", "TTG"})
cMet = set({"CAT", "CAC", "CAA"})

Stop_codons = set({"TAA", "TAG", "TGA"})
cStop_codons = set({"TTA", "CTA", "TCA"})



class GeneEndsFinder():

    def __init__(self, g):
        self.graph = g.graph
        self.edges = g.edges
        self.g = g

    def restore_subpath(self, color, s, f):
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

    def find_ends(self, p, edge, find_set, stop_set, max_length):
        queue = deque()
        color = {}
        if p == len(self.edges[edge]):
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
            cur_str = cur_state["prev"] + self.edges[cur_state["edge"]][cur_state["pos"]]
            if len(cur_str) == 3 and cur_str in find_set:
                final_state = cur_state["edge"] + "_" + str(cur_state["pos"]) + "_" + cur_state["prev"] 
                res.append({"path": self.restore_subpath(color, start_state, final_state), "d": cur_state["dist"]/3})

            if len(cur_str) == 3:
                prev = ""
            else:
                prev = cur_str
            if (len(cur_str) != 3 or cur_str not in stop_set) and cur_state["dist"] < max_length:
                if cur_state["pos"] + 1 > len(self.edges[cur_state["edge"]]) - 1:
                    for a_edge in self.graph[cur_state["edge"]]:
                        if a_edge + "_" + str(self.g.K) + "_" + prev not in color: 
                            queue.append({"edge": a_edge, "prev": prev, "pos": self.g.K, "dist": cur_state["dist"] + 1})
                            color[a_edge + "_" + str(self.g.K) + "_" + prev] = cur_state["edge"] + "_" + str(cur_state["pos"]) + "_" + cur_state["prev"]
                else:
                    if cur_state["edge"] + "_" + str(cur_state["pos"] + 1) + "_" + prev not in color:
                        queue.append({"edge": cur_state["edge"], "prev": prev, "pos": cur_state["pos"] + 1, "dist": cur_state["dist"] + 1})
                        color[cur_state["edge"] + "_" + str(cur_state["pos"] + 1) + "_" + prev] = \
                                cur_state["edge"] + "_" + str(cur_state["pos"]) + "_" + cur_state["prev"]
        return res

    def check_for_stopcodon(self, p, edge, find_set, stop_set, max_length):
        queue = deque()
        color = {}
        if p == len(self.edges[edge]):
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
            cur_str = cur_state["prev"] + self.edges[cur_state["edge"]][cur_state["pos"]]
            if len(cur_str) == 3 and cur_str in stop_set:
                found_StopCodon = True

            if len(cur_str) == 3:
                prev = ""
            else:
                prev = cur_str
            if (len(cur_str) != 3 or cur_str not in find_set) and cur_state["dist"] < max_length:
                if cur_state["pos"] + 1 > len(self.edges[cur_state["edge"]]) - 1:
                    for a_edge in self.graph[cur_state["edge"]]:
                        if a_edge + "_" + str(self.g.K) + "_" + prev not in color: 
                            queue.append({"edge": a_edge, "prev": prev, "pos": self.g.K, "dist": cur_state["dist"] + 1})
                            color[a_edge + "_" + str(self.g.K) + "_" + prev] = cur_state["edge"] + "_" + str(cur_state["pos"]) + "_" + cur_state["prev"]
                else:
                    if cur_state["edge"] + "_" + str(cur_state["pos"] + 1) + "_" + prev not in color:
                        queue.append({"edge": cur_state["edge"], "prev": prev, "pos": cur_state["pos"] + 1, "dist": cur_state["dist"] + 1})
                        color[cur_state["edge"] + "_" + str(cur_state["pos"] + 1) + "_" + prev] = \
                                cur_state["edge"] + "_" + str(cur_state["pos"]) + "_" + cur_state["prev"]
        return found_StopCodon

    def find_start_codons(self, s_p, s_edge, max_length, only_longest):
        cs_edge = self.g.revert_edge(s_edge)
        if len(self.edges[cs_edge]) - (s_p + 3) > 0:
            cs_p = len(self.edges[cs_edge]) - (s_p + 3)
        else:
            cs_p = len(self.edges[cs_edge]) - s_p
        cstart_codon_pos = self.find_ends(cs_p, cs_edge, cMet, cStop_codons, max_length)
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
                start_codon_pos.append({"path": [prefix, len(self.edges[prefix[0]]) - cpos - 1], "d": potential_start_pos["d"]})
            else:
                pos = len(self.edges[cs_edge]) - cpos - 1
                start_codon_pos.append({"path":[[], pos], "d": potential_start_pos["d"]})
        start_codon_pos_filtered = []
        for p in start_codon_pos:
            if len(p["path"][0]) == 0:
                e = self.g.revert_edge(s_edge)
            else:
                e = self.g.revert_edge(p["path"][0][0])

            if len(self.edges[e]) - p["path"][1] > 0:
                cs_p = len(self.edges[e]) - (p["path"][1] + 1)
            else:
                cs_p = len(self.edges[e]) - p["path"][1]

            if self.check_for_stopcodon(cs_p, e, cStop_codons, cMet, max_length):
                p["after_stopcodon"] = True
                start_codon_pos_filtered.append(p)
            else:
                p["after_stopcodon"] = False
                if not only_longest:
                    start_codon_pos_filtered.append(p)


        return start_codon_pos_filtered

    def find_stop_codons(self, f_p, f_edge, max_length):
        stop_codon_pos = self.find_ends(f_p - 3 + 1, f_edge, Stop_codons, Stop_codons, max_length)
        return stop_codon_pos


    def find_ends_positions(self, s_p, f_p, inner_path, startcodon_dist, only_longest):
        s_edge = inner_path[0]
        f_edge = inner_path[-1]

        if len(startcodon_dist) == 0:
            max_length = 3000
        else:
            max_length = 3*max(startcodon_dist) + 300

        start_codon_pos = self.find_start_codons(s_p, s_edge, max_length, only_longest)

        stop_codon_pos = self.find_stop_codons(f_p, f_edge, 3000)
        return start_codon_pos, stop_codon_pos