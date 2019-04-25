from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def load_sequences(seq_file, fmt = "map"):
    if fmt == "map":
        res = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))
    else:
        res = list(SeqIO.parse(seq_file, "fasta"))
    return res

def make_record(record_ini, seq, suf):
    return SeqRecord(seq, id=record_ini.id + "/" + suf, name="protein_"+record_ini.name, description=record_ini.description)

def load_seqedges(edges_file):
    seq_nuc = SeqIO.to_dict(SeqIO.parse(edges_file, "fasta"))
    seq = {}
    seq_nuc_new = {}
    for s in seq_nuc:
        seq[s + "/0"] = make_record(seq_nuc[s], seq_nuc[s].seq.translate(), "0")
        seq[s + "/1"] = make_record(seq_nuc[s], seq_nuc[s].seq[1:].translate(), "1")
        seq[s + "/2"] = make_record(seq_nuc[s], seq_nuc[s].seq[2:].translate(), "2")
        seq[s + "/0"] = make_record(seq_nuc[s], seq_nuc[s].seq.translate(), "0")
        seq[s + "/1"] = make_record(seq_nuc[s], seq_nuc[s].seq[1:].translate(), "1")
        seq[s + "/2"] = make_record(seq_nuc[s], seq_nuc[s].seq[2:].translate(), "2")
    return seq_nuc, seq

def convert_edge(p, edges):
    p_id = p + "+"
    if p_id not in edges:
        p_id = str(int(p) - 1) + "-"
    return p_id            

def convert_pos(s_pos, f_pos, shift, path, edges, K):
    # convert to nucs
    res_s = shift + s_pos * 3
    res_f = shift + f_pos * 3 + 2

    #convert to edges    
    for p_id in path[:-1]:
        res_f -= len(edges[p_id]) - K

    if res_f > 0 and res_s < len(edges[path[0]]):
        new_path = path
    else:
        new_path = [] 
        
    return res_s, res_f, new_path

def merge(edges, K, seq_name):
    path = seq_name.split("/")[0].split("_")
    offset = int(seq_name.split("/")[1])
    res = ""
    for p in path:
        res += edges[convert_edge(p, edges)][:-K]
    res += edges[convert_edge(path[-1], edges)][-K:]
    res = res[offset:]
    return SeqRecord(res, id=seq_name)

def hmmer_results_parser(hmm_filename, edges, len_th, evalue_th, K, sequences = []):
    sequences_hits = {}
    hmm_hits = {}
    hits_lst = []
    with open (hmm_filename,'rU') as handle: 
        for record in SearchIO.parse(handle, 'hmmscan3-domtab'):
            hmm_name = record.id
            hmm_len = record.seq_len
            for h in record.hits:
                seq_name = h.id
                hit_evalue = h.evalue
                for f in h.fragments:
                    hit_len = f.hit_end - f.hit_start
                    seq_s, seq_e = f.query_start, f.query_end - 1
                    if hit_len > len_th*hmm_len and hit_evalue < evalue_th:
                        if len(sequences) > 0:
                            cur_seq = sequences[seq_name]
                        else:
                            cur_seq = merge(edges, K, seq_name)
                        hits_lst.append({"seq_name": seq_name,"hmm_name": hmm_name, "start": seq_s, "end": seq_e, \
                                         "seq": cur_seq.seq[seq_s:seq_e], "e-val": hit_evalue })
                        if seq_name not in sequences_hits:
                            sequences_hits[seq_name] = []
                        if hmm_name not in hmm_hits:
                            hmm_hits[hmm_name] = []
                        sequences_hits[seq_name].append({"hmm_name": hmm_name, "start": seq_s, "end": seq_e})
                        hmm_hits[hmm_name].append({"seq_name": seq_name, "start": seq_s, "end": seq_e, "len": len(cur_seq.seq), "hmm_len": hmm_len})
    return hits_lst, sequences_hits, hmm_hits

def load_true_crygenes(genes_file, hmm_filename, edges, K):
    genes = load_sequences(genes_file)
    hits_lst, genes_hits, hmm_hits = hmmer_results_parser(hmm_filename, edges, 0.9, 0.000000001, K, genes)
    return hmm_hits

def load_pathracer_mapping(domtblfile, edges, len_th, evalue_th, K):
    hits_lst, _, _ = hmmer_results_parser(domtblfile, edges, len_th, evalue_th, K)
    pathracer_lst = []
    for hit in hits_lst:
        path, shift = hit["seq_name"].split("/")[0].split("_"), int(hit["seq_name"].split("/")[1])
        for i in range(len(path)):
            path[i] = convert_edge(path[i], edges)
        seq_s, seq_e, path = convert_pos(hit["start"], hit["end"], shift, path, edges, K) 
        if len(path) > 0:
            pathracer_lst.append({"seq_name": hit["seq_name"], "name": hit["hmm_name"], "start": seq_s, \
                                    "end": seq_e, "path": path})
    return pathracer_lst

def load_spaligner_mapping(mappings_fasta):
    res = []
    seqs = load_sequences(mappings_fasta, "list")
    for s in seqs:
        lst = s.description.split("|")
        ind = 0
        while not lst[ind].startswith("Edges="):
            ind += 1
        path, start_g, end_g, start_s = lst[ind][len("Edges="):].split("_")[:-1], \
                                        int(lst[ind+1][len("start_g="):]), int(lst[ind+2][len("end_g="):]) - 1, int(lst[ind+3][len("start_s="):])
        for i in range(len(path)):
            if not path[i].endswith("-") and not path[i].endswith("+"):
                path[i] = path[i] + "+" 
        res.append({"name": "_".join(lst[:ind]), "start": start_g, "end": end_g, "path": path, "d":start_s})
    return res
