# ORFs search in assembly graphs

Pipeline for generating potential gene sequences, ORFs (Open Reading Frames), from assembly graphs.
Our pipeline incorporates the power of two graph alignment tools (PathRacer and SPAligner) and uses their output as initial anchors to search for full gene sequences in assembly graphs.

## Installation

The main pipeline is written in Python3 and uses several libraries described below.
Used libraries and tools: 
- python3:
    - [biopython](https://biopython.org/wiki/Download)
    - [pyyaml](https://pyyaml.org/wiki/PyYAMLDocumentation)
- [Edlib](https://pypi.org/project/edlib/)
- [HMMer](http://hmmer.org) 
- [Mummer4](https://github.com/mummer4/mummer/releases)
- [PathRacer](http://cab.spbu.ru/software/pathracer/)
- [SPAligner](http://cab.spbu.ru/software/spaligner/)

Paths to local versions of HMMer and Mummer4 can be set in config.yaml (or leave it empty if they were installed globaly).
Execution files of PathRacer and SPAligner must be in `aligners/` folder.

## Running

Search for potential Cry and Vip proteins in assembly graph of *Brevibacillus laterosporus* strain MG64(SRR8316749):
    
    run_ORFs_search.py  -m test_data/pfamA.of.interest_pfam.hmm  # list of HMMs in HMMer format that represent domains for PathRacer input
                        -s test_data/toxins.fasta                # list of known Cry and Vip sequences (either -m or -s has to be set)
                        -g test_data/graph.gfa                   # path to assembly graph
                        -k 55                                    # assembly graph k-mer size
                        -c test_data/contigs.fasta               # contig sequences (optional)
                        -t 16                                    # number of threads to use
                        -o test                                  # output folder


Try `run_ORFs_search.py -h` for more options. Test data for this example can be downloaded from [figshare](https://figshare.com/s/28de3bac33d6f0156998).

## Output

Results can be found in user defined folder `-o test`:
    
    test/pfamA.of.interest_pfam/            PathRacer run results
    test/toxins/                            SPAligner run results
    test/orfs_raw.fasta                     Full list of potential ORFs that were found in assembly graph
    test/orfs_final.fasta                   List of potential ORFs after initial filtering
    test/orfs_final_clustered.fasta         Final list of potential ORFs clustered by their identity


## Contacts

For any questions or suggestions please do not hesitate to contact Tatiana Dvorkina <tanunia@gmail.com>.