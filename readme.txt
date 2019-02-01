ORFs search in assembly graphs

Prerelease version.

Used libraries: 
* python3:
   - biopython https://biopython.org/wiki/Download
   - pyyaml https://pyyaml.org/wiki/PyYAMLDocumentation  
* HMMer http://hmmer.org 
* Mummer4 https://github.com/mummer4/mummer/releases

Paths to local versions of HMMer and Mummer4 can be set in config.yaml (or leave it empty if they were installed globaly)

********

Running:

Search for hmms from MON.pfam.20171025.hmm in assembly graph of MTG000543 dataset:
./run_ORFs_search.py  -m test_data/MON.pfam.20171025.hmm -g test_data/graph.gfa -k 55 -s test_data/toxins.fasta -c test_data/contigs.fasta -t 16 -o test

Try ./run_ORFs_search.py -h for more options.

********

Main pipeline:
1) Searches for HMMs alignments in given assembly graph
2) Searches for Protein alignments in assembly graph (optional)
3) Identifies ORFs near each alignment and prints all ORFs in fasta-format
4) Filters resulting ORFs presented in contigs and in known proteins