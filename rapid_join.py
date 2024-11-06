#!/usr/bin/env python3

import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", required = True)
parser.add_argument("-t", "--tpf", required = True)
parser.add_argument("-c", "--csv")
parser.add_argument("-o", "--out", required = True)
args = parser.parse_args()

def read_chromosome_file(f):
    ids = []
    for l in open(f,"r"):
        items = l.strip().split(",")
        ids.extend(items)
    return ids

# get sequences from the fasta
seqs = SeqIO.index(args.fasta,"fasta")

# get tpf_lines and build new seqs
records = []
last_id = ""
seq = ""
for line in open(args.tpf,"r"):
    s = Seq("")
    columns = line.split()
    
    if columns[0] == "GAP":
        s = Seq("N"*int(columns[2]))
    else:
        (_,idline,scaffold_id,orientation)=columns
        (seq_id,start,end)=re.split("[:-]",idline)

        s = seqs[seq_id][int(start)-1 : int(end)].seq
        if orientation == "MINUS" :
            s = s.reverse_complement()
        if last_id != scaffold_id:
            if last_id != "":
                records.append(SeqRecord(Seq(seq),id=last_id,description=""))
            last_id = scaffold_id
            seq = ""
    seq = seq+str(s)
records.append(SeqRecord(Seq(seq),id=last_id,description=""))

records.sort(key=lambda i: len(i.seq), reverse=True)

# if there is a CSV
if args.csv:
    chromosome_ids = read_chromosome_file(args.csv)
    file_stem = ".".join(args.out.split(".")[0:-1])
    rl_to_super = {}
    super_to_order = {}
    super_unlocs = {}
    counter = 1
    for record in records:
        if record.id in chromosome_ids:
            new_id = ""
            if "unloc" in record.id : # unloc fragment
                ids = record.id.split("_")
                parent_id = rl_to_super["_".join(ids[0:-2])]
                new_id = parent_id + "_" + ids[-2] + "_" + ids[-1]
                super_unlocs[parent_id].append(new_id)
                super_to_order[new_id]=super_to_order[parent_id]
            elif "RL_" in record.id : # autosome
                new_id = f"SUPER_{counter}"
                rl_to_super[record.id]=new_id
                super_unlocs[new_id]=[]
                super_to_order[new_id]=counter
                counter+=1
            else: # sex chromosome
                new_id = f"SUPER_{record.id}"
                rl_to_super[record.id]=new_id
                super_unlocs[new_id]=[]
                super_to_order[new_id] = record.id
            record.id = new_id

    super_to_rl = {v: k for k, v in rl_to_super.items()}

    inter_csv = open(file_stem+".inter.csv","w")
    chromosome_csv = open(file_stem+".chromosome.list.csv","w")

    for k,v in super_to_order.items():
        flag = "yes"
        old_id = ""
        if "unloc" in k:
            flag = "no"
            ids = k.split("_")
            rl_id = super_to_rl["_".join(ids[0:-2])]
            old_id = rl_id + "_" + ids[-2] + "_" + ids[-1]
        else:
            old_id = super_to_rl[k]

        chromosome_csv.write(f"{k},{v},{flag}\n")
        inter_csv.write(f"{old_id},{v},{flag}\n")

SeqIO.write(records, args.out, "fasta")
