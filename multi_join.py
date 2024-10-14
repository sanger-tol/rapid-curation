#!/bin/env python3

import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", required = True)
parser.add_argument("-t", "--tpf")
parser.add_argument("-c", "--csv")
parser.add_argument("-o", "--out", required = True)
parser.add_argument("-t2", "--tpf2")
parser.add_argument("-c2", "--csv2")
parser.add_argument("-d", "--hap")
args = parser.parse_args()

out1 = f'{args.out}.hap1.1'
out2 = f'{args.out}.hap2.1'
outt = f'{args.out}.1'

def read_chromosome_file(f):
    ids = []
    for l in open(f,"r"):
        items = l.strip().split(",")
        ids.extend(items)
    return ids

# get sequences from the fasta
seqs = SeqIO.index(args.fasta,"fasta")

# get tpf_lines and build new seqs

def join_full(tpf, csv, out):

    records = []
    last_id = ""
    seq = ""
    for line in open(tpf,"r"):
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

    chromosome_ids = read_chromosome_file(csv)
    file_stem = ".".join(out.split(".")[0:-1])
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

    inter_csv = open(out+".inter.csv","w")
    chromosome_csv = open(out+".primary.chromosome.list.csv","w")

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

    SeqIO.write(records, out+".primary.curated.fa", "fasta")

def join_short(tpf,out):

    records = []
    last_id = ""
    seq = ""
    for line in open(tpf,"r"):
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

    SeqIO.write(records, out+".primary.curated.fa", "fasta")

def join_hap(tpf,out):

    records = []
    last_id = ""
    seq = ""
    for line in open(tpf,"r"):
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

    SeqIO.write(records, out+".1.additional_haplotigs.curated.fa", "fasta")

if args.csv2:
    join_full(args.tpf, args.csv, out1)
    join_full(args.tpf2, args.csv2, out2)
    with open(out1+".all_haplotigs.curated.fa", 'w') as fp:
        pass
    with open(out2+".all_haplotigs.curated.fa", 'w') as fp:
        pass
elif args.tpf2:
    join_full(args.tpf, args.csv, out1)
    join_short(args.tpf2, out2)
    with open(out1+".all_haplotigs.curated.fa", 'w') as fp:
        pass
    with open(out2+".all_haplotigs.curated.fa", 'w') as fp:
        pass
    with open(out2+".primary.chromosome.list.csv", 'w') as fp:
        pass
else:
    join_full(args.tpf, args.csv, outt)

if args.hap:
    join_hap(args.hap, args.out)