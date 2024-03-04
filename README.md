``# Rapid curation - Curation tools guide

This section will cover the curation tools and how they work.



## Support 

Slack Channel for curation assembly help:
```
https://join.slack.com/t/assemblycuration/shared_invite/zt-1kx2ww71y-823ruaAxswgQGypgofBaOA
```

PretextView tutorial:

PretextView - Tutorial.pdf
```
https://youtu.be/3IL2Q4f3k3I
https://youtu.be/LWy6pwCQNDU
```

PretextView can be obtained here
```
https://github.com/wtsi-hpag/PretextView/releases
```

Please also see documentation:

RAPID CURATION TRAINING MANUAL.pdf



## Training data

Pretext maps and assembly files can be found in the test_data/ repository.

A training environment can be accessed (requires github account) that contains all the scripts and data, feel free to have a play!

http://gitpod.io/#https://github.com/thebgacademy/manual-curation


 
## Suggested workflow

Split assembly fasta to create a TPF (rapid_split.pl).

Manipulate the assembly in PretextView app using Edit mode (E).
Tag unlocalised scaffolds, haplotigs, contaminats, haplotypes etc, using meta-data tag mode (M) in PretextView.
Unlocs (linked to chromosome but no specific location found) should be placed at the start/and/or/end 
of a chromosome to be included when the chromosomes are painted in the subsequent step.

Once the map has been rearranged and meta-data tags added paint the chromosomes using Scaffold Edit mode (S).

Export the map to AGP. (Some input scaffolds won't feature in the output AGP due to PretextView resolution.  
These will be recovered in the next step).

Run pretext-to-tpf , pointing it at the TPF and the Pretext AGP. 
This will output TPF(s) that mirror the chromosomes that have been built in PretextView, and which will have the same number 
of basepairs as the original TPF.The script will produce stats (joins/breaks/hap-dup removals).

Turn the new TPF back into a fasta file using rapid_join.pl.

Produce a new Pretext map from the curated fasta file to check it.



## Scripts

rapid_split.pl
rapid_join.pl
pretext-to-tpf.py
telo_finder.py

**Scripts and documentation for Rapid curation:**

## rapid_split.pl

This script takes a fasta file and produces a TPF on contig level, i.e. it splits at all Ns

Usage:

```
perl rapid_split.pl
              -fa <fasta>

              -printfa # if you want to print out the split fasta

              -h/help  # this message
	      
```


##  pretext-to-tpf.py	

git clone https://github.com/sanger-tol/agp-tpf-utils.git

This script takes the AGP output from PretextView and the TPF from rapid_split to generate new assembly tpfs

Usage: 

```
pretext-to-tpf [OPTIONS]

Options:
  -a, --assembly PATH             Assembly file from before curation, which is
                                  usually a TPF.  [required]
  -p, --pretext PATH              Assembly file from Pretext, which is usually
                                  an AGP.  [required]
  -o, --output FILE               Output file, usually a TPF. If not given,
                                  prints to STDOUT in 'STR' format.
  -c, --autosome-prefix TEXT      Prefix for naming autosomal chromosomes.
                                  [default: RL_]
  -f, --clobber / --no-clobber    Overwrite an existing output file.
                                  [default: no-clobber]
  -l, --log-level [DEBUG|INFO|WARNING|ERROR|CRITICAL]
                                  Diagnostic messages to show.  [default:
                                  INFO]
  -w, --write-log / -W, --no-write-log
                                  Write messages into a '.log' file alongside
                                  the output file  [default: no-write-log]
  --help                          Show this message and exit.
```

See https://github.com/sanger-tol/agp-tpf-utils# for complete documentation



## rapid_join.pl

This script takes an original assembly fasta file, a one-line per chromosome pre-csv file and the TPF file(s) generated from pretext-to-tpf and creates the finalised assembly fasta from the TPF file(s)

Usage:

```
perl rapid_join.pl
             -fa <fasta>
             -tpf <tpf>
             -csv <pre-csv>
             -out <outfile_fasta_prefix>
             -hap # optional use only if generating haplotigs fasta

```




## Other scripts

Telomere identification script (telo_finder.py):

The telo_finder script aims to detect likely telomere motifs in assemblies.
It assumes the assemblies are high quality (hi-fi or similar) and that there are a number of chromosomes (it won't return a result with just one chromosome for example).
It uses a hard-coded set of non-redundant telomere motifs with which it attempts to distinguish true telomere sequences from background noise.


Usage: 
```
telo_finder.py [-h] [--size size] [--klo klo] [--khi khi] [--ends ends] fasta

Finds most likely telomere motif in Hi-fi or equivalent quality assembly where
telomeres are expected to occur at the ends of multiple scaffolds

positional arguments:
  fasta        assembly fasta

optional arguments:
  -h, --help   show this help message and exit
  --size size  top and tail scf bp (default: 200)
  --klo klo    min kmer (default: 4)
  --khi khi    max kmer (default: 15)
  --ends ends  ends to scan (default: 1000)
  
```

