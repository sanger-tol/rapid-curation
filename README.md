# Rapid curation 

## Support 
Slack Channel for curation assembly help:
```
https://join.slack.com/t/assemblycuration/shared_invite/zt-yezxfd4w-P0xJdV1TJg47OaQKJaqlhA
```

PretextView tutorial
```
https://youtu.be/3IL2Q4f3k3I
```

PretextView can be obtained here
```
https://github.com/wtsi-hpag/PretextView/releases
```

# Software guide
This section will cover the software side of things and how it works.
Tested on lsf and slurm

## Set up

First directory variables must be set up, this will tell the singularity image
where the working directory and the data are located.
```
WORKDIR=/Where/you/will/be/working
DESTDIR=/Where/the/output/will/go
```

For Example:
```
WORKDIR=/user/grit/curation/organism/
DESTDIR=/user/grit/curation/organism/output
```

Next, we must group these variables together and bind them to a location inside
the singularity image.
```
export SINGULARITY_BIND="
/nfs:/nfs,\ <-- This will need changing depending on your working environment.
/lustre:/lustre,\ <-- This will need chaniging depending on where data is stored.
$WORKDIR:/data,\ <-- This does not need changing
$DESTDIR:/output,\ <-- This does not need changing
"
``` 

$SINGULARITY_BIND mounts the specified paths to paths inside the singularity image.
Files on paths not noted here will not be be visible and cause the pipeline to fail.
Only the text before the ':' should be changed to the base directory e.g. `/User, /nfs, /lustre`


For Example:
```
export SINGULARITY_BIND="
/user:/nfs,\
/user:/lustre,\
$WORKDIR:/data,\
$DESTDIR:/output,\
"
```

## Running the individual steps
### Command Flags
EXPLAIN

EXPLAIN bsub is LSF specific, if you use SLURM you need to change your job command

### Extra Notes
*_done files found in the output directories are empty files which aid in the control of
the underlying snakemake pipeline.

*.conf files in the same folders detail the commands used by the pipeline inside of the Singularity environment.

### HiC - Generates the HiGlass files and Pretext Map
<details open>
<summary> Open </summary>

EXPLANATION

The output of this program includes:
- ref.fa		<-- Original fasta file
- ref.fa.amb
- ref.fa.ann
- ref.fa.bwt		<-- Burrows Wheeler Transform file
- ref.fa.fai		<-- Fasta index file
- ref.fa.pac
- ref.fa.sa
- MORE

Example, used in GRIT:
```
echo '/software/singularity-v3.6.4/bin/singularity run /lustre/scratch123/tol/teams/grit/yy5/sigs/runHiC.sif -q 0 -s ilWatBina1_dr'|bsub -J test_sig -q basement -o test_sig.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```
</details>

### Coverage Track - Coverage plot used by HiGlass and Pretext
<details open>
<summary> Open </summary>

EXPLANATION

The output of this program includes:
- coverage.bed			<--
- geval.bed			<--
- {sample}.bw			<-- HiGlass and Pretext
- {input Pacbio fasta}.bam	<-- Input fasta converted to bam
- merged.bam			<-- Multiple bams are merged
- merged_sort.bam		<-- Merged bam is then sorted
- my.genome			<--
- pre.bam			<--
- ref.mmi			<-- Minimap Index File
- sort_tmp/			<--

Example, used in GRIT:
```
echo "/software/singularity-v3.6.4/bin/singularity run /nfs/team135/yy5/sif/runCoverage.sif -t ilWatBina1_dr" | bsub -J coverage -q basement -o coverage.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```
</details>

### Repeat Density Track - Needed for Pretext
<details open>
<summary> Open </summary>

EXPLANATION

The output of this program includes:
- bin.bed		<--
- density_nodot.bed	<--
- ref.fa.density.bed	<-- 
- ref.fa.genome		<--
- ref.fa.intersect.fa	<--
- ref.fa.repeat.bed	<--
- ref.fa.stage1		<--
- ref.fa.wm		<--
- sorted.genome		<--
- sorted_intersect.bed	<--
- {sample}_repeat_density.bw <-- HiGlass and Pretext

Example, used in GRIT:
```
echo "/software/singularity-v3.6.4/bin/singularity run /nfs/team135/yy5/sif/runRepeat.sif -t ilWatBina1_dr" | bsub -J repeat -q basement -o repeat.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```
</details>

### Gap Track - Needed for HiGlass and Pretext
<details open>
<summary> Open </summary>

The Gap Track utilises SeqTK to identity strings of N's wich indicate a gap in the sequence. This can then be used to inform curators on possible locations they can break scaffolds.
The output of this program is:

- {sample}_gap.bed	<-- HiGlass
- {sample}_gap.bedgraph <-- Pretext

seqtk can be found on GitHub [here](https://github.com/lh3/seqtk).

Example, used in GRIT:
```
echo "/software/singularity-v3.6.4/bin/singularity run /nfs/team135/yy5/sif/runRepeat.sif -t ilWatBina1_dr" | bsub -J gap -q basement -o gap.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```
</details>

### Telomere Track - Needed for HiGlass and Pretext
<details open>
<summary> Open </summary>

The Telomere Track uses the find_telomere program created by the VGP and modified to accept a custom telomere motif to search a given genome.
The output of this program include:

- {sample}_telomere.bed      <-- HiGlass
- {sample}_telomere.bedgraph <-- Pretext
- ref.telomere		     <-- 
- ref.windows		     <-- 

The original VGP's find_telomere files are found [here](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere).	

Example, used in GRIT:
```
echo "/software/singularity-v3.6.4/bin/singularity run /nfs/team135/yy5/sif/runTelo.sif -t ilWatBina1 -s ilWatBina1_dr" | bsub -J telo -q basement -o telo.o -n 10 -M50000 -R'select[mem>50000] rusage[mem=50000] span[hosts=1]'
```
</details>


# Curation guide
This section will cover the curation tools and how they work.

## Suggested workflow

Split assembly fasta to create a TPF (rapid_split.pl).
Manipulate the assembly in PretextView app.  Action any necessary breaks in TPF *contigs* in order to mirror breaks in PretextView.  This involves editing coordinates and introducing new gaps in the TPF.
Once the map has been rearranged, paint the chromosomes.
Export the map to AGP.  (Some input scaffolds won't feature in the output AGP due to PretextView resolution.  These will be recovered in the next step)
Run rapid_pretext2tpf.py, pointing it at the edited TPFyes (where all necessary gaps have been added) and the Pretext AGP.
This will output a TPF that mirrors the chromosomes that have been built in PretextView, and which will have the same number of basepairs as the original TPF.
Turn the new TPF back into a fasta file (rapid_join.pl)
A count of the number of breaks/joins and haplotypic duplicate removals can be produced using rapid_stats.py


Scripts and documentation for Rapid curation

## rapid_split.pl

This script takes a fasta file and produces a TPF on contig level, i.e. it splits at all Ns
Usage:
perl split.pl -fa <fasta>
              -printfa # if you want to print out the split fasta
              -h/help  # this message
	      
##  rapid_pretext2tpf.py	      

usage: rapid_pretext2tpf.py [-h] tpf agp

Designed to take pretext generated AGP and fit your assembly TPF to it.

positional arguments:
  tpf         assembly TPF with gaps as needed to allow rearrangement to match
              the edited PretextView map.
  agp         Pretext agp
  
  
## rapid_join.pl

This script takes an original assembly fasta file, a one-line per chromosome pre-csv file and a TPF file generated with split.pl and creates the assembly from the TPF file
Usage:
perl join.pl -fa <fasta>
             -tpf <tpf>
             -csv <pre-csv>
             -out <outfile_fasta_prefix> # unless specified the output will be written to <fasta>.curated.fasta
             -hap # optional use only if generating haplotigs fasta
             -h/help # this message
	     
## rapid_stats.py

Calculates number of breaks, joins and haplotypic duplicate removals.

positional arguments:
  tpf1           original tpf
  tpf2           curated tpf
  tpf3           haplotigs tpf

optional arguments:
  -h, --help     show this help message and exit
  --rapid RAPID  tpf is generated by rapid_split.pl
usage: rapid_stats.py [-h] [--rapid RAPID] tpf1 tpf2 tpf3



	      

## Getting started

To make it easy for you to get started with GitLab, here's a list of recommended next steps.

Already a pro? Just edit this README.md and make it your own. Want to make it easy? [Use the template at the bottom](#editing-this-readme)!

## Add your files

- [ ] [Create](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://gitlab.com/wtsi-grit/rapid-curation.git
git branch -M main
git push -uf origin main
```

## Integrate with your tools

- [ ] [Set up project integrations](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/user/project/integrations/)

## Collaborate with your team

- [ ] [Invite team members and collaborators](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Automatically merge when pipeline succeeds](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing(SAST)](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://docs.gitlab.com/ee/user/application_security/sast/)

***

# Editing this README

When you're ready to make this README your own, just edit this file and use the handy template below (or feel free to structure it however you want - this is just a starting point!).  Thank you to [makeareadme.com](https://gitlab.com/-/experiment/new_project_readme_content:98d9f6afce7ffb152eed371c794ed7b5?https://www.makeareadme.com/) for this template.

## Suggestions for a good README
Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.

## Name
Choose a self-explaining name for your project.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.

