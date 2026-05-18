# Early embryonic transcription patterns of the marker gene *Vasa* in the Mediterranean mussel confirm divergent routes to germline specification in bivalves

In this repository you will find data and codes used to perform the analyses for the research paper:

> **[Nicolini F](https://github.com/filonico), Nuzhdin S, [Ghiselli F](https://github.com/fghiselli), [Luchetti A](https://github.com/andluche), Milani L**. *Early embryonic transcription patterns of the marker gene* Vasa *in the Mediterranean mussel confirm divergent routes to germline specification in bivalves.*
>
> **Astract.** TBA.

Visit our research group website, **[EVO·COM](https://sites.google.com/view/evo-com-unibo)**!

## What you'll find here
In this repository there are all the code/scripts, input files, and intermediate results that have been used and generated for the paper.

Note that the pipeline can be run either through Snakemake (see [`Snakefile`](Snakefile)) or by manually executing commands from the command-line (see [`pipeline.sh`](pipeline.sh)).

Here's the description of each folder content/file:

* [`00_input/`](00_input) contains starting input files used throughout the analysis;
    * [`00_input/OG0007197_vasa.fa`](00_input/OG0007197_vasa.fa) is the fasta file of bivalve Vasa (Ddx3) orthogroup from [Nicolini et al., 2025](https://github.com/filonico/bivalvia_SRGs/);
    * [`00_input/OG0007978_ddx3_OUT.fa`](00_input/OG0007978_ddx3_OUT.fa) is the fasta file of bivalve Ddx4 orthogroup from [Nicolini et al., 2025](https://github.com/filonico/bivalvia_SRGs/);
    * [`00_input/PF00270.alignment.full.hmm`](00_input/PF00270.alignment.full.hmm) is the Pfam hmm profile of the DEAD/DEAH box helicase;
    * [`00_input/PF00270.alignment.full.stk.gz`](00_input/PF00270.alignment.full.stk.gz) is the Pfam stockholm alignment of the DEAD/DEAH box helicase (same as before, but different format);
    * [`00_input/SRA_metadata.tsv`](00_input/SRA_metadata.tsv) is the metadata table of RNA-seq libraries from [Miglioli et al., 2024](https://doi.org/10.1242/dev.202256);
    * [`00_input/SRA_toDownload.ls`](00_input/SRA_toDownload.ls) is the SRA list of RNA-seq libraries from [Miglioli et al., 2024](https://doi.org/10.1242/dev.202256);
    * [`00_input/contaminants2trimm.fa`](00_input/contaminants2trimm.fa) is a fasta file containing most commons Illumina adapters to trim from fastq files;
    * [`00_input/edesign_object.tsv`](00_input/edesign_object.tsv) is the eDesign object used as input to perform the maSigPro time-series analysis;
    * [`00_input/edesign_object_24hpf.tsv`](00_input/edesign_object_24hpf.tsv) is the eDesign object used as input to perform the maSigPro time-series analysis, but restricted to the first 24h samples;
    * [`00_input/mgal_genome.tar.gz`](00_input/mgal_genome.tar.gz) is the compressed directory containing the *Mytilus galloprovincialis* genome to which reads were mapped (NCBI ID: [GCA_900618805.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_900618805.1/); processed as in [Nicolini et al., 2025](https://github.com/filonico/bivalvia_SRGs/));
    * [`00_input/mgal_proteome_GOannotation.tsv`](00_input/mgal_proteome_GOannotation.tsv) is the GO term annotation of the *M. galloprovincialis* predicted proteome, as returned by the [OMA browser](https://omabrowser.org/oma/functions/);
    * [`00_input/vasa_piwi_nanos_spPHI.faa`](00_input/vasa_piwi_nanos_spPHI.faa) is the fasta file with sequences of some germline markers;
    * [`00_input/vasa_ref.faa`](00_input/vasa_ref.faa) is the fasta file with *Vasa* sequences from reference species;
* [`Snakefile`](Snakefile) is the snakefile to run the pipeline;
* [`config.yaml`](config.yaml) is the snakemake configuration file;
* [`conda_envs/`](conda_envs/) contains the YAML files of the conda environments used in the analyses;
* [`intermediate_results/`](intermediate_results/) contains some of the intermediate results obtained throughout the analyses;
* [`scripts/`](scripts/) contains all the scripts used in the analyses, each with extensive code comments;
* [`pipeline.sh`](pipeline.sh) is the pipeline documenting the analyses, containing extensive code comments; mind that every command is supposed to be run from the current directory;
* [`protocols/`](protocol/) contains the wet-lab protocols used for embryo rearing, fixation and HCR.