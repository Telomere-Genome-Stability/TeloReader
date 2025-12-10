# TeloReader

@author: Clotilde Garrido, Sorbonne Université - CNRS

TeloReader detects all telomere sequences in individual long reads and genome assemblies. It uses k-mer scoring in sliding windows based on similarity to a canonical telomeric motif.  
TeloReader identifies sequences ≥16 bp, either at chromosome ends (terminal) or internal (ITS), in a fasta file.

# Dependencies

  - python=3.14.2
  - numpy=2.3.5
  - pandas=2.3.3
  - seqkit=2.12.0

# Installation

### Clone the repository
```bash
git clone https://github.com/Telomere-Genome-Stability/TeloReader.git
cd TeloReader
```

### Create and activate the environment
```bash
conda create -n teloreader python=3.10
conda activate teloreader
```

### Install Python dependencies
```bash
conda install -c conda-forge numpy pandas
```

### Install SeqKit
```bash
conda install -c bioconda seqkit
```
# Usage
```bash
python Teloreader.py <strain> <fastafile> <motifs>
```
# Arguments

| Argument                   | Description                                                                             | Default value   |
|----------------------------|-----------------------------------------------------------------------------------------|-----------------|
| strain                     | Strain name. (sample label)                                                             | required        |
| fastafile                  | FASTA file used to search for telomeric sequences.                                      | required        |
| telomeric motif            | Telomeric motif. (options: Yeast / TTAGGG / TTTAGGG / TTAGG / TTGGG / TTTGGG / TTTGGGG) | required        |
| `-t`, `--threads`          | Number of threads.                                                                      | 4               |
| `-ms`, `--max_split_seq`   | Maximal number of sequences per split FASTA chunk.                                      | 1000            |
| `-w`, `--size_window`      | Sliding window size for mean score calculation. (must be < max_size_gap)                | 15              |
| `-g`, `--max_size_gap`     | Maximum size of a non-telomeric gap allowed inside a telomere.                          | 20              |
| `-mw`, `--min_mean_window` | Minimum average score per window.                                                       | motif length-1  |
| `-mint`, `--min_mean_telo` | Minimum mean score of the full telomere.                                                | motif length*0.8|
| `-minl`, `--min_len`       | Minimum telomere length. (must be > motif length)                                       | 16              |
| `-maxt`, `--max_dist_term` | Maximum distance from sequence end to be considered terminal.                           | 50              |
| `-o`, `--out_pwd`          | Output directory.                                                                       | FASTA directory |

# Outputs

### Directory created: 
**`out_teloreader/`**

### Files generated:
- **`merged_output_<strain>_<fastaname>.tsv`**: Tab-separated table listing all telomeric sequences detected in the input FASTA.

- **`merged_output_<strain>_<fastaname>.fasta`**: FASTA file containing all detected telomeric sequences.
  
# Description of output table columns
| Column     | Description                                                |
|------------|------------------------------------------------------------|
| `strain`     | Strain name (label)                                        |
| `name`       | Sequence name                                              |
| `N`          | Total number of telomeric sequences found in this sequence |
| `type`       | `C` for C-rich motif, `G` for G-rich                           |
| `len`        | Telomere length                                            |
| `start`      | Start position of the telomeric sequence                   |
| `end`        | End position of the telomeric sequence                     |
| `Loc`        | `term` for terminal telomere, `intern` for ITS                 |
| `Score_Kmer` | Mean k-mer score (motif length = perfect score)            |
| `reads_len`  | Total length of the parent sequence                        |
