# OrganelleCleaner

OrganelleCleaner identifies likely organelle-derived contigs in a nuclear genome assembly and writes a cleaned nuclear assembly after those contigs are removed. It uses assembly-graph structure, contig-level sequence features, and optional internal BLAST searches against plastid/chloroplast and mitochondrial reference sequences.

This tool is intended for assemblies in which plastid- and mitochondrial-derived contigs co-assemble with the nuclear genome and need to be screened before downstream analyses such as curation, scaffolding, assembly statistics, or release preparation.

## Table of Contents

- [Overview](#overview)
- [When to use each mode](#when-to-use-each-mode)
- [Inputs](#inputs)
- [Recommended organelle reference generation](#recommended-organelle-reference-generation)
- [How it works](#how-it-works)
- [Installation](#installation)
- [Usage examples](#usage-examples)
- [Outputs](#outputs)
- [BLAST support policy](#blast-support-policy)
- [Optional expert settings](#optional-expert-settings)
- [Interpretation and caveats](#interpretation-and-caveats)

## Overview

Organelle Cleaner is designed to answer a practical question: which contigs in a nuclear assembly are likely derived from plastid/chloroplast or mitochondrial genomes rather than from the nuclear genome itself?

Why this matters:

- organelle contigs can inflate nuclear assembly size
- they can distort GC and coverage summaries
- they can complicate scaffolding and curation
- they often need to be removed before reporting a final nuclear assembly

At a high level, the tool:

- reads an assembly graph in GFA format
- extracts graph and sequence-based evidence for each contig
- optionally runs internal `makeblastdb` and `blastn` against organelle reference FASTA files
- scores evidence for organelle origin
- writes a cleaned assembly FASTA and per-contig reports

Typical inputs:

- required: assembly graph (`.gfa`)
- optional for BLAST-supported modes:
  - plastid / chloroplast reference FASTA
  - mitochondrial reference FASTA

Main outputs:

- `cleaned_assembly.fa`
- `organelle_contigs.txt`
- `nuclear_contigs.txt`
- `report.tsv`

## When to Use Each Mode

### `graph-only`

Uses:

- assembly-graph structure
- contig length, GC content, and coverage-related signals

Choose this when:

- you do not have plastid or mitochondrial reference sequences
- you want a lightweight first-pass screen
- you want to avoid BLAST entirely

Requirements:

- GFA required
- assembly FASTA optional
- no BLAST tools required

Pros:

- simple and fast
- no organelle references needed

Cons:

- usually less sensitive than hybrid mode

### `hybrid`

Uses:

- graph-based evidence
- internal BLAST evidence against plastid and/or mitochondrial references

Choose this when:

- organelle reference sequences are available
- you want the most practical default mode for real datasets

Requirements:

- GFA required
- at least one of `--plastid-fasta` or `--mit-fasta`

Pros:

- usually the best balance of sensitivity and specificity
- can use plastid only, mitochondrial only, or both

Cons:

- requires BLAST+ and suitable reference FASTA input

### `blast-only`

Uses:

- internal BLAST evidence only

Choose this when:

- you want to evaluate BLAST evidence without graph-aware refinement
- you want a simpler BLAST-based comparison run

Requirements:

- GFA required
- assembly FASTA required
- at least one of `--plastid-fasta` or `--mit-fasta`

Pros:

- useful for benchmarking or ablation-style comparisons

Cons:

- does not use graph context to refine calls
- generally less robust than hybrid mode for practical filtering

## Inputs

### Required input

- assembly graph in GFA format

### Optional but often useful

- organelle reference FASTA via `--plastid-fasta` and/or `--mit-fasta`

Advanced fallback note:
If the GFA does not contain usable contig sequences, you can still provide `--assembly-fasta`. Most users do not need this because GFA sequence content can be used directly for internal BLAST subject preparation.

### Optional organelle reference FASTA inputs

- `--plastid-fasta`
- `--mit-fasta`

These FASTA files are used only for internal BLAST-based evidence. You may provide:

- plastid only
- mitochondrial only
- both plastid and mitochondrial references

Important workflow rule:

- `graph-only` does not require BLAST and does not require organelle FASTA input
- `hybrid` and `blast-only` fail only if neither `--plastid-fasta` nor `--mit-fasta` is provided
- when the GFA contains contig sequences, internal BLAST can use them directly as the assembly subject source

## Recommended Organelle Reference Generation

For best BLAST-based evidence, users are encouraged to provide high-quality organelle reference sequences, ideally assembled from the same sample or from a closely related accession.

When HiFi data are available, a practical approach is to assemble plastid and mitochondrial genomes first and then use those assemblies as Organelle Cleaner references. HiMT is recommended for this purpose when suitable:

- HiMT: https://github.com/aiPGAB/HiMT

This is a recommendation, not a requirement. Organelle Cleaner does not depend on HiMT directly.

## How It Works

At a user level, the pipeline is:

1. Parse the assembly graph.
   The tool reads the GFA and records the contigs and their graph connectivity.

2. Extract contig-level features.
   It summarizes graph structure and sequence-derived features such as contig length, GC content, and coverage.

3. Optionally run internal BLAST.
   If plastid and/or mitochondrial FASTA references are provided in `hybrid` or `blast-only` mode, Organelle Cleaner builds a BLAST nucleotide database from the available assembly contig sequences and runs `blastn` internally. It uses GFA sequence content when available and falls back to `--assembly-fasta` only when needed.

4. Classify BLAST support.
   BLAST hits are summarized per assembly contig, with overlapping or split HSPs merged on the subject contig before support coverage is measured.

5. Combine evidence.
   Depending on the selected mode, graph evidence, BLAST evidence, or both are combined to identify likely organelle-derived contigs.

6. Write cleaned outputs.
   The tool writes a cleaned nuclear assembly FASTA, contig ID lists, and a detailed report table.

## Installation

Create and activate the conda environment, then install the package in editable mode:

```bash
conda env create -f environment.yml
conda activate OrganelleCleaner
pip install -e .
```

This installs the `organelle-cleaner` command and the dependencies listed in `environment.yml`, including BLAST+.

Requirements:

- Python 3.11 in the provided conda environment
- `networkx`
- BLAST+ tools: `makeblastdb` and `blastn`

## Usage Examples

### Minimal graph-only run

```bash
organelle-cleaner assembly.gfa \
  --output-dir results
```

### Hybrid mode with plastid FASTA only

```bash
organelle-cleaner assembly.gfa \
  --plastid-fasta plastid_refs.fa \
  --mode hybrid \
  --output-dir hybrid_plastid_results \
  --threads 8
```

### Hybrid mode with plastid and mitochondrial FASTA

```bash
organelle-cleaner assembly.gfa \
  --plastid-fasta PI652213_chl.fa \
  --mit-fasta PI652213_mit.fa \
  --mode hybrid \
  --output-dir results \
  --threads 65
```

### Blast-only mode

```bash
organelle-cleaner assembly.gfa \
  --plastid-fasta plastid_refs.fa \
  --mode blast-only \
  --output-dir blast_only_results \
  --threads 8
```

Advanced fallback note:
If your GFA does not include usable contig sequences, you can still provide `--assembly-fasta` explicitly. If you need a separate candidate-only report file, `--all-candidates-name` remains available as an advanced option.

## Outputs

Main output files:

- `cleaned_assembly.fa`
  FASTA containing contigs retained as likely nuclear.

- `organelle_contigs.txt`
  One contig ID per line for contigs flagged as likely organelle-derived.

- `nuclear_contigs.txt`
  One contig ID per line for contigs retained as likely nuclear.

- `report.tsv`
  Per-contig summary table with graph features, sequence features, BLAST support summaries, scores, and final classification.

- candidate report specified by `--all-candidates-name`
  Per-contig report for all contigs flagged at any non-default confidence tier. This is written only when `--all-candidates-name` is provided, using the filename you supply.

### `report.tsv`: useful columns

The report contains many columns, but most users will focus on:

- `contig_id`
  Contig name from the assembly.

- `mode`
  The analysis mode used for scoring.

- `contig_length`, `coverage`, `gc_content`
  Basic contig-level properties.

- `blast_support_level`
  Overall BLAST support tier for that contig.

- `blast_support_sources`
  Which organelle source contributed support (`plastid`, `mitochondrial`, or both).

- `blast_selected_source`
  The stronger BLAST source selected for that contig.

- `plastid_support_level`, `mit_support_level`
  Per-source BLAST support summaries.

- `graph_score`, `blast_score`, `final_score`
  The component scores that drive the final classification.

- `confidence_tier`
  Confidence category assigned by the scoring system.

- `classification`
  Final class used in the output workflow.

- `supporting_evidence`, `decision_explanation`
  Short text summaries explaining the call.

### Internal BLAST intermediate files

When BLAST is used, intermediate files are written under:

```text
OUTPUT_DIR/blast_intermediates/
```

This directory includes:

- `assembly_blast_db*`
- `plastid_hits.tsv` when plastid BLAST is run
- `mitochondrial_hits.tsv` when mitochondrial BLAST is run

## BLAST Support Policy

The user-facing workflow now uses FASTA inputs and internal BLAST, but the downstream BLAST summarization logic remains the same in spirit:

- BLAST hits are attached to assembly contigs through `sseqid`
- subject-coordinate HSPs are merged before coverage is calculated
- support coverage is recalculated on identity-filtered hit subsets
- plastid and mitochondrial evidence are summarized separately
- the stronger source is used as the selected BLAST source

Current default support thresholds:

- moderate support: `90.0` identity and `0.8` merged subject coverage
- strong support: `98.0` identity and `0.8` merged subject coverage
- exceptional support: `99.0` identity and `0.95` merged subject coverage

Threshold-policy note:

- the moderate-support default changed in this refactor
- previous default: `95.0` identity and `0.5` merged subject coverage
- current default: `90.0` identity and `0.8` merged subject coverage

This was an intentional policy change to match the requested default minimum support rule of 90% identity over 80% of the subject contig. Strong and exceptional defaults were not changed.

## Optional Expert Settings

Most users can run the tool with only the main inputs shown above. Advanced settings are available for users who need to tune behavior for unusual assemblies.

Conceptually, these settings fall into a few groups:

### Baseline overrides

- `--nuclear-coverage`
- `--gc-baseline`

By default, Organelle Cleaner infers nuclear coverage baseline and GC baseline automatically from the input assembly. These options are overrides for cases where you already know the correct baseline and want to force it explicitly.

### Graph confidence thresholds

- `--high-confidence-threshold`
- `--medium-confidence-threshold`
- `--low-confidence-threshold`

These control how much graph-derived evidence is required before a contig is promoted into a higher confidence tier.

### Size-bin settings

- `--very-small-max-length`
- `--small-fragment-max-length`
- `--min-organelle-length`
- `--preferred-max-length`
- `--intermediate-max-length`
- `--large-contig-length`
- `--very-large-contig-length`

These define the contig length bins used by the scoring rules. They matter because very small fragments, preferred-size organelle contigs, and unusually large contigs are treated differently.

### Large-contig penalties and signal requirements

- `--isolated-length-cap`
- `--large-contig-penalty`
- `--very-large-contig-penalty`
- `--large-contig-high-min-signals`
- `--large-contig-medium-min-signals`

These settings control how cautious the tool is with large or very large contigs, which often need stronger supporting evidence to avoid false positives.

### GC and coverage sensitivity

- `--gc-deviation-threshold`
- `--coverage-multiplier-threshold`

Practical interpretation:

- `--gc-deviation-threshold` controls how different a contig's GC content must be from the inferred nuclear baseline before it counts as positive evidence.
- `--coverage-multiplier-threshold` controls how much higher a contig's coverage must be than the inferred nuclear baseline before it counts as positive evidence.

In most runs, the defaults are appropriate. These settings are best treated as expert tuning knobs rather than routine required inputs.

## Interpretation and Caveats

Organelle-derived contig detection is evidence-based, not perfect. Results should be interpreted in the context of the assembly and the available references.

Practical guidance:

- `hybrid` mode is generally recommended when organelle reference FASTA sequences are available
- `graph-only` remains useful when BLAST references are unavailable
- `blast-only` is mainly useful for comparison, benchmarking, or simpler BLAST-focused runs

Be cautious when interpreting results if:

- the assembly is highly fragmented
- organelle references are incomplete or distant from the sample
- the assembly contains transferred organelle-derived nuclear segments
- true biological structure makes graph or BLAST evidence ambiguous

As with any evidence-based screen, the report table should be reviewed if a contig appears borderline or biologically unexpected.
