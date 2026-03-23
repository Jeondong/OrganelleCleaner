# OrganelleCleaner

Detect organelle-derived contigs from GFA assemblies using graph-only, BLAST-only, or hybrid mode.

Graph-only remains the default CLI mode for backward compatibility, so a plain invocation works without BLAST inputs. Hybrid is the recommended mode for practical use and benchmarking because it combines BLAST recall with graph-based confidence refinement and large-contig suppression. Hybrid requires at least one plastid/chloroplast or mitochondrial BLAST TSV input.

Blast-only mode also requires at least one BLAST TSV input and uses only BLAST-derived scoring, but the current CLI still uses the GFA as the assembly contig inventory/input source.

## Installation

```bash
conda env create -f environment.yml
conda activate organelle-cleaner
pip install -e .
```

The current codebase requires Python 3.10+ syntax. BLAST itself is not required for installation or for the default `graph-only` mode. It is only needed externally if you want to generate BLAST TSV inputs for `hybrid` or `blast-only` runs.

## Usage

```bash
organelle-cleaner assembly.gfa \
  --output-dir results
```

You can also invoke the existing script entry point directly:

```bash
python organelle_cleaner.py assembly.gfa \
  --output-dir results
```

Use hybrid mode when BLAST evidence is available:

```bash
organelle-cleaner assembly.gfa \
  --mode hybrid \
  --plastid-blast-tsv plastid_hits.tsv \
  --mit-blast-tsv mitochondrial_hits.tsv \
  --output-dir results
```

Accepted BLAST TSV columns are:

```text
qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send
```

BLAST TSV orientation is important: the organelle reference must be the query and the assembly contig must be the subject. In other words, `sseqid` must be the assembly contig ID, and `slen`, `sstart`, and `send` must describe the aligned span on that assembly contig.

Example `blastn` shape:

```bash
blastn -query organelle_refs.fa -subject assembly_contigs.fa \
  -outfmt '6 qseqid sseqid pident length slen mismatch gapopen qstart qend sstart send'
```

The tool computes BLAST support from merged subject spans on each assembly contig, so split or overlapping HSPs on the same contig are aggregated before support tiers are assigned. Strong and moderate BLAST support are assigned from merged coverage contributed by HSPs that meet that tier's identity threshold, which prevents a tiny high-identity HSP from certifying a much larger low-identity merged span.

`--plastid-blast-tsv` is the preferred plastid/chloroplast option name. The older `--chl-blast-tsv` flag is still accepted as a backward-compatible alias.

If you run `--mode hybrid` or `--mode blast-only`, you must provide at least one of `--plastid-blast-tsv`/`--chl-blast-tsv` or `--mit-blast-tsv`.
