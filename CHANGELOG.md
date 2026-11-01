# Changelog

## 1.1.0

### Added
- `--utg-gfa`: unitig-isolation rules for `graph-only` and `hybrid` modes using the hifiasm unitig graph (`p_utg.gfa`, or `r_utg.gfa`) of the same assembly. Contig reads are matched to unitigs via GFA `A` lines; `iso_frac` is the fraction of a contig's unitigs in components isolated from the rest of the graph.
  - hybrid: strong BLAST (identity >= 80%, HSP covering >= 20% of the contig) OR (weak BLAST: identity >= 80%, alignment >= 500 bp AND `iso_frac >= 0.5` AND length <= L).
  - graph-only: `iso_frac >= 0.5` AND length <= L.
  - HSPs with e-value > 1e-20 are ignored.
  - L = min(2 x the longest organelle reference FASTA total length, 2 Mb); 2 Mb without references.
  - New `report.tsv` columns: `utg_used`, `iso_frac`, `utg_max_len`, `utg_rule`.
- `organelle_assembly.fa` output with the contigs classified as organelle.

### Fixed
- Depth parsing: segment depth is now read from assembler-specific GFA tags (`RC:i`/LN, hifiasm `rd:i`, `dp`/`DP`, `KC:i`/LN). Previously only `RC:i` was read, so hifiasm graphs had no coverage and the coverage signal never fired.
- Internal BLAST target limit: `blastn -max_target_seqs` is now max(500, number of assembly contigs). Contigs are the BLAST subjects, so with the blastn default of 500 an organelle reference matching more than 500 contigs silently lost support for the rest.
- `organelle_contigs.txt`, `nuclear_contigs.txt` and the FASTA outputs now follow the `classification` column of `report.tsv`; contigs without sequence are excluded from FASTA output with a warning.

## 0.1.0
- Initial release.
