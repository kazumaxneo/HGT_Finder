# LKSZ1 / GTDB R232 example

Download this directory and open `sankey.html` in a browser. All five taxonomy ranks are available offline. Each ORF contributes one vote; the default quality floor is 0.5.

These are derived results for LKSZ1 (`GCF_040436315.1`), with 3,306 ORFs searched against 199,923 GTDB R232 representative genomes (618,638,921 proteins). The 2026-10-10 taxonomy display applies `(pident / 100) × query coverage × subject coverage ≥ 0.5` to the original top bit-score hits. It does not rerun DIAMOND or recalculate HGT scores.

`report.html`, `orf_scores.tsv` and `regions.tsv` retain the original 2026-10-09 HGT analysis before the quality filter. `run.json` distinguishes the original analysis parameters from the display quality threshold. `summary.json` contains counts for each display rank. Raw hits, reference sequences and private machine paths are not distributed here.

The nif interval is `contig_1:291513–313504` (27 ORFs). The original `candidate_7` spans a larger interval with 33 ORFs.

These plots show similarity-hit taxonomy. They do not establish HGT events, donor direction or contamination fractions.

Version 0.1.3 allows selecting taxa at any of the five ranks and comparing their species-best histograms (bit score on x, reference-species count on y) in `hit_distributions/` for five selected ORFs (33,388 reference-species hits). Use the taxon search/checklist to choose a host genus and another genus; selections persist separately for each rank. Overlay and panel layouts are available. ORF IDs in the Sankey and original score report link to these pages. Legacy hits have no stored identity/E-value; their quality is explicitly unavailable, while the Sankey top-hit quality was recovered separately by the original audit.
