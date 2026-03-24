# Default Monomers Library

This directory is auto-generated from `examples/batch_test_monomers/` using the
current monomer role detector. Successful autodetections are regrouped into explicit
role/count libraries, while failures and ambiguous cases are recorded separately.

It is already grouped detector output. Use this directory directly with the batch
commands and do not add `--auto-detect-libraries`.

- source_dir: `examples/batch_test_monomers`
- intended batch templates from this snapshot: `hydrazone_bridge`, `imine_bridge`, `keto_enamine_bridge`
- registered monomers: 489
- failed or ambiguous detections: 23
- source-label mismatches among registered monomers: 39

## Registered Libraries

- `aldehydes_count_2.txt`: 149
- `aldehydes_count_3.txt`: 30
- `aldehydes_count_4.txt`: 7
- `aldehydes_count_6.txt`: 7
- `amines_count_2.txt`: 197
- `amines_count_3.txt`: 49
- `amines_count_4.txt`: 7
- `amines_count_6.txt`: 5
- `hydrazides_count_2.txt`: 4
- `keto_aldehydes_count_1.txt`: 4
- `keto_aldehydes_count_2.txt`: 22
- `keto_aldehydes_count_3.txt`: 6
- `keto_aldehydes_count_4.txt`: 2

## Metadata Files

- `registry.jsonl`: registered monomers with detector metadata and source provenance
- `failures.jsonl`: failed or ambiguous autodetections
