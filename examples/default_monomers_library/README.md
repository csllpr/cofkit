# Default Monomers Library

This directory is auto-generated using the current monomer role detector.

- source_dir: `/share/data/csllpr/rctest/cofkit/examples/batch_test_monomers`
- registered monomers: 440
- failed or ambiguous detections: 0
- source-label mismatches among registered monomers: 26
- registered via filename fallback (multi-role monomers): 27
- omitted monofunctional (connectivity-1) monomers: 3

## Registered Libraries

- `aldehydes_count_2.txt`: 134
- `aldehydes_count_3.txt`: 28
- `aldehydes_count_4.txt`: 6
- `aldehydes_count_6.txt`: 8
- `amines_count_2.txt`: 180
- `amines_count_3.txt`: 46
- `amines_count_4.txt`: 6
- `amines_count_6.txt`: 5
- `hydrazides_count_2.txt`: 4
- `keto_aldehydes_count_2.txt`: 17
- `keto_aldehydes_count_3.txt`: 4
- `keto_aldehydes_count_4.txt`: 2

## Metadata Files

- `registry.jsonl`: registered monomers with detector metadata and source provenance
- `failures.jsonl`: failed or ambiguous autodetections, plus monofunctional (connectivity-1) monomers omitted by default (status `skipped`; register them with `--include-monofunctional`)
