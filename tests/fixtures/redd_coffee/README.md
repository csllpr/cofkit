# ReDD-COFFEE imine regressions

`A.cif.gz`, `B.cif.gz`, and `C.cif.gz` are unchanged CIF inputs from the
binary Imine_2D rerun, compressed losslessly. They were converted from
`redd-coffee/Imine_2D.tar.xz` with `chk_to_cofkit_cif.py --bond-orders kekule`.
`manifest.json` records the full source names and SHA-256 of each uncompressed
CIF. No runtime access to the external archive or `/tmp` is needed by tests.

The `expected_cuts` oracle is independent of the repair: these are bonds between
different `xx-xx-xx` force-field serials in each original CHK. Atom indices are
zero-based and preserved by the converter. Monomer SMILES in the manifest are
the repaired regression expectations; boundary-edge comparisons also check
against the independent CHK labels.

| Label | Atoms | CHK/CIF/RDKit bonds | Real imines | Previously missed | Expected precursor multiplicities |
|---|---:|---:|---:|---:|---|
| A: hnb, 30-02-04 + 01-03-04 | 528 | 584 | 24 | 12 | 8 trialdehyde + 12 diamine |
| B: bew, 29-02-04 + 02-03-04 | 444 | 486 | 24 | 10 | 6 tetraaldehyde + 12 diamine |
| C: dhc, 24-02-04 + 01-03-04 | 6996 | 7788 | 264 | 90 | 66 tetraaldehyde + 132 diamine |

The failure is a quinoid bond-order assignment from whole-framework valence
perception. **No periodic edges are lost in these fixtures.** Both the raw CHK
and the parsed CIF have the bond counts above. The previous 660-versus-584
comparison accidentally compared different structures.
