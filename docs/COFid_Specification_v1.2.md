# COFid: A Canonical Naming Schema for Covalent Organic Frameworks

**Version:** 1.2.0
**Status:** Draft Specification
**Last Updated:** 2026-04-07

---

## 1. Overview

The **COFid** system provides a unique, machine-readable identifier for Covalent Organic Frameworks (COFs). Unlike MOFs, COFs lack distinct metal nodes, creating ambiguity in how components are ordered.

This schema resolves that ambiguity by enforcing a strict **Topological Hierarchy** for component ordering and utilizing a **Linkage Discriminator** to distinguish between materials formed from identical precursors (e.g., Imine vs. Amine COFs).

### Scope

This specification is designed for:
- Standard COFs with well-defined RCSR topologies
- Single-linkage-type frameworks
- Non-interpenetrated structures

**Out of Scope:**
- Mixed-linkage COFs (multiple linkage types in one framework)
- Non-RCSR topologies
- Interpenetrated (catenated) structures

---

## 2. Structure Format

The COFid is a single ASCII string composed of three mandatory blocks separated by double ampersands (`&&`).

```
[Monomer_Block]&&[Topology_Block]&&[Linkage_Block]
```

Each monomer in the Monomer Block has the format:

```
[connectivity]:[reactive_group]:[Canonical_SMILES]
```

### Complete Example

**TpPa-1 (Beta-ketoenamine COF):**

```
3:keto_aldehyde:O=Cc1c(O)c(C=O)c(O)c(C=O)c1O.2:amine:Nc1ccc(N)cc1&&hcb&&bken
```

**Breaking it down:**
- `3:keto_aldehyde:O=Cc1c(O)c(C=O)c(O)c(C=O)c1O` → Triformylphloroglucinol (3-connected, keto_aldehyde reactive group)
- `2:amine:Nc1ccc(N)cc1` → p-Phenylenediamine (2-connected, amine reactive group)
- `hcb` → Hexagonal topology
- `bken` → Beta-ketoenamine linkage

---

## 3. Block 1: The Monomer Block

**Format:** `[conn_1]:[group_1]:[SMILES_1].[conn_2]:[group_2]:[SMILES_2]...`

This block contains the annotated **Canonical SMILES** of the starting precursors. Each monomer includes:

1. **Connectivity count** (integer)
2. **Reactive group name** (string)
3. **Canonical SMILES** (RDKit-generated)

These three components are joined by colons (`:`), and multiple monomers are joined by periods (`.`).

### 3.1 SMILES Canonicalization

**CRITICAL:** All SMILES strings MUST be generated using **RDKit** with default canonicalization settings to ensure reproducibility across different implementations.

**Implementation:**
```python
from rdkit import Chem

def get_canonical_smiles(smiles_or_mol):
    """Generate canonical SMILES using RDKit."""
    if isinstance(smiles_or_mol, str):
        mol = Chem.MolFromSmiles(smiles_or_mol)
    else:
        mol = smiles_or_mol
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
```

**Settings:**
- Use `Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)`
- No isomeric SMILES (stereochemistry excluded)
- No explicit hydrogens unless required by the structure

### 3.2 Ordering Rules

To generate the COFid, monomers must be sorted using this two-step logic:

#### Step 1: Primary Sort by Connectivity (Degree)

- Calculate the number of reactive functional groups (connection points) for each monomer.
- **Descending Order:** The monomer with the **higher** connectivity (the Node) is placed **first**.
- *Example:* A 3-connected trialdehyde comes before a 2-connected diamine.

#### Step 2: Secondary Sort Alphabetically (Tie-Breaker)

- If connectivity is equal (e.g., both are 2-connected, or both are 4-connected), sort the SMILES strings **alphabetically (A-Z)**.
- *Example:* `Nc1...` comes before `O=Cc1...`

**Note:** Sorting is performed on the SMILES string only, not on the full annotated monomer string.

### 3.3 Reactive Group Definitions

The following table defines standard reactive groups and their connectivity counting rules:

| Reactive Group | Functional Group | Connectivity per Group | Example Monomer |
|----------------|------------------|------------------------|-----------------|
| `aldehyde` | -CHO | 1 | Benzaldehyde |
| `activated_methylene` | Activated methylene or methyl donor for vinylene-forming condensation | 1 | 2,4,6-Trimethyl-1,3,5-triazine |
| `amine` | -NH₂ | 1 | Aniline |
| `boronic_acid` | -B(OH)₂ | 1 | Phenylboronic acid |
| `catechol` | Two adjacent -OH | 1 (pair counts as 1) | 1,2-dihydroxybenzene |
| `carboxylic_acid` | -COOH | 1 | Benzoic acid |
| `hydrazide` | -C(=O)NHNH₂ | 1 | Benzohydrazide |
| `hydrazine` | -NH-NH₂ | 1 | Phenylhydrazine |
| `keto_aldehyde` | Ortho-hydroxylated aldehyde used in beta-ketoenamine COFs | 1 | Triformylphloroglucinol |
| `nitrile` | -CN | 1 | Benzonitrile |
| `hydroxyl` | -OH (isolated) | 1 | Phenol |
| `thiol` | -SH | 1 | Thiophenol |
| `anhydride` | Cyclic anhydride | 1 | Pyromellitic dianhydride |

**Connectivity Counting Example:**
- 1,3,5-Triformylbenzene has **3** aldehyde groups → **connectivity = 3**
- Triformylphloroglucinol has **3** keto_aldehyde groups → **connectivity = 3**
- p-Phenylenediamine has **2** amine groups → **connectivity = 2**
- Tetrakis(4-formylphenyl)methane has **4** aldehyde groups → **connectivity = 4**

### 3.4 Precursor Normalization

Before generating SMILES:

- **Desalting:** Remove counter-ions (e.g., Cl⁻, Na⁺) before generating SMILES.
- **Isotopes:** Use non-isotopic SMILES unless the isotope is critical to the material's definition.
- **Protonation State:** Use the neutral form of the molecule under standard conditions.
- **Tautomers:** Use the RDKit canonical tautomer.

---

## 4. Block 2: The Topology Block

**Format:** `[RCSR_Code]`

This block defines the underlying net topology of the framework.

### 4.1 Topology Code

- **Source:** Must correspond to a valid 3-letter code from the [RCSR Database](http://rcsr.net).
- **Case:** Lowercase (e.g., `hcb`, `sql`, `kgm`, `dia`).
- **Examples:**
  - `hcb` - Hexagonal (honeycomb)
  - `sql` - Square lattice
  - `dia` - Diamond
  - `kgm` - Kagome
  - `fxt` - Fxt net
  - `bor` - Boracite

**Note:** Stoichiometry is implicitly defined by the topology's node geometry and does not require explicit encoding.

---

## 5. Block 3: The Linkage Block (Discriminator)

**Format:** `[Linkage_Code]`

This block resolves the "Chemical Ambiguity" where the same precursors can form different linkages (e.g., redox states or tautomers). It acts as a metadata tag to differentiate the final material.

### 5.1 Linkage Code Guidelines

The linkage code is an implementation-defined discriminator that distinguishes COFs formed from identical precursors but different final linkages (e.g., due to redox state or tautomerization).

This specification does **not** define a controlled terminology list for linkage codes. Implementations MAY maintain their own curated lists and evolve them over time.

**Recommended constraints (for interoperability):**
- Use lowercase ASCII tokens (e.g., `imine`, `amine`, `bken`) as stable codes within a dataset.
- Do not include whitespace.
- Avoid reserved separators used elsewhere in COFid (`&&`, `.`, `:`).

---

## 6. Implementation Examples

### Scenario A: The Standard Imine COF

- **Precursors:**
  - Terephthalaldehyde (2-connected, aldehyde)
  - p-Phenylenediamine (2-connected, amine)
- **Topology:** Hexagonal (hcb)
- **Linkage:** Imine
- **Logic:** Degrees are equal (both 2-connected). Sort alphabetically: `Nc1...` < `O=Cc1...`

**COFid:**
```
2:amine:Nc1ccc(N)cc1.2:aldehyde:O=Cc1ccc(C=O)cc1&&hcb&&imine
```

### Scenario B: The Reduced Amine COF

- **Precursors:** Same as Scenario A
- **Reaction:** Post-synthetic reduction of the imine bond
- **Logic:** The precursors remain the identifying "parents." Only the linkage code changes.

**COFid:**
```
2:amine:Nc1ccc(N)cc1.2:aldehyde:O=Cc1ccc(C=O)cc1&&hcb&&amine
```

### Scenario C: The Keto-Enamine COF (TpPa-1)

- **Precursors:**
  - Triformylphloroglucinol (3-connected, keto_aldehyde)
  - p-Phenylenediamine (2-connected, amine)
- **Topology:** Hexagonal (hcb)
- **Linkage:** Beta-ketoenamine
- **Logic:** Degrees are unequal (3 vs 2). The 3-connected monomer goes first.

**COFid:**
```
3:keto_aldehyde:O=Cc1c(O)c(C=O)c(O)c(C=O)c1O.2:amine:Nc1ccc(N)cc1&&hcb&&bken
```

### Scenario D: Boronate Ester COF

- **Precursors:**
  - 1,3,5-Benzenetriboronic acid (3-connected, boronic_acid)
  - 2,3,6,7,10,11-Hexahydroxytriphenylene (3-connected, catechol)
- **Topology:** Hexagonal (hcb)
- **Linkage:** Boronate ester
- **Logic:** Both 3-connected. Sort alphabetically on SMILES.

**COFid (example):**
```
3:boronic_acid:OB(O)c1cc(B(O)O)cc(B(O)O)c1.3:catechol:Oc1cc2cc(O)c(O)cc2cc1&&hcb&&boest
```

### Scenario E: Diamond Topology COF

- **Precursors:**
  - Tetrakis(4-formylphenyl)methane (4-connected, aldehyde)
  - p-Phenylenediamine (2-connected, amine)
- **Topology:** Diamond (dia)
- **Linkage:** Imine
- **Logic:** 4-connected node comes first (4 > 2)

**COFid:**
```
4:aldehyde:O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1.2:amine:Nc1ccc(N)cc1&&dia&&imine
```

---

## 7. Implementation Guidelines

### 7.1 COFid Generation Algorithm

```python
from rdkit import Chem

def generate_cofid(monomers_data, topology, linkage):
    """
    Generate a COFid string.

    Args:
        monomers_data: List of dicts with keys:
            - 'smiles': Input SMILES
            - 'connectivity': Integer connectivity count
            - 'reactive_group': Reactive group name (e.g., 'aldehyde', 'amine')
        topology: RCSR code (e.g., 'hcb', 'dia')
        linkage: Linkage code (e.g., 'imine', 'bken')

    Returns:
        COFid string
    """
    # Step 1: Canonicalize SMILES
    for monomer in monomers_data:
        mol = Chem.MolFromSmiles(monomer['smiles'])
        if mol is None:
            raise ValueError(f"Invalid SMILES: {monomer['smiles']}")
        monomer['canonical_smiles'] = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)

    # Step 2: Sort by connectivity (descending), then alphabetically by SMILES
    sorted_monomers = sorted(
        monomers_data,
        key=lambda x: (-x['connectivity'], x['canonical_smiles'])
    )

    # Step 3: Build annotated monomer strings
    annotated_monomers = [
        f"{m['connectivity']}:{m['reactive_group']}:{m['canonical_smiles']}"
        for m in sorted_monomers
    ]

    # Step 4: Build COFid string
    monomer_block = '.'.join(annotated_monomers)
    cofid = f"{monomer_block}&&{topology}&&{linkage}"

    return cofid


# Example usage
monomers = [
    {
        'smiles': 'Oc1c(C=O)c(O)c(C=O)c(O)c1C=O',
        'connectivity': 3,
        'reactive_group': 'keto_aldehyde'
    },
    {
        'smiles': 'Nc1ccc(N)cc1',
        'connectivity': 2,
        'reactive_group': 'amine'
    }
]

cofid = generate_cofid(monomers, topology='hcb', linkage='bken')
print(cofid)
# Output: 3:keto_aldehyde:O=Cc1c(O)c(C=O)c(O)c(C=O)c1O.2:amine:Nc1ccc(N)cc1&&hcb&&bken
```

### 7.2 COFid Parsing Algorithm

```python
def parse_cofid(cofid_string):
    """
    Parse a COFid string into its components.

    Args:
        cofid_string: Valid COFid string

    Returns:
        Dictionary with keys:
            - 'monomers': List of monomer dictionaries
            - 'topology': RCSR code
            - 'linkage': Linkage code
    """
    # Split into three blocks
    blocks = cofid_string.split('&&')
    if len(blocks) != 3:
        raise ValueError("COFid must have exactly 3 blocks separated by &&")

    monomer_block, topology, linkage = blocks

    # Parse monomer block
    monomer_strings = monomer_block.split('.')
    monomers = []

    for mon_str in monomer_strings:
        parts = mon_str.split(':', 2)  # Split at most 2 times (colon delimiter)
        if len(parts) != 3:
            raise ValueError(f"Invalid monomer format: {mon_str}")

        connectivity, reactive_group, smiles = parts
        monomers.append({
            'connectivity': int(connectivity),
            'reactive_group': reactive_group,
            'canonical_smiles': smiles
        })

    return {
        'monomers': monomers,
        'topology': topology,
        'linkage': linkage
    }
```

### 7.3 Validation Checklist

Before accepting a COFid as valid:

- [ ] All SMILES are RDKit-canonical
- [ ] Monomers are sorted by connectivity (descending), then alphabetically by SMILES
- [ ] Each monomer has format `[int]:[string]:[SMILES]`
- [ ] Connectivity values are positive integers
- [ ] Reactive group names are in the allowed vocabulary
- [ ] Topology code exists in RCSR database
- [ ] Topology code is lowercase
- [ ] Linkage code is a non-empty discriminator string (no whitespace)
- [ ] All three blocks are present and separated by `&&`
- [ ] No whitespace in the COFid string

---

## 8. Versioning and Future Extensions

### Version History

- **v1.0.0** (Initial): Basic three-block structure with ordering rules
- **v1.1.0**:
  - Added RDKit specification
  - Embedded connectivity and reactive group into namespace
  - Removed interpenetration (out of scope)
  - Clarified scope boundaries
- **v1.1.1**:
  - Made RDKit canonicalization explicit (`isomericSmiles=False`)
  - Removed the linkage-code controlled terminology section (linkage codes are implementation-defined)
- **v1.2.0** (Current):
  - Corrected beta-ketoenamine examples to use the explicit `keto_aldehyde` reactive group
  - Expanded the reactive-group vocabulary to cover all currently implemented builtin monomer types
  - Aligned the published examples and quick reference with the current cofkit COFid implementation

### Planned Extensions

Future versions may include:
- Extended linkage vocabulary based on production data
- Extended reactive group vocabulary based on emerging COF chemistries
- Guidelines for handling edge cases as they arise in production

---

## 9. References

- **RCSR Database:** [http://rcsr.net](http://rcsr.net)
- **RDKit Documentation:** [https://www.rdkit.org/docs/](https://www.rdkit.org/docs/)
- **COF Literature:** Côté et al., Science 2005, 310, 1166-1170 (first COF report)

---

## Appendix A: Quick Reference

### COFid Structure

```
[conn_1]:[group_1]:[SMILES_1].[conn_2]:[group_2]:[SMILES_2]&&[topology]&&[linkage]
```

### Sorting Priority

1. **Connectivity** (high to low)
2. **Alphabetical by SMILES** (A to Z)

### Tools Required

- **RDKit** (SMILES canonicalization)
- **RCSR Database** (topology validation)

### Format Rules

- Monomer components separated by `:` (colon)
- Multiple monomers separated by `.` (period)
- Blocks separated by `&&` (double ampersand)
- No whitespace allowed
- Topology code in lowercase
- SMILES must be RDKit canonical

### Example Complete COFid

```
3:keto_aldehyde:O=Cc1c(O)c(C=O)c(O)c(C=O)c1O.2:amine:Nc1ccc(N)cc1&&hcb&&bken
```

**Breakdown:**
- Monomer 1: connectivity=3, group=keto_aldehyde, SMILES=O=Cc1c(O)c(C=O)c(O)c(C=O)c1O
- Monomer 2: connectivity=2, group=amine, SMILES=Nc1ccc(N)cc1
- Topology: hcb (hexagonal)
- Linkage: bken (beta-ketoenamine)

---

**Document Maintained By:** COFid Development Team
**Contact:** [Your contact information]
**License:** [Specify license if applicable]
