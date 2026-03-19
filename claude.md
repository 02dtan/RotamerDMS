# RotamerDMS Implementation Documentation

This document tracks implementation choices, differences from `design.md`, and provides usage instructions.

## Overview

RotamerDMS is a two-phase algorithm for maximizing binding pocket volume in protein structures through rotamer sampling:
1. **Mutation Phase**: Identify residues whose mutation to Alanine increases pocket volume
2. **Rotamer Sampling Phase**: Sample rotamer states for impactful residues to maximize cavity volume

## Implementation Choices

### Project Structure

```
RotamerDMS/
├── src/                          # Source modules
│   ├── __init__.py
│   ├── preparation.py            # Structure loading, sequence validation, binding site identification
│   ├── cavity.py                 # pyKVfinder cavity detection and volume measurement
│   ├── mutation.py               # Alanine mutation phase
│   └── rotamer_sampling.py       # Rotamer sampling phase
├── data/                         # Input structures
├── output/                       # Modified structures
├── checkpoints/                  # Intermediate checkpoints and statistics
├── run_rotamer_dms.py            # Main runner script
├── run_rotamer_dms.sh            # Shell wrapper with environment setup
├── design.md                     # Original design specifications
└── claude.md                     # This documentation file
```

### Key Design Decisions

| Decision | design.md Specification | Implementation |
|----------|------------------------|----------------|
| Rotamer frequency | Extract from Schrodinger rotamer library | Using `rotamer.percentage` attribute |
| Cavity identification | Cavities contacting ≥2 binding site residues | Implemented as specified |
| Residue selection | Within 7Å of ligand | Using Schrodinger atom distance calculations |
| Volume metric | Sum of relevant cavity volumes | Implemented as specified |

---

## Module Documentation

### 1. `src/preparation.py`

*Status: Pending implementation*

**Purpose**: Load structures, validate sequence identity, identify binding site residues.

**Key Functions**:
- `load_structure(pdb_path)`: Load PDB using Schrodinger StructureReader
- `validate_sequence_identity(ref_struct, sample_struct)`: Ensure 100% sequence identity
- `get_binding_site_residues(structure, ligand_code, distance_threshold)`: Find residues within threshold of ligand

---

### 2. `src/cavity.py`

*Status: Implemented*

**Purpose**: Cavity detection and volume measurement using pyKVfinder.

**Parameters**: Matched to ChimeraX "Find Cavities" defaults (from [ChimeraX docs](https://www.cgl.ucsf.edu/chimerax/docs/user/tools/findcavities.html)):
| Parameter | Value | pyKVFinder Name |
|-----------|-------|-----------------|
| Grid spacing | 0.60 Å | `step` |
| Inner probe radius | 1.4 Å | `probe_in` |
| Outer probe radius | 4.0 Å | `probe_out` |
| Exterior trim distance | 2.4 Å | `removal_distance` |
| Minimum cavity volume | 5.0 ų | `volume_cutoff` |

**Key Functions**:
- `detect_cavities(struct, params)`: Run pyKVFinder on a structure, returns volumes and residue contacts
- `get_cavities_contacting_residues(results, binding_residues, min_contacts=4)`: Filter cavities by binding site contact
- `measure_binding_site_volume(struct, binding_residues)`: Main function for volume measurement in optimization loop
- `get_atom_cavity_distances(struct, results, binding_residues)`: Cross-reference atom coordinates with cavity points

**Implementation Note**: pyKVFinder residue format is "CHAIN_RESNAME_RESNUM" (e.g., "A_MET_102"). The module includes parsing to convert to our standard (chain, resnum, resname) tuple format.

---

### 3. `src/mutation.py`

*Status: Implemented*

**Purpose**: Mutation phase - mutate residues to Alanine and measure volume impact.

**Key Functions**:
- `mutate_residue_to_ala(struct, chain, resnum)`: Mutate single residue to Alanine using Schrodinger `build.mutate()`
- `run_mutation_phase(sample_struct, binding_site_residues, checkpoint_dir)`: Main mutation phase driver
- `save_mutation_checkpoint(result, checkpoint_dir)`: Save results as JSON checkpoint
- `select_top_residues(sorted_residues, volume_changes, top_percent)`: Filter to top N% of volume-affecting residues

**Algorithm** (per design.md):
1. Measure initial binding site volume
2. For each binding site residue:
   - Create copy of structure
   - Mutate to Alanine
   - Measure new volume
   - If volume INCREASED: keep residue and record delta
   - If volume unchanged/decreased: remove from list
3. Sort kept residues by volume increase (descending)
4. Save checkpoint with all statistics

**Implementation Notes**:
- Skips ALA and GLY residues (no sidechain to remove)
- Uses `struct.copy()` to avoid modifying original structure during testing
- Checkpoint saved as JSON with timestamp

---

### 4. `src/rotamer_sampling.py`

*Status: Implemented*

**Purpose**: Rotamer sampling phase - sample states and select volume-maximizing rotamers.

**Key Functions**:
- `get_rotamer_library(struct, chain, resnum)`: Load Schrodinger backbone-dependent rotamer library
- `sample_rotamers_for_residue(...)`: Test all rotamer states for one residue, rank by volume impact
- `run_rotamer_sampling_phase(...)`: Main driver - process residues in priority order
- `apply_best_rotamers(struct, best_rotamers)`: Apply selected rotamers to structure
- `finalize_structure(...)`: Verify volume increase, check cavity count, save output

**Algorithm** (per design.md):
1. Process residues in order of volume impact (from mutation phase)
2. For each residue:
   - Load backbone-dependent rotamer library
   - Test each rotamer state, measure volume change
   - Record `rotamer.percentage` for each state
   - Sort states by volume contribution (descending)
3. Apply best rotamer for each residue to working structure
4. Final verification with pyKVfinder
5. Save with appropriate filename (flag if multiple pockets persist)

**Implementation Notes**:
- Uses `rotamer.percentage` attribute for frequency proxy
- Accumulates best rotamers progressively (each residue sampled in context of prior changes)
- Saves checkpoint JSON with per-residue rotamer details

---

## Usage Tutorial

### Prerequisites

1. Activate the Schrodinger environment:
```bash
export SB_BASE_OVERRIDE=/lustre/fs4/ruit/scratch/rebecca/LYU_SBGRID/INST
source /lustre/fs4/ruit/scratch/rebecca/LYU_SBGRID/INST/sbgrid.shrc
```

2. Activate the conda environment:
```bash
conda activate rotamer_dms
```

### Running the Algorithm

**Method 1: Using the shell wrapper (recommended)**

The shell wrapper handles environment setup automatically:

```bash
cd /lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS

# Make executable (first time only)
chmod +x run_rotamer_dms.sh

# Basic usage
./run_rotamer_dms.sh \
    --reference data/4W59_n-hexylbenzene.pdb \
    --sample data/s0067_rec.pdb \
    --ligand 3GZ

# With custom options
./run_rotamer_dms.sh \
    --reference data/4W59_n-hexylbenzene.pdb \
    --sample data/s0067_rec.pdb \
    --ligand 3GZ \
    --distance 5.0 \
    --top-percent 50
```

**Method 2: Direct Python execution**

If environment is already set up:

```bash
$SCHRODINGER/run python3 run_rotamer_dms.py \
    --reference data/4W59_n-hexylbenzene.pdb \
    --sample data/s0067_rec.pdb \
    --ligand 3GZ \
    --top-percent 50
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reference`, `-r` | Path to experimental reference structure with bound ligand | Required |
| `--sample`, `-s` | Path to sample structure to modify | Required |
| `--ligand`, `-l` | 3-letter ligand code in reference structure | Required |
| `--distance`, `-d` | Distance threshold (Å) for binding site residue selection | 7.0 |
| `--top-percent`, `-t` | Percentage of top volume-affecting residues to sample | 100 |
| `--output-dir`, `-o` | Output directory | `./output` |
| `--checkpoint-dir`, `-c` | Checkpoint directory | `./checkpoints` |
| `--quiet`, `-q` | Reduce output verbosity | False |

### Output Files

**Output directory (`output/`):**
- `<sample_name>_optimized.pdb` - Final optimized structure (single cavity achieved)
- `<sample_name>_optimized_still_multiple_pockets.pdb` - Structure if multiple cavities persist

**Checkpoint directory (`checkpoints/`):**
- `mutation_phase_<timestamp>.json` - Mutation phase results:
  - Initial volume
  - Sorted residues by volume impact
  - Volume changes per residue
  - Removed residues
- `rotamer_sampling_<timestamp>.json` - Rotamer sampling results:
  - Best rotamer index per residue
  - Volume changes per rotamer
  - Rotamer percentages (frequencies)

### Example Workflow

```bash
# 1. Navigate to RotamerDMS directory
cd /lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS

# 2. Run with default settings
./run_rotamer_dms.sh \
    --reference data/4W59_n-hexylbenzene.pdb \
    --sample data/s0067_rec.pdb \
    --ligand 3GZ

# 3. Check output
ls -la output/
# s0067_rec_optimized.pdb  (or s0067_rec_optimized_still_multiple_pockets.pdb)

# 4. Review checkpoints for statistics
cat checkpoints/mutation_phase_*.json
cat checkpoints/rotamer_sampling_*.json
```

### Algorithm Summary

1. **Preparation**: Load structures, validate sequence identity, identify binding site residues within 7Å of ligand
2. **Mutation Phase**: Mutate each residue to Ala, keep only those that increase pocket volume, sort by impact
3. **Rotamer Sampling**: For each sorted residue, test all rotamer states, select volume-maximizing states
4. **Finalization**: Apply best rotamers, verify volume increase, save output

---

## Changelog

- **v0.1** (Complete): Initial implementation following design.md specifications
  - Preparation module with sequence validation
  - Cavity detection using pyKVFinder with ChimeraX-consistent parameters
  - Mutation phase with Ala scanning
  - Rotamer sampling with percentage extraction
  - Main runner script and shell wrapper
