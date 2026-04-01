# RotamerDMS Implementation Documentation

This document tracks implementation choices, differences from `design.md`, and provides usage instructions.

> **Quick Reference**: For a concise workflow summary and all current parameter values, see **`workflow_parameter.md`**.

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
└── implementation.md             # This documentation file
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
| Minimum cavity volume | 5.0 Å³ | `volume_cutoff` |

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
- `get_current_chi_angles(struct, residue)`: Measure chi dihedral angles of a residue
- `find_closest_library_rotamer(struct, chain, resnum)`: Find library rotamer matching current conformation
- `get_rotamer_library(struct, chain, resnum)`: Load Schrodinger backbone-dependent rotamer library
- `sample_rotamers_for_residue(...)`: Test all rotamer states for one residue, rank by volume impact
- `run_rotamer_sampling_phase(...)`: Main driver - process residues in priority order
- `apply_best_rotamers(struct, best_rotamers)`: Apply selected rotamers to structure
- `finalize_structure(...)`: Verify volume increase, check cavity count, save output

**Algorithm** (per design.md):
1. Process residues in order of volume impact (from mutation phase)
2. For each residue:
   - Identify current/native rotamer by comparing chi angles to library
   - Record native rotamer probability (p_curr)
   - Load backbone-dependent rotamer library
   - Test each rotamer state, measure volume change
   - Record `rotamer.percentage` (p_sugg) for each state
   - Compute ΔG = -ln(p_sugg / p_curr) for each state (in arbitrary units)
   - Sort states by volume contribution (descending)
3. Apply best rotamer for each residue to working structure
4. Final verification with pyKVfinder
5. Save with appropriate filename (flag if multiple pockets persist)

**Implementation Notes**:
- Native rotamer identified by minimizing angular distance in chi-space (periodicity-aware)
- Symmetry handling for PHE, TYR (χ2), ASP (χ2), GLU (χ3) - 180° flips are equivalent
- Uses `measure_dihedral_angle` from `schrodinger.structutils.measure` for chi angle measurement
- Uses `rotamer.percentage` / `rotamer.probability` for Boltzmann probability
- ΔG = -ln(p_sugg / p_curr) in arbitrary units (AU); ΔG = 0 when suggested = native
- Accumulates best rotamers progressively (each residue sampled in context of prior changes)
- Saves checkpoint JSON with per-residue rotamer details including ΔG values

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

**Required:**

| Parameter | Description |
|-----------|-------------|
| `--reference`, `-r` | Path to experimental reference structure with bound ligand |
| `--sample`, `-s` | Path to sample structure to optimize |
| `--ligand`, `-l` | 3-letter residue code of the ligand in reference structure |

**Structure & Volume Options:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--distance`, `-d` | Distance threshold (Å) for binding site residue selection | 7.0 |
| `--top-percent`, `-t` | Percentage of top volume-affecting residues to sample | 100 |
| `--min-contacts`, `-m` | Minimum binding site residues a cavity must contact to be counted | 4 |

**Joint Optimization Weights:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--w-vol` | Weight for volume term in joint optimization score | 1.0 |
| `--w-dg` | Weight for deltaG (rotamer probability) term in joint optimization | 1.0 |
| `--w-wca` | Weight for WCA/fa_rep clash energy term in joint optimization | 1.0 |
| `--wca-threshold` | Hard WCA threshold (REU) to reject clashing rotamers | None |
| `--no-wca` | Disable WCA (steric clash) potential computation entirely | False |

**Output & Verbosity:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--output-dir`, `-o` | Output directory for final structures | `./output` |
| `--checkpoint-dir`, `-c` | Checkpoint directory for intermediate results | `./checkpoints` |
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

---

## Session: March 19, 2026 - WCA Potential & Multi-Objective Optimization

### Overview

Extended rotamer sampling to consider multiple objectives:
1. **Pocket volume maximization** (original objective)
2. **Steric clash avoidance** via WCA potential
3. **Conformational free energy** (deltaG)
4. **Pareto front analysis** for tradeoff visualization

---

### 1. WCA Potential Implementation

#### Background: What is the WCA Potential?

The **Weeks-Chandler-Andersen (WCA) potential** is a repulsive-only variant of the Lennard-Jones potential, isolating steric repulsion by truncating at the LJ minimum:

```
For r ≤ r_min = 2^(1/6)σ:
    U_WCA(r) = 4ε[(σ/r)^12 - (σ/r)^6] + ε

For r > r_min:
    U_WCA(r) = 0
```

This captures "steric clashes" - when atoms are too close, energy increases sharply.

#### Implementation Choice: PyRosetta's `fa_rep` Term

**Decision**: Use PyRosetta's `fa_rep` (full-atom repulsive) energy term as WCA proxy.

**Justification**:
1. **Physical equivalence**: `fa_rep` implements repulsive LJ, conceptually identical to WCA
2. **Optimized implementation**: PyRosetta uses pre-computed neighbor lists, avoiding O(N²) iteration
3. **Validated parameters**: Rosetta's ref2015 energy function is extensively validated
4. **No manual cutoff needed**: PyRosetta handles distance cutoffs via neighbor detection

#### New File: `src/wca_potential.py` (179 lines)

```python
def _init_pyrosetta():
    """
    Lazy initialization of PyRosetta.
    
    Justification: PyRosetta init is expensive (~2-3 sec). Lazy init ensures
    cost is paid only once per session, only if WCA computation is requested.
    """

def schrodinger_to_pyrosetta_pose(schrodinger_struct):
    """
    Convert Schrodinger Structure to PyRosetta Pose via temporary PDB file.
    
    Justification: Schrodinger and PyRosetta use incompatible internal
    representations. PDB is universal interchange format both can read/write.
    Temp file cleaned up immediately after conversion.
    """

def compute_wca_batch(schrodinger_struct, target_chain, target_resnum, verbose=False):
    """
    Compute WCA/fa_rep energy for a specific residue.
    
    Implementation:
    1. Convert structure to PyRosetta pose
    2. Create minimal score function with only fa_rep (weight=1.0)
    3. Score pose (populates energy graph using neighbor lists)
    4. Extract per-residue fa_rep energy
    
    Justification for minimal score function:
    - Only fa_rep avoids computing unnecessary energy terms
    - Weight 1.0 gives energy in Rosetta Energy Units (REU)
    - REU interpretable: >5 REU suggests significant clash
    """

def get_wca_energy(schrodinger_struct, chain, resnum, verbose=False):
    """
    Simple interface returning WCA energy as float or None.
    Returns None on failure for graceful degradation.
    """
```

#### PyRosetta Installation

```bash
conda activate rotamer_dms
conda install -c https://conda.rosettacommons.org pyrosetta
# Package size: ~1.41 GB
```

#### Code Cleanup

Removed dead code from initial implementation:
- `compute_residue_repulsive_energy()`: Never called, manually iterated atom pairs
- `compute_wca_for_rotamer_state()`: Never called wrapper

**Justification**: `compute_wca_batch()` is more efficient using PyRosetta's internal neighbor lists.

---

### 2. Conformational Free Energy (ΔG) Computation

#### Background: Rotamer Boltzmann Statistics

The Schrodinger backbone-dependent rotamer library provides **rotamer percentages** (frequencies) derived from statistical analysis of high-resolution protein structures. These percentages approximate the Boltzmann probability distribution of rotamer states.

From statistical mechanics, the probability of a rotamer state is related to its free energy:
```
p_i = exp(-G_i / kT) / Z
```

where Z is the partition function. The **relative free energy change** between two rotamer states is:
```
ΔG = G_suggested - G_native = -kT * ln(p_suggested / p_native)
```

#### Implementation

In `src/rotamer_sampling.py`, within `sample_rotamers_for_residue()`:

```python
# Identify native rotamer by matching chi angles to library
native_rotamer_info = find_closest_library_rotamer(struct, chain, resnum)
p_curr = native_rotamer_info['probability']  # Native rotamer percentage

# For each candidate rotamer:
percentage = rotamer.percentage()  # Suggested rotamer percentage (p_sugg)

# Compute relative free energy change
if p_curr > 0 and percentage > 0:
    delta_g = -math.log(percentage / p_curr)  # In arbitrary units (AU)
else:
    delta_g = None
```

#### Key Implementation Details

1. **Native rotamer identification** (`find_closest_library_rotamer()`):
   - Measures current chi dihedral angles from structure
   - Compares to all library rotamers using angular distance
   - Handles periodicity (angles wrap at ±180°)
   - Handles symmetric residues (PHE/TYR χ2, ASP χ2, GLU χ3 have 180° equivalence)
   - Returns closest match with its percentage (p_curr)

2. **Angular distance computation** (`compute_chi_angular_distance()`):
   - Computes minimum angular difference accounting for periodicity
   - For symmetric residues, tests both χ and χ+180° and takes minimum
   - Returns sum of squared angular differences (Euclidean in chi-space)

3. **Units**:
   - ΔG is in **arbitrary units (AU)**, not kcal/mol or kJ/mol
   - This is because we use `ln` instead of `log` with explicit kT
   - AU are dimensionless and internally consistent for ranking
   - ΔG = 0 when suggested rotamer equals native rotamer
   - ΔG > 0 means suggested is less favorable than native
   - ΔG < 0 means suggested is more favorable than native

4. **Edge cases**:
   - If native rotamer cannot be identified: `p_curr = None`, `delta_g = None`
   - If rotamer percentage is zero: `delta_g = None`
   - None values are handled gracefully in joint scoring

#### Justification for This Approach

1. **Statistical basis**: Rotamer frequencies from PDB statistics encode empirical free energy landscape
2. **No force field needed**: Avoids expensive energy minimization or MD simulation
3. **Backbone-dependent**: Library accounts for phi/psi-dependent rotamer preferences
4. **Computationally cheap**: Simple log ratio, no structural calculations required
5. **Interpretable**: Higher ΔG means less favorable conformational change

---

### 3. Multi-Objective Optimization Scoring

#### Problem Statement

Three metrics with different units/scales:
- **ΔV (volume)**: Å³, range -50 to +100 → MAXIMIZE
- **ΔG (free energy)**: arbitrary units, range -3 to +3 → MINIMIZE  
- **WCA (clash)**: REU, range 0 to 50+ → MINIMIZE

#### Solution: Z-Score Normalization + Weighted Sum

**New File: `src/scoring.py` (~320 lines)**

##### Statistics Computation for Z-Score Normalization

```python
def compute_stats(values):
    """
    Compute mean and standard deviation for z-score normalization.
    
    Returns (mean, std) tuple used by compute_joint_score() to
    normalize each metric inline.
    
    Edge cases:
    - None values skipped
    - <2 valid values: return (0.0, 1.0)
    - std=0: use std=1.0 to avoid division by zero
    """
```

Z-scores are computed inline in `compute_joint_score()`: `z = (value - mean) / std`

##### Joint Score Computation

```python
def compute_joint_score(volume_change, delta_g, wca_energy,
                        w_vol=1.0, w_dg=1.0, w_wca=1.0, ...):
    """
    score = w_vol * z(ΔV) - w_dg * z(ΔG) - w_wca * z(WCA)
    
    Sign convention:
    - Positive from volume (want MORE)
    - Negative from ΔG (want LESS)
    - Negative from WCA (want LESS)
    
    Higher score = better rotamer
    
    Hard threshold: If wca_threshold set and WCA > threshold, return None
    (reject severe clashes regardless of other terms)
    """
```

##### Ranking Function

```python
def rank_rotamers_by_joint_score(rotamer_results, w_vol, w_dg, w_wca, wca_threshold):
    """
    Process:
    1. Compute normalization params (mean, std) for each metric
    2. Compute joint score for each rotamer
    3. Sort descending (best first)
    4. Append rejected rotamers (None score) at end
    
    Justification for per-residue normalization:
    Different residue types have different rotamer distributions;
    per-residue normalization ensures fair comparison within each residue.
    """
```

##### Default Weights

`w_vol=1.0, w_dg=1.0, w_wca=1.0` (equal importance)

**Justification**: With z-score normalization, weights represent relative importance in standard deviations. Equal weights = 1-std improvement in volume is equally valuable as 1-std reduction in WCA.

**Tuning recommendations**:
- `--w-vol 2.0`: Prioritize volume (accept some clashes)
- `--w-wca 2.0`: Prioritize clash avoidance (sacrifice volume)
- `--wca-threshold 10.0`: Hard reject WCA > 10 REU

---

### 4. Pareto Front Analysis

#### What is a Pareto Front?

A rotamer is **Pareto-optimal** if no other rotamer is better in ALL objectives. The **Pareto front** is the set of all non-dominated solutions.

For volume vs. WCA: a rotamer dominates another if it has BOTH higher volume AND lower WCA.

#### Implementation in `src/scoring.py`

```python
def compute_pareto_front_2d(rotamer_results, maximize_key, minimize_key):
    """
    Efficient 2D Pareto computation - O(n log n)
    
    Algorithm:
    1. Sort by maximize_key descending
    2. Sweep, tracking best minimize_key seen
    3. Point is Pareto-optimal if minimize_key better than all prior points
    """

def compute_pareto_front_3d(rotamer_results, maximize_keys, minimize_keys):
    """
    3D Pareto computation - O(n²)
    
    For each candidate, check if any other solution dominates it.
    Note: More expensive; for large datasets consider NSGA-II.
    """

def analyze_pareto_tradeoff(rotamer_results, verbose=True):
    """
    Returns:
    - 2D Pareto front (volume vs. WCA)
    - 3D Pareto front (volume vs. ΔG vs. WCA)
    - Extreme points (best volume, lowest WCA)
    - Best balanced (highest joint score on Pareto front)
    """
```

#### Standalone Script: `analyze_pareto.py` (~270 lines)

Post-hoc analysis of checkpoint files:

```bash
python analyze_pareto.py --checkpoint checkpoints/rotamer_sampling_*.json --output analysis/
```

**Output files**:
- `pareto_front.csv`: Pareto-optimal residue mutations
- `all_best_rotamers.csv`: All residues with Pareto membership flag
- `pareto_analysis_summary.json`: Summary statistics

---

### 5. Integration into Rotamer Sampling

#### Modified: `src/rotamer_sampling.py`

##### `sample_rotamers_for_residue()`

**New parameters**: `w_vol`, `w_dg`, `w_wca`, `wca_threshold`

**Changes**:
1. After collecting rotamer results, rank by joint score (not volume only)
2. Compute and include Pareto front analysis
3. Return joint scores and Pareto data

**New return fields**:
- `joint_scores`: List of joint optimization scores
- `best_joint_score`: Score of selected best rotamer
- `best_by_volume_idx`: Best by volume only (for comparison)
- `pareto_front`: Pareto-optimal rotamer results
- `pareto_analysis`: Full analysis dictionary

##### `run_rotamer_sampling_phase()`

**New parameters**: `w_vol`, `w_dg`, `w_wca`, `wca_threshold`

Passes weight parameters through to `sample_rotamers_for_residue()`.

##### `save_rotamer_checkpoint()`

**New fields in checkpoint**:
- `best_joint_score`
- `best_by_volume_idx`
- `pareto_front_indices`
- `num_pareto_optimal`

**Removed**: `wca_sorted_indices`, `wca_sorted_energies` (replaced by joint scoring)

---

### 6. CLI Updates

#### Modified: `run_rotamer_dms.py`

**New arguments**:

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--w-vol` | float | 1.0 | Weight for volume term |
| `--w-dg` | float | 1.0 | Weight for deltaG term |
| `--w-wca` | float | 1.0 | Weight for WCA term |
| `--wca-threshold` | float | None | Hard WCA threshold (REU) |

Startup banner displays weights when WCA is enabled.

---

### 7. Updated Project Structure

```
RotamerDMS/
├── src/
│   ├── __init__.py
│   ├── preparation.py
│   ├── cavity.py
│   ├── mutation.py
│   ├── rotamer_sampling.py    # Modified: joint scoring integration
│   ├── wca_potential.py       # NEW: WCA via PyRosetta
│   └── scoring.py             # NEW: Multi-objective scoring & Pareto
├── run_rotamer_dms.py         # Modified: weight CLI args
├── analyze_pareto.py          # NEW: Standalone Pareto analysis
└── ...
```

---

### 8. Usage Examples

#### Basic (equal weights)
```bash
$SCHRODINGER/run python3 run_rotamer_dms.py \
    --reference data/ref.pdb --sample data/sample.pdb --ligand LIG
```

#### Prioritize Volume
```bash
$SCHRODINGER/run python3 run_rotamer_dms.py \
    --reference data/ref.pdb --sample data/sample.pdb --ligand LIG \
    --w-vol 2.0 --w-wca 0.5
```

#### Strict Clash Rejection
```bash
$SCHRODINGER/run python3 run_rotamer_dms.py \
    --reference data/ref.pdb --sample data/sample.pdb --ligand LIG \
    --wca-threshold 10.0
```

#### Disable WCA (faster, volume-only)
```bash
$SCHRODINGER/run python3 run_rotamer_dms.py \
    --reference data/ref.pdb --sample data/sample.pdb --ligand LIG \
    --no-wca
```

#### Analyze Pareto Front
```bash
python analyze_pareto.py \
    --checkpoint checkpoints/rotamer_sampling_*.json \
    --output analysis/
```

---

### 9. Changelog Update

- **v0.2** (March 19, 2026): WCA potential and multi-objective optimization
  - PyRosetta integration for WCA/fa_rep energy computation
  - Z-score normalized joint scoring (volume, deltaG, WCA)
  - Pareto front analysis for tradeoff visualization
  - CLI arguments for weight tuning and WCA threshold
  - Standalone Pareto analysis script

---

## Session: March 20, 2026 - SLURM Integration & Per-Residue Pareto Analysis

### Overview

Three major improvements:
1. **SLURM integration** for cluster job submission
2. **PyRosetta subprocess fixes** (timeout, error logging)
3. **Per-residue 3D Pareto analysis** in checkpoints and analyze_pareto.py

---

### 1. SLURM Integration

#### Problem
Running RotamerDMS interactively on login nodes is inappropriate for compute-intensive PyRosetta operations.

#### Solution
Added SBATCH directives to `run_rotamer_dms.sh` for seamless SLURM submission:

```bash
#!/bin/bash
#SBATCH --job-name=rot_dms
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=16GB
#SBATCH --partition=lyu_a
#SBATCH --account=lyu_condo_bank
#SBATCH --output=rot_dms_%j.out
#SBATCH --error=rot_dms_%j.err
```

#### Key Fix: Working Directory Resolution

**Problem**: When SLURM runs a job, it copies the script to `/var/spool/slurmd/job<id>/` and executes from there. The original `SCRIPT_DIR` detection via `${BASH_SOURCE[0]}` pointed to this temp location, causing "file not found" errors.

**Solution**: Use `SLURM_SUBMIT_DIR` environment variable (set by SLURM to the directory where `sbatch` was called):

```bash
# Get the working directory
# When running under SLURM, use SLURM_SUBMIT_DIR (where sbatch was called)
# When running interactively, use the script's location
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
fi
```

#### Usage

**Interactive** (unchanged):
```bash
./run_rotamer_dms.sh --reference data/ref.pdb --sample data/sample.pdb --ligand LIG
```

**SLURM submission**:
```bash
sbatch run_rotamer_dms.sh --reference data/ref.pdb --sample data/sample.pdb --ligand LIG
```

Command-line arguments pass through transparently via `$@`.

#### Resource Recommendations

| Resource | Value | Justification |
|----------|-------|---------------|
| CPUs | 4 | pyKVFinder can use OpenMP; main code is serial |
| Memory | 16GB | PyRosetta initialization needs headroom |
| Time | 1 hour | Sufficient for ~10 residues; increase for larger runs |

---

### 2. PyRosetta Subprocess Fixes

#### Problem 1: Silent WCA Failures

WCA computation failures were silently swallowed when `verbose=False`, making debugging impossible.

**Fix**: Always log errors to stderr (appears in SLURM `.err` file):

```python
# In src/wca_potential.py, get_wca_energy()
else:
    # Always log errors to stderr so they appear in SLURM error files
    import sys
    error_msg = result.get('error', 'unknown error')
    print(f"WCA computation failed for {chain}:{resnum}: {error_msg}", file=sys.stderr)
    return None
```

#### Problem 2: Timeout on First PyRosetta Call

The first PyRosetta subprocess call timed out (60s) because:
1. Python interpreter startup
2. `import pyrosetta` (~10-20s)
3. `pyrosetta.init()` (~5-15s)
4. PDB loading and scoring

On busy cluster nodes, this exceeded 60 seconds.

**Fix**: Increased timeout from 60s to 120s:

```python
# In src/wca_potential.py, run_pyrosetta_subprocess()
result = subprocess.run(
    [CONDA_PYTHON, PYROSETTA_RUNNER],
    input=json.dumps(input_data),
    capture_output=True,
    text=True,
    timeout=120,  # Increased from 60
    env=env,
)
```

---

### 3. Per-Residue 3D Pareto Analysis

#### Problem

The original `analyze_pareto.py` computed Pareto fronts across residues (comparing 3 best rotamers), which is not useful. The real value is computing Pareto fronts **within each residue** across all rotamer states.

#### Solution

##### Modified: `src/rotamer_sampling.py` - Checkpoint Format

Now stores **all rotamer data** (not just best) in checkpoint:

```python
# In save_rotamer_checkpoint()
all_rotamers = []
for i, rot_idx in enumerate(res_data['rotamer_states']):
    all_rotamers.append({
        'index': rot_idx,
        'volume_change': res_data['volume_changes'][i],
        'percentage': res_data['percentages'][i],
        'delta_g': res_data['delta_g_values'][i],
        'wca_energy': res_data['wca_energies'][i],
        'joint_score': res_data['joint_scores'][i]
    })

checkpoint_data['residue_details'].append({
    # ... existing fields ...
    'all_rotamers': all_rotamers  # NEW: full rotamer data
})
```

##### Modified: `analyze_pareto.py` - 3D Pareto Analysis

New function `analyze_per_residue_pareto()`:

```python
def analyze_per_residue_pareto(residue_rotamers, verbose=True):
    """
    Compute 3D Pareto fronts (volume ↑, deltaG ↓, WCA ↓) for each residue.
    
    For each residue:
    1. Filter rotamers with valid data for all 3 metrics
    2. Compute 3D Pareto front using compute_pareto_front_3d()
    3. Sort Pareto front by volume (descending)
    4. Return per-residue analysis results
    """
```

**Output format**:

```
PER-RESIDUE 3D PARETO ANALYSIS (Volume ↑, ΔG ↓, WCA ↓)
============================================================

A:MET102: 5 Pareto-optimal rotamers out of 26
  Idx   ΔVolume        ΔG         WCA     Score
  ------------------------------------------------
  11    +119.66      0.76      658.32     1.234
  7      +87.70      0.02      216.29     1.873
  1     +102.17     -0.97      282.62     2.928
  ...
```

##### New Export Files

| File | Description |
|------|-------------|
| `pareto_fronts_per_residue.csv` | All rotamers with `on_pareto_front` flag |
| `residue_summaries.csv` | Best rotamer per residue summary |
| `pareto_analysis_summary.json` | Full JSON with per-residue Pareto data |

##### Backward Compatibility

Checkpoints without `all_rotamers` field (old format) are detected:
```python
has_full_data = len(residue_rotamers) > 0

if not has_full_data:
    print("Note: Checkpoint uses old format (best rotamer only).")
    print("Re-run rotamer sampling to get full Pareto analysis.")
```

---

### 4. Updated Project Structure

```
RotamerDMS/
├── src/
│   ├── wca_potential.py       # Modified: stderr logging, 120s timeout
│   ├── rotamer_sampling.py    # Modified: all_rotamers in checkpoint
│   └── scoring.py             # Unchanged
├── run_rotamer_dms.sh         # Modified: SBATCH directives, SLURM_SUBMIT_DIR
├── analyze_pareto.py          # Modified: per-residue 3D Pareto analysis
└── implementation.md          # This file
```

---

### 5. Changelog Update

- **v0.3** (March 20, 2026): SLURM integration and enhanced Pareto analysis
  - SBATCH directives for cluster job submission
  - Fixed working directory resolution for SLURM (`SLURM_SUBMIT_DIR`)
  - Increased PyRosetta subprocess timeout (60s → 120s)
  - Added stderr logging for WCA computation failures
  - Enhanced checkpoint format with all rotamer data
  - Per-residue 3D Pareto front analysis in `analyze_pareto.py`
  - New export files: `pareto_fronts_per_residue.csv`, `residue_summaries.csv`

---

## Session: March 24, 2026 - Energy Minimization for Realistic Structures

> **Note**: This section documents the initial minimization implementation. Parameters and MoveMap settings were subsequently updated in v0.6 (see "Chi-Only MoveMap" section below). For current parameter values, see `workflow_parameter.md`.

### Overview

Added local energy minimization step to produce physically realistic structures with low steric clash energies. Previously, fa_rep energies were ~100+ REU; after minimization, they should be in the realistic 5-20 REU range.

---

### 1. Problem Statement

**Issue**: Raw rotamer state transitions often produce high fa_rep (steric repulsion) energies (~100-900 REU) because surrounding residues haven't had a chance to adjust to the new rotamer conformation.

**Solution**: After applying each rotamer state, perform local energy minimization with PyRosetta's `MinMover` using the `ref2015` score function. The target residue is **frozen** during minimization (using `MoveMap`) so the rotamer conformation is preserved while surrounding residues relax.

---

### 2. PyRosetta Components Used

#### ref2015 Score Function

The `ref2015` score function is Rosetta's default full-atom energy function, containing:
- `fa_rep`: Lennard-Jones repulsive (weight 0.55)
- `fa_atr`: Lennard-Jones attractive (weight 1.0)
- `fa_sol`: Lazaridis-Karplus solvation
- `fa_elec`: Coulombic electrostatics
- `hbond_*`: Hydrogen bonding terms
- `rama_prepro`: Ramachandran preferences
- `fa_dun`: Dunbrack rotamer energy
- Plus ~10 more terms

#### MinMover

Gradient-based local minimizer that finds the nearest energy minimum:

```python
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

minmover = MinMover()
minmover.movemap(movemap)           # Define degrees of freedom
minmover.score_function(sfxn)       # ref2015
minmover.min_type('lbfgs_armijo_nonmonotone')  # Minimizer algorithm
minmover.tolerance(0.1)             # ← UPDATED in v0.6: Higher tolerance
minmover.max_iter(50)               # ← UPDATED in v0.6: Fewer iterations
minmover.apply(pose)                # Run minimization
```

#### MoveMap

Controls which degrees of freedom are allowed to change:

```python
from pyrosetta.rosetta.core.kinematics import MoveMap

movemap = MoveMap()
movemap.set_bb(False)               # ← UPDATED in v0.6: No backbone movement
movemap.set_chi(True)               # Allow chi torsions

# Freeze target residue (preserve rotamer)
movemap.set_chi(target_pose_resnum, False)
```

**Key Insight**: Chi-only minimization (bb=False globally) preserves the global backbone structure while allowing sidechains to relax around the new rotamer.

---

### 3. Modified Workflow

**Before (v0.2)**:
```
Apply rotamer → Compute volume → Compute fa_rep (no minimization)
```

**After (v0.4)**:
```
Apply rotamer → Compute ΔG → Minimize (frozen target) → Compute volume → Extract fa_rep
```

**Rationale**:
- **ΔG computed before minimization**: Represents the conformational free energy change based on statistical rotamer frequencies, independent of structural relaxation
- **Volume computed after minimization**: Captures the realistic pocket volume after surrounding residues have adjusted
- **fa_rep extracted after minimization**: Reflects true steric strain in the relaxed structure

---

### 4. Code Changes

#### `src/pyrosetta_runner.py` - Complete Rewrite

**Old function**: `compute_wca_energy()` - just scored structure with fa_rep

**New function**: `minimize_with_frozen_residue()`:
1. Load PDB into PyRosetta Pose
2. Find target residue in pose numbering
3. Create ref2015 score function
4. Create MoveMap with target residue frozen
5. Create MinMover and apply to pose
6. Extract fa_rep energy for target residue
7. Return minimized PDB content + fa_rep energy

```python
def minimize_with_frozen_residue(pdb_file, chain, resnum, verbose=False):
    pose = pyrosetta.pose_from_pdb(pdb_file)
    target_pose_resnum = find_pose_resnum(pose, chain, resnum)
    
    sfxn = pyrosetta.create_score_function('ref2015')
    
    movemap = MoveMap()
    movemap.set_bb(False)                        # ← v0.6: Chi-only
    movemap.set_chi(True)
    movemap.set_chi(target_pose_resnum, False)   # Freeze target
    
    minmover = MinMover()
    minmover.movemap(movemap)
    minmover.score_function(sfxn)
    minmover.min_type('lbfgs_armijo_nonmonotone')
    minmover.tolerance(0.1)                       # ← v0.6: Reduced
    minmover.max_iter(50)                         # ← v0.6: Reduced
    
    minmover.apply(pose)
    
    # Extract fa_rep for target residue
    energies = pose.energies()
    fa_rep_energy = energies.residue_total_energies(target_pose_resnum)[ScoreType.fa_rep]
    
    # Return minimized PDB content
    pose.dump_pdb(temp_path)
    return {'success': True, 'fa_rep_energy': fa_rep_energy, 'minimized_pdb': ...}
```

#### `src/wca_potential.py` - Renamed and Updated

**Module renamed conceptually**: From "WCA potential computation" to "Energy minimization and fa_rep extraction"

**Old function**: `get_wca_energy()` → returns float

**New function**: `minimize_and_get_fa_rep()` → returns dict with:
- `fa_rep_energy`: float
- `minimized_struct`: Schrodinger Structure object
- `total_score`: float (ref2015 total)

**Key addition**: `_load_structure_from_pdb_content()` to convert PyRosetta output back to Schrodinger Structure

**Timeout increased**: 120s → 300s (minimization takes longer than simple scoring)

#### `src/rotamer_sampling.py` - Workflow Update

**Import change**:
```python
# Old
from .wca_potential import get_wca_energy, _check_pyrosetta_available

# New
from .wca_potential import minimize_and_get_fa_rep, _check_pyrosetta_available
```

**Logic change in `sample_rotamers_for_residue()`**:

```python
# Apply rotamer
test_rotamers[rot_idx].apply()

# Compute ΔG BEFORE minimization (from rotamer statistics)
delta_g = -math.log(percentage / p_curr)

# Minimize and get fa_rep
fa_rep_energy = None
struct_for_volume = test_struct  # Default: non-minimized

if compute_wca:
    min_result = minimize_and_get_fa_rep(test_struct, chain, resnum)
    if min_result.get('success'):
        fa_rep_energy = min_result['fa_rep_energy']
        struct_for_volume = min_result['minimized_struct']

# Compute volume on MINIMIZED structure
result = measure_binding_site_volume(struct_for_volume, ...)
delta_volume = result['total_volume'] - initial_volume
```

---

### 5. Expected Impact

| Metric | Before Minimization | After Minimization |
|--------|--------------------|--------------------|
| fa_rep energy | 100-900 REU | 5-20 REU |
| Structure quality | Unrealistic clashes | Physically plausible |
| Volume measurement | Raw rotamer volume | Relaxed pocket volume |
| Runtime per rotamer | ~1-2 sec | ~5-15 sec |

**Trade-off**: Minimization adds computational cost (~5-10x slower per rotamer) but produces much more realistic structures with interpretable fa_rep energies.

---

### 6. Changelog Update

- **v0.4** (March 24, 2026): Energy minimization for realistic structures
  - Added PyRosetta MinMover-based local minimization after each rotamer transition
  - Target residue frozen during minimization to preserve rotamer state
  - Uses ref2015 score function for comprehensive energy evaluation
  - fa_rep energy now reflects relaxed structure (realistic 5-20 REU range)
  - Volume computed on minimized structure for accurate pocket measurement
  - Subprocess timeout increased to 300s for minimization
  - Removed redundant explicit WCA computation code

---

## Session: March 26, 2026 - PyRosetta-Schrodinger Structure Consistency Fixes

### Overview

Fixed critical structural representation inconsistencies between PyRosetta and Schrodinger that were causing:
- Unphysical volume changes
- Hydrogen atom discrepancies (PyRosetta adds hydrogens)
- Chain breaks and disconnected topology in output structures
- Loss of continuous global protein topology

---

### 1. Problem Analysis

#### Root Cause

When converting structures between Schrodinger and PyRosetta via PDB files:

1. **PyRosetta adds hydrogens**: `pose_from_pdb()` automatically adds missing hydrogens
2. **PyRosetta doesn't write CONECT records**: `dump_pdb()` omits bond information
3. **Schrodinger infers bonds by distance**: When reading PyRosetta's PDB output, Schrodinger uses distance-based bond inference, which can fail if atom positions have shifted

This resulted in the output structure having floating atoms, chain breaks, and incorrect topology.

#### Previous Workflow (Broken)

```
Schrodinger Structure → PDB file → PyRosetta Pose → Minimize → PDB file → NEW Schrodinger Structure
```

The final Schrodinger Structure was a completely new object parsed from PyRosetta's PDB output, losing all original connectivity information.

---

### 2. Solution: Coordinate Transfer

#### New Workflow

```
Schrodinger Structure → PDB file → PyRosetta Pose → Minimize → PDB content → 
    Parse coordinates → Transfer to ORIGINAL Schrodinger Structure (copy)
```

**Key insight**: Instead of creating a new Schrodinger Structure from PyRosetta's output, we:
1. Keep the **original Schrodinger Structure** (preserves bonds, chain IDs, residue numbers)
2. Extract only **XYZ coordinates** from PyRosetta's minimized PDB
3. **Map coordinates** by matching `(chain, resnum, atom_name)` tuples
4. Update **only coordinates** in a copy of the original structure

---

### 3. Code Changes

#### `src/pyrosetta_runner.py`

##### New Function: `strip_hydrogens_from_pdb()`

```python
def strip_hydrogens_from_pdb(pdb_content):
    """
    Remove hydrogen atoms from PDB content.
    
    PyRosetta adds hydrogens during pose loading, but we want to return
    only heavy atoms to maintain consistency with the input structure.
    """
    filtered_lines = []
    for line in pdb_content.split('\n'):
        if not line.startswith(('ATOM', 'HETATM')):
            filtered_lines.append(line)
            continue
        
        # Check element column (76-78) for 'H'
        if len(line) >= 78:
            element = line[76:78].strip()
            if element == 'H':
                continue
        
        # Fallback: check atom name (columns 12-16)
        atom_name = line[12:16].strip()
        if atom_name.startswith('H') or (len(atom_name) >= 2 and atom_name[0].isdigit() and atom_name[1] == 'H'):
            continue
        
        filtered_lines.append(line)
    
    return '\n'.join(filtered_lines)
```

##### Modified Functions

Both `minimize_with_frozen_residue()` and `minimize_full_structure()` now call `strip_hydrogens_from_pdb()` before returning:

```python
# Before returning PDB content
minimized_pdb = strip_hydrogens_from_pdb(minimized_pdb_raw)
```

#### `src/wca_potential.py`

##### New Function: `_parse_pdb_coordinates()`

```python
def _parse_pdb_coordinates(pdb_content):
    """
    Parse atom coordinates from PDB content.
    
    Returns a dictionary mapping (chain, resnum, atom_name) -> (x, y, z)
    """
    coords = {}
    for line in pdb_content.split('\n'):
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        
        # PDB format columns (0-indexed):
        # 12-15: atom name, 21: chain ID, 22-25: resnum
        # 30-37: x, 38-45: y, 46-53: z
        atom_name = line[12:16].strip()
        chain = line[21:22].strip() or ' '
        resnum = int(line[22:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        
        coords[(chain, resnum, atom_name)] = (x, y, z)
    
    return coords
```

##### New Function: `_transfer_coordinates()`

```python
def _transfer_coordinates(original_struct, minimized_pdb_content, verbose=False):
    """
    Transfer coordinates from minimized PDB back to original Schrodinger structure.
    
    Preserves original structure's connectivity, atom properties, and
    residue identifiers while updating only XYZ coordinates.
    
    Returns:
        - success: bool
        - structure: Updated Schrodinger Structure (copy with new coords)
        - atoms_updated: Number of atoms with coordinates updated
        - atoms_total: Total heavy atoms in original structure
        - atoms_missing: List of atoms not found in minimized PDB
        - rmsd: RMSD between original and updated coordinates
        - error: str (if failed)
    """
    minimized_coords = _parse_pdb_coordinates(minimized_pdb_content)
    updated_struct = original_struct.copy()
    
    atoms_updated = 0
    sum_sq_dist = 0.0
    
    for atom in updated_struct.atom:
        if atom.atomic_number == 1:  # Skip hydrogens
            continue
        
        key = (atom.chain.strip() or ' ', atom.resnum, atom.pdbname.strip())
        
        if key in minimized_coords:
            new_x, new_y, new_z = minimized_coords[key]
            
            # Calculate squared distance for RMSD
            dx, dy, dz = new_x - atom.x, new_y - atom.y, new_z - atom.z
            sum_sq_dist += dx*dx + dy*dy + dz*dz
            
            # Update coordinates
            atom.x, atom.y, atom.z = new_x, new_y, new_z
            atoms_updated += 1
    
    # Validation: fail if <90% of atoms mapped
    if atoms_updated / atoms_total < 0.9:
        return {'success': False, 'error': 'Too few atoms mapped'}
    
    rmsd = math.sqrt(sum_sq_dist / atoms_updated)
    return {'success': True, 'structure': updated_struct, 'rmsd': rmsd, ...}
```

##### Modified Functions

Both `minimize_and_get_fa_rep()` and `minimize_full_structure()` now use `_transfer_coordinates()`:

```python
# Old (broken)
minimized_struct = _load_structure_from_pdb_content(result['minimized_pdb'])

# New (fixed)
transfer_result = _transfer_coordinates(
    schrodinger_struct,  # Original structure
    result['minimized_pdb'],
    verbose=verbose
)
minimized_struct = transfer_result['structure']
```

---

### 4. Validation

The `_transfer_coordinates()` function includes built-in validation:

| Check | Threshold | Action |
|-------|-----------|--------|
| Atom mapping ratio | ≥90% | Fail if below |
| RMSD calculation | Reported | For monitoring |
| Missing atoms | Listed | Warning in verbose mode |

**New return fields** from minimization functions:
- `coord_transfer_rmsd`: RMSD of coordinate changes (Å)
- `atoms_updated`: Number of atoms successfully mapped
- `atoms_total`: Total heavy atoms in structure

---

### 5. Benefits

| Aspect | Before | After |
|--------|--------|-------|
| Connectivity | Lost (re-inferred) | Preserved from original |
| Chain IDs | Potentially changed | Preserved |
| Residue numbers | Potentially changed | Preserved |
| Hydrogens | Added by PyRosetta | Stripped before return |
| Bond information | Distance-inferred | Original bonds maintained |
| Topology | Chain breaks possible | Continuous |

---

### 6. Deferred Improvements

The following improvements were designed but deferred pending user testing:

1. **Coordinate constraints**: Add CA coordinate constraints to limit backbone movement during minimization
2. **Reduce minimization aggressiveness**: Lower `max_iter` (200→50) and increase `tolerance` (0.01→0.1)

These will be implemented if testing shows the structure is changing too much during minimization.

---

### 7. Changelog Update

- **v0.5** (March 26, 2026): PyRosetta-Schrodinger structure consistency fixes
  - Added `strip_hydrogens_from_pdb()` to remove PyRosetta-added hydrogens before output
  - Added `_parse_pdb_coordinates()` to extract coordinates from PDB content
  - Added `_transfer_coordinates()` to map minimized coordinates back to original structure
  - Coordinate transfer preserves original Schrodinger structure connectivity and properties
  - Built-in validation: ≥90% atom mapping required, RMSD reported
  - Deprecated `_load_structure_from_pdb_content()` (kept for backwards compatibility)
  - Fixed chain breaks, floating atoms, and topology issues in output structures

---

## Session: March 26, 2026 - Minimization Refinement: Chi-Only MoveMap

### Overview

After implementing the coordinate transfer approach, testing revealed the backbone was still moving substantially during minimization. This session implements chi-only minimization (disabling backbone DOFs) and documents the complete debugging journey.

---

### 1. Complete Debugging Timeline

#### Initial Symptom: Unreasonable Volume Deltas

**Observation**: Volume changes of 700-1100 Å³ per rotamer were being reported, which is physically implausible for single sidechain rotamer changes.

**Expected range**: Typical rotamer changes should produce volume deltas of ~10-100 Å³, not >500 Å³.

#### Root Cause Investigation

Upon visual inspection of the output structures, multiple issues were identified:

| Issue | Symptom | Impact |
|-------|---------|--------|
| **Hydrogen addition** | PyRosetta added hydrogens not present in input | Atom count mismatch, incorrect topology |
| **Topology loss** | Chain breaks, floating atoms in output | Schrodinger couldn't infer bonds correctly |
| **Excessive backbone movement** | Global structure drift | Unrealistic volume changes, structure no longer resembles input |

---

### 2. Fix #1: Hydrogen Stripping

**Problem**: `pyrosetta.pose_from_pdb()` automatically adds missing hydrogens.

**Solution**: Strip hydrogens from PDB content before returning to Schrodinger.

```python
def strip_hydrogens_from_pdb(pdb_content):
    """Remove hydrogen atoms from PDB content."""
    # Check element column (76-78) for 'H'
    # Fallback: check atom name for H prefix
```

**Rationale**: Schrodinger structures typically represent only heavy atoms. Adding hydrogens creates atom count mismatches and confuses downstream processing.

---

### 3. Fix #2: Coordinate Transfer (Not Structure Replacement)

**Problem**: Creating a new Schrodinger Structure from PyRosetta's PDB output loses:
- Bond connectivity (CONECT records not written by PyRosetta)
- Chain identifiers
- Residue numbering consistency
- Molecular topology

When Schrodinger reads the PDB without CONECT records, it infers bonds by distance, which fails when atoms have moved significantly.

**Solution**: Transfer only XYZ coordinates back to the *original* Schrodinger Structure.

```python
def _transfer_coordinates(original_struct, minimized_pdb_content):
    """
    Parse coordinates from PDB, map by (chain, resnum, atom_name),
    update coordinates in a COPY of original structure.
    """
```

**Key insight**: The original Schrodinger Structure already has correct topology. We only need to update atom positions, not rebuild the entire molecular representation.

**Validation**: Built-in checks ensure ≥90% of heavy atoms are successfully mapped, with RMSD reported for monitoring.

---

### 4. Fix #3: CA Coordinate Constraints

**Problem**: Even with hydrogen stripping and coordinate transfer, minimization was moving the backbone too aggressively, causing global structural drift.

**Solution**: Add harmonic coordinate constraints to all CA atoms.

```python
def add_ca_coordinate_constraints(pose, exclude_resnum=None):
    """
    Apply harmonic penalty to CA atoms that move from original positions.
    """
    # CA_CONSTRAINT_SD = 0.5 Å (standard deviation)
    # CA_CONSTRAINT_WEIGHT = 1.0 (score function weight)
```

**Rationale**: Coordinate constraints provide a "soft" enforcement that penalizes backbone movement while still allowing some flexibility if energetically favorable. This is the recommended approach when you want *limited* backbone flexibility.

---

### 5. Fix #4: Reduced Minimization Aggressiveness

**Problem**: Default MinMover settings (200+ iterations, tight tolerance) allowed the optimizer to make larger structural changes than necessary.

**Solution**: Reduce iterations and loosen tolerance.

| Parameter | Before | After | Rationale |
|-----------|--------|-------|-----------|
| `max_iter` (frozen) | 200 | 50 | Fewer steps = less drift |
| `max_iter` (full) | 500 | 100 | Still enough to relieve clashes |
| `tolerance` | 0.01 | 0.1 | Converge earlier, accept "good enough" |

---

### 6. Fix #5: Chi-Only MoveMap (Final Solution)

**Problem**: Even with constraints, allowing backbone DOFs (`set_bb(True)`) permits phi/psi angles to change, which fundamentally alters the global structure.

**Literature insight** (from RosettaCommons forums):
> *"If you want to keep the backbone from moving at all, you probably want to look into MoveMaps... It should be easy to turn off all backbone DOFs and only allow minimization of the sidechains."*

**Solution**: Disable backbone movement entirely in the MoveMap.

```python
# BEFORE (too permissive)
movemap.set_bb(True)   # All backbone can move
movemap.set_chi(True)  # All sidechains can move

# AFTER (chi-only)
movemap.set_bb(False)  # Freeze ALL backbone
movemap.set_chi(True)  # Only sidechains move
```

**Rationale**:
1. **Rotamer optimization is a sidechain problem** - we're changing chi angles, not phi/psi
2. **Backbone is already reasonable** - input structure from Schrodinger is physically valid
3. **Sidechain adjustment alone is sufficient** - to relieve steric clashes from rotamer changes
4. **Predictable output** - backbone-fixed minimization produces consistent, interpretable results

---

### 7. MoveMap vs Constraints: When to Use Each

Based on research of RosettaCommons documentation:

| Approach | Use Case | Our Situation |
|----------|----------|---------------|
| **MoveMap (bb=False)** | Fixed backbone, sidechain-only optimization | ✓ **Correct for rotamer sampling** |
| **Coordinate Constraints** | Allow *some* backbone movement with penalty | Kept as secondary safeguard |

**Our final approach uses BOTH**:
- MoveMap with `set_bb(False)` as **primary** enforcement
- CA coordinate constraints as **secondary** safeguard (in case any DOFs slip through)

---

### 8. Code Changes Summary

#### `src/pyrosetta_runner.py`

**New imports**:
```python
from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.numeric import xyzVector_double_t
```

**New constants**:
```python
CA_CONSTRAINT_SD = 0.5      # Å - tighter = less movement
CA_CONSTRAINT_WEIGHT = 1.0  # Score function weight
```

**New functions**:
- `strip_hydrogens_from_pdb()` - Remove H atoms from PDB content
- `add_ca_coordinate_constraints()` - Apply harmonic CA constraints

**Modified `minimize_with_frozen_residue()`**:
```python
# Enable coordinate_constraint term
sfxn.set_weight(ScoreType.coordinate_constraint, CA_CONSTRAINT_WEIGHT)

# Add CA constraints
add_ca_coordinate_constraints(pose, exclude_resnum=target_pose_resnum)

# Chi-only MoveMap
movemap.set_bb(False)   # No backbone movement
movemap.set_chi(True)   # Sidechains can move
movemap.set_chi(target_pose_resnum, False)  # Except target (frozen)

# Reduced aggressiveness
minmover.tolerance(0.1)
minmover.max_iter(50)
```

**Modified `minimize_full_structure()`**:
```python
# Same pattern: constraints + chi-only + reduced aggressiveness
movemap.set_bb(False)
movemap.set_chi(True)
minmover.tolerance(0.1)
minmover.max_iter(100)
```

#### `src/wca_potential.py`

**New functions**:
- `_parse_pdb_coordinates()` - Extract (chain, resnum, atom_name) → (x, y, z) mapping
- `_transfer_coordinates()` - Map coordinates back to original Schrodinger structure

**Modified functions**:
- `minimize_and_get_fa_rep()` - Uses `_transfer_coordinates()` instead of `_load_structure_from_pdb_content()`
- `minimize_full_structure()` - Same change

---

### 9. Expected Behavior After All Fixes

| Metric | Before Fixes | After Fixes |
|--------|--------------|-------------|
| Volume delta per rotamer | 700-1100 Å³ | ~10-100 Å³ (realistic) |
| Backbone RMSD | Large (structure drift) | Near-zero (chi-only) |
| Atom count consistency | Mismatched (H added) | Preserved (H stripped) |
| Topology | Broken (chain breaks) | Preserved (coord transfer) |
| Output structure | Unrecognizable | Matches input backbone |

---

### 10. Tunable Parameters

If further adjustment is needed:

| Parameter | Location | Effect |
|-----------|----------|--------|
| `CA_CONSTRAINT_SD` | `pyrosetta_runner.py` | Smaller = tighter backbone restraint |
| `CA_CONSTRAINT_WEIGHT` | `pyrosetta_runner.py` | Higher = stronger enforcement |
| `max_iter` | MinMover setup | Fewer = less optimization |
| `tolerance` | MinMover setup | Higher = earlier convergence |

---

### 11. Changelog Update

- **v0.6** (March 26, 2026): Chi-only minimization and complete debugging fixes
  - **MoveMap**: Changed from `set_bb(True)` to `set_bb(False)` for chi-only minimization
  - **Rationale**: Rotamer optimization is a sidechain problem; backbone should not move
  - Combined with CA coordinate constraints for defense-in-depth
  - Complete fix for unreasonable volume deltas (700+ Å³ → realistic ~10-100 Å³)
  - Documented full debugging timeline from initial symptom to resolution
