# RotamerDMS Workflow & Parameters

## Workflow Summary

```
1. LOAD STRUCTURES
   Reference (with ligand) + Sample (to optimize)
   ↓
2. IDENTIFY BINDING SITE
   Residues within distance_threshold of ligand
   ↓
3. MUTATION PHASE
   For each binding site residue:
     - Mutate to Alanine
     - Measure cavity volume change
     - Keep if volume increased
   Sort by volume impact (descending)
   ↓
4. SELECT TOP RESIDUES
   Keep top_percent of volume-affecting residues
   ↓
5. ROTAMER SAMPLING PHASE
   For each selected residue (in priority order):
     a. Load backbone-dependent rotamer library
     b. Identify current/native rotamer (by chi angle matching)
     c. For each candidate rotamer:
        - Apply rotamer to structure
        - Compute ΔG = -ln(p_suggested / p_native)
        - Minimize surrounding sidechains (chi-only, target frozen)
        - Measure cavity volume
        - Extract fa_rep energy
     d. Rank rotamers by joint score:
        score = w_vol*z(ΔV) - w_dg*z(ΔG) - w_wca*z(fa_rep)
     e. Select best rotamer, apply permanently
   ↓
6. FINAL MINIMIZATION
   Chi-only minimization of full structure
   ↓
7. OUTPUT
   Optimized PDB structure + checkpoint JSON
```

---

## All Parameters

### CLI Parameters (User-Configurable)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--distance` | 7.0 Å | Ligand-residue distance threshold for binding site |
| `--top-percent` | 100 | % of volume-affecting residues to sample |
| `--min-contacts` | 4 | Min binding site contacts for cavity counting |
| `--w-vol` | 1.0 | Joint score weight for volume |
| `--w-dg` | 1.0 | Joint score weight for ΔG |
| `--w-wca` | 1.0 | Joint score weight for fa_rep |
| `--wca-threshold` | None | Hard fa_rep cutoff (REU) to reject rotamers |
| `--no-wca` | False | Disable minimization/fa_rep computation |

---

### Cavity Detection (pyKVFinder)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `step` | 0.6 Å | Grid spacing |
| `probe_in` | 1.4 Å | Inner probe radius |
| `probe_out` | 4.0 Å | Outer probe radius |
| `removal_distance` | 2.4 Å | Exterior trim distance |
| `volume_cutoff` | 5.0 Å³ | Minimum cavity volume |

*Source: ChimeraX "Find Cavities" defaults*

---

### PyRosetta Minimization

#### Score Function
| Parameter | Value |
|-----------|-------|
| Score function | `ref2015` |
| `coordinate_constraint` weight | 1.0 |

#### CA Coordinate Constraints
| Parameter | Value | Description |
|-----------|-------|-------------|
| `CA_CONSTRAINT_SD` | 0.5 Å | Harmonic std dev (tighter = less movement) |
| `CA_CONSTRAINT_WEIGHT` | 1.0 | Score function weight |

#### MoveMap (Degrees of Freedom)
| DOF | Setting | Description |
|-----|---------|-------------|
| Backbone (phi/psi) | **False** | No backbone movement |
| Chi angles | **True** | Sidechain-only minimization |
| Target residue chi | **False** | Preserve applied rotamer |

#### MinMover Settings
| Parameter | Frozen Residue | Full Structure |
|-----------|----------------|----------------|
| `min_type` | `lbfgs_armijo_nonmonotone` | `lbfgs_armijo_nonmonotone` |
| `tolerance` | 0.1 | 0.1 |
| `max_iter` | 50 | 100 |

#### Subprocess
| Parameter | Value |
|-----------|-------|
| Timeout | 300 seconds |

---

### Structure Transfer (PyRosetta → Schrodinger)

| Parameter | Value | Description |
|-----------|-------|-------------|
| Hydrogen handling | Strip all | PyRosetta adds H; we remove before transfer |
| Coordinate mapping | (chain, resnum, atom_name) | Key for atom matching |
| Min mapping ratio | 90% | Fail if <90% atoms mapped |
| RMSD | Computed | Reported for monitoring |

---

### Rotamer Library

| Parameter | Value |
|-----------|-------|
| Library type | Backbone-dependent (Schrodinger) |
| Native ID method | Chi angle distance minimization |
| Symmetry handling | PHE/TYR χ2, ASP χ2, GLU χ3 (180° equiv) |

---

### Joint Scoring

| Component | Direction | Z-Score Normalized |
|-----------|-----------|-------------------|
| Volume change (ΔV) | Maximize (+) | Yes |
| Free energy (ΔG) | Minimize (-) | Yes |
| Steric clash (fa_rep) | Minimize (-) | Yes |

**Formula**: `score = w_vol*z(ΔV) - w_dg*z(ΔG) - w_wca*z(fa_rep)`

---

### Output Files

| File | Location | Description |
|------|----------|-------------|
| `<name>_optimized.pdb` | `output/` | Final structure |
| `mutation_phase_*.json` | `checkpoints/` | Ala scan results |
| `rotamer_sampling_*.json` | `checkpoints/` | Rotamer results + Pareto |

---

### SLURM Resources

| Resource | Value |
|----------|-------|
| CPUs | 4 |
| Memory | 16 GB |
| Time | 1 hour |
| Partition | `lyu_a` |
