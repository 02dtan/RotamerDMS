# Forward Outlook: Physically Rigorous Enhancements to RotamerDMS

**Author:** Generated with literature review  
**Date:** March 24, 2026  
**Purpose:** Document proposed algorithmic improvements to address the problem of over-opened binding pockets leading to loss of docking selectivity.

---

## Table of Contents

1. [Problem Statement](#1-problem-statement)
2. [Current Algorithm Limitations](#2-current-algorithm-limitations)
3. [Proposed Enhancement 1: Reference Ligand Volume Constraint](#3-proposed-enhancement-1-reference-ligand-volume-constraint)
4. [Proposed Enhancement 2: Anchor Residue Exclusion](#4-proposed-enhancement-2-anchor-residue-exclusion)
5. [Proposed Enhancement 3: Hot Spot Preservation via Computational Solvent Mapping](#5-proposed-enhancement-3-hot-spot-preservation-via-computational-solvent-mapping)
6. [Proposed Enhancement 4: MD-Derived Conformational Energy Penalties](#6-proposed-enhancement-4-md-derived-conformational-energy-penalties)
7. [Proposed Enhancement 5: Shape Complementarity Scoring](#7-proposed-enhancement-5-shape-complementarity-scoring)
8. [Proposed Enhancement 6: Desolvation Penalty](#8-proposed-enhancement-6-desolvation-penalty)
9. [Integrated Multi-Objective Function](#9-integrated-multi-objective-function)
10. [Implementation Priority](#10-implementation-priority)
11. [References](#11-references)

---

## 1. Problem Statement

### Observed Behavior

When the RotamerDMS algorithm is applied to sample structures from BioEmu, the resulting modified structures exhibit a critical failure mode during downstream molecular docking:

> **All compounds, including property-matched decoy ligands, dock favorably to the modified protein structure.**

This observation indicates that the binding pocket has been opened to such an extent that it loses **selectivity**—the ability to discriminate between true binders and non-binders.

### Root Cause Analysis

The current algorithm optimizes for **maximum cavity volume** subject to:
1. Conformational free energy penalty (ΔG from backbone-dependent rotamer library frequencies)
2. Steric clash avoidance (WCA potential / PyRosetta `fa_rep`)

However, **volume maximization is fundamentally anti-selective**. The literature on protein-ligand binding consistently demonstrates that:

1. Binding affinity correlates with **shape complementarity**, not pocket size (Lawrence & Colman, 1993; Nicholls et al., 1991)
2. Large, solvent-exposed pockets incur **desolvation penalties** that reduce binding affinity (Mysinger & Shoichet, 2010)
3. Selective binding sites contain **"hot spots"**—regions of high ligand-binding propensity that must be preserved (Zerbe et al., 2012; Vajda et al., 2018)
4. Using too many receptor conformations or overly open conformations in flexible docking **reduces enrichment** of true ligands over decoys (Fischer et al., 2014)

---

## 2. Current Algorithm Limitations

### 2.1 Volume Maximization Without Upper Bound

The current objective function rewards any increase in cavity volume:

```
Score_current = w_vol × z(ΔV) - w_dg × z(ΔG) - w_wca × z(WCA)
```

Where:
- `ΔV` = change in cavity volume (Å³)
- `ΔG` = conformational free energy change (arbitrary units)
- `WCA` = Weeks-Chandler-Andersen repulsive energy (REU)
- `z(·)` = z-score normalization

**Problem:** There is no constraint preventing the volume from growing arbitrarily large. Any rotamer state that increases volume will be favored, even if it creates an unrealistically open pocket.

### 2.2 Context-Independent Rotamer Frequencies

The current ΔG calculation uses Schrodinger's backbone-dependent rotamer library:

```python
delta_g = -math.log(p_suggested / p_native)  # Arbitrary units
```

Where `p_suggested` and `p_native` are library frequencies.

**Problem:** These frequencies are derived from statistical analysis of the PDB and represent **average** rotamer preferences for a given φ/ψ backbone configuration. They do not account for:
- The specific local environment of the target protein
- Interactions with neighboring residues
- The actual conformational landscape accessible to the protein

### 2.3 No Preservation of Binding Specificity Determinants

The algorithm treats all binding site residues equally, with the only distinction being their impact on cavity volume when mutated to alanine.

**Problem:** Some residues are **critical for binding specificity**:
- Anchor residues that make multiple contacts with ligands
- Residues forming druggable "hot spots"
- Conserved residues across homologs

Flipping these residues outward destroys the structural features that enable selective ligand recognition.

---

## 3. Proposed Enhancement 1: Reference Ligand Volume Constraint

### 3.1 Physical Basis

The binding pocket needs to be large enough to accommodate physiologically relevant ligands, but not larger. Shoichet and colleagues have demonstrated that **constraining receptor flexibility** to physically reasonable states dramatically improves docking selectivity (Fischer et al., 2014; Merski et al., 2021).

From thermodynamic principles, the optimal pocket volume `V_opt` should satisfy:

```
V_ligand < V_opt < V_ligand × α
```

Where:
- `V_ligand` = molecular volume of the reference ligand
- `α` = expansion factor (typically 1.3–1.8) accounting for:
  - Imperfect packing between ligand and pocket
  - Space for conserved water molecules
  - Breathing room for ligand flexibility

### 3.2 Mathematical Formulation

Instead of maximizing volume, we define a **Gaussian penalty** centered on a target volume:

$$V_{target} = V_{ligand} \times \alpha$$

$$P_{volume}(V) = -\frac{(V - V_{target})^2}{2\sigma^2}$$

Where:
- `V` = current cavity volume
- `V_target` = target cavity volume
- `σ` = standard deviation controlling penalty steepness (recommended: 30–50 Å³)

This penalty:
- Equals 0 when `V = V_target`
- Becomes increasingly negative as `V` deviates from `V_target`
- Penalizes both under-opened and over-opened pockets

### 3.3 Implementation

```python
import math
from schrodinger.structutils.measure import get_volume

def compute_ligand_volume(ligand_struct):
    """
    Compute molecular volume of the reference ligand.
    
    Uses Schrodinger's volume calculation based on van der Waals radii.
    Alternative: Use convex hull volume for more conservative estimate.
    
    Parameters
    ----------
    ligand_struct : schrodinger.structure.Structure
        Structure object containing only the ligand atoms.
    
    Returns
    -------
    float
        Ligand volume in Å³.
    """
    # Schrodinger approach
    volume = get_volume(ligand_struct)
    return volume


def compute_target_volume(ligand_volume, alpha=1.5):
    """
    Compute target pocket volume based on ligand size.
    
    Parameters
    ----------
    ligand_volume : float
        Molecular volume of the reference ligand (Å³).
    alpha : float
        Expansion factor. Default 1.5 based on empirical observations
        that well-packed protein-ligand interfaces have ~65-75% packing
        efficiency (Liang et al., 1998).
    
    Returns
    -------
    float
        Target pocket volume (Å³).
    """
    return ligand_volume * alpha


def volume_penalty(V_pocket, V_target, sigma=50.0):
    """
    Gaussian penalty for deviation from target volume.
    
    Parameters
    ----------
    V_pocket : float
        Current pocket volume (Å³).
    V_target : float
        Target pocket volume (Å³).
    sigma : float
        Standard deviation controlling penalty steepness (Å³).
        Recommended: 30-50 Å³ based on typical pocket volume variations.
    
    Returns
    -------
    float
        Penalty score (0 at target, negative elsewhere).
    
    Notes
    -----
    The Gaussian form is chosen because:
    1. It penalizes both under- and over-opening
    2. The quadratic form near the optimum allows smooth optimization
    3. σ provides intuitive control over acceptable volume range
    """
    return -((V_pocket - V_target) ** 2) / (2 * sigma ** 2)


def hard_volume_constraint(V_pocket, V_max):
    """
    Hard upper bound on pocket volume.
    
    Parameters
    ----------
    V_pocket : float
        Current pocket volume (Å³).
    V_max : float
        Maximum allowed volume (Å³). Recommended: 2.0 × V_ligand.
    
    Returns
    -------
    bool
        True if volume is acceptable, False if rejected.
    """
    return V_pocket <= V_max
```

### 3.4 Recommended Parameters

| Parameter | Recommended Value | Justification |
|-----------|-------------------|---------------|
| `α` (expansion factor) | 1.5 | Packing efficiency ~67%, consistent with Liang et al. (1998) |
| `σ` (penalty steepness) | 50 Å³ | Allows ~100 Å³ variation before significant penalty |
| `V_max` (hard constraint) | 2.0 × V_ligand | Prevents grossly over-opened pockets |

### 3.5 Literature Support

- **Fischer et al. (2014):** "High-energy decoy conformations may dominate the docking" if conformational states are not properly weighted. The paper demonstrates that constraining receptor states improves ligand enrichment.
- **Liang et al. (1998):** Analysis of protein-ligand interfaces shows packing efficiency of 65-75%, suggesting optimal pocket volumes are ~1.3-1.5× ligand volume.

---

## 4. Proposed Enhancement 2: Anchor Residue Exclusion

### 4.1 Physical Basis

Certain binding site residues are **critical for molecular recognition**. These "anchor residues" (Rajamani et al., 2004) or "hot spot residues" (Clackson & Wells, 1995) contribute disproportionately to binding free energy.

From alanine scanning mutagenesis studies:
- A small subset of interface residues (typically 5-10) contribute >2 kcal/mol each to binding
- These residues are often evolutionarily conserved
- Their spatial arrangement defines the binding specificity

**If rotamer sampling flips these residues outward, the binding site loses its ability to discriminate between ligands.**

### 4.2 Identification Criteria

Anchor residues can be identified through multiple approaches:

#### 4.2.1 Contact-Based Identification

Residues making multiple contacts with the crystallographic ligand are likely important for binding:

$$N_{contacts}(r) = \sum_{a \in ligand} \mathbb{1}[d(r, a) < d_{cutoff}]$$

Where:
- `r` = residue
- `a` = ligand atom
- `d(r, a)` = minimum distance from any residue atom to ligand atom
- `d_cutoff` = contact distance threshold (typically 4.0 Å)

**Criterion:** Residues with `N_contacts ≥ 3` are designated as anchors.

#### 4.2.2 Conservation-Based Identification

Evolutionarily conserved residues are often functionally important:

$$C(r) = \frac{\text{number of homologs with identical residue at position } r}{\text{total number of homologs}}$$

**Criterion:** Residues with `C(r) > 0.8` (>80% conservation) are designated as anchors.

#### 4.2.3 Energy-Based Identification (if alanine scanning data available)

$$\Delta\Delta G_{binding}(r) = \Delta G_{binding}^{A} - \Delta G_{binding}^{WT}$$

**Criterion:** Residues with `|ΔΔG| > 2.0 kcal/mol` are designated as anchors.

### 4.3 Implementation

```python
from schrodinger.structutils.measure import measure_distance
from collections import defaultdict

def count_ligand_contacts(residue, ligand_struct, cutoff=4.0):
    """
    Count number of ligand atoms within contact distance of residue.
    
    Parameters
    ----------
    residue : schrodinger.structure._Residue
        Residue object to analyze.
    ligand_struct : schrodinger.structure.Structure
        Structure containing ligand atoms.
    cutoff : float
        Contact distance cutoff in Å. Default 4.0 Å corresponds to
        van der Waals contact distance for typical atom pairs.
    
    Returns
    -------
    int
        Number of ligand atoms in contact with residue.
    """
    contacts = 0
    residue_atoms = list(residue.atom)
    
    for lig_atom in ligand_struct.atom:
        for res_atom in residue_atoms:
            dist = measure_distance(res_atom, lig_atom)
            if dist < cutoff:
                contacts += 1
                break  # Count each ligand atom only once
    
    return contacts


def identify_anchor_residues(
    ref_struct,
    ligand_struct,
    binding_site_residues,
    contact_threshold=3,
    conservation_scores=None,
    conservation_threshold=0.8,
    alanine_scan_data=None,
    ddg_threshold=2.0
):
    """
    Identify anchor residues that should be excluded from rotamer sampling.
    
    Parameters
    ----------
    ref_struct : schrodinger.structure.Structure
        Reference structure with bound ligand.
    ligand_struct : schrodinger.structure.Structure
        Ligand structure extracted from reference.
    binding_site_residues : list of tuple
        List of (chain, resnum, resname) tuples for binding site residues.
    contact_threshold : int
        Minimum ligand contacts to designate as anchor. Default 3.
    conservation_scores : dict, optional
        Mapping of (chain, resnum) -> conservation score [0, 1].
    conservation_threshold : float
        Minimum conservation to designate as anchor. Default 0.8.
    alanine_scan_data : dict, optional
        Mapping of (chain, resnum) -> ΔΔG in kcal/mol.
    ddg_threshold : float
        Minimum |ΔΔG| to designate as anchor. Default 2.0 kcal/mol.
    
    Returns
    -------
    set of tuple
        Set of (chain, resnum, resname) tuples for anchor residues.
    
    Notes
    -----
    A residue is designated as an anchor if it meets ANY of the criteria.
    This conservative approach ensures important residues are not perturbed.
    """
    anchors = set()
    
    for chain, resnum, resname in binding_site_residues:
        residue = get_residue(ref_struct, chain, resnum)
        
        # Criterion 1: Multiple ligand contacts
        contacts = count_ligand_contacts(residue, ligand_struct)
        if contacts >= contact_threshold:
            anchors.add((chain, resnum, resname))
            continue
        
        # Criterion 2: High conservation
        if conservation_scores is not None:
            score = conservation_scores.get((chain, resnum), 0)
            if score >= conservation_threshold:
                anchors.add((chain, resnum, resname))
                continue
        
        # Criterion 3: Large ΔΔG from alanine scanning
        if alanine_scan_data is not None:
            ddg = alanine_scan_data.get((chain, resnum), 0)
            if abs(ddg) >= ddg_threshold:
                anchors.add((chain, resnum, resname))
    
    return anchors


def filter_residues_for_sampling(binding_site_residues, anchor_residues):
    """
    Remove anchor residues from the list of residues to sample.
    
    Parameters
    ----------
    binding_site_residues : list of tuple
        Full list of binding site residues.
    anchor_residues : set of tuple
        Set of anchor residues to exclude.
    
    Returns
    -------
    list of tuple
        Filtered list excluding anchors.
    """
    return [r for r in binding_site_residues if r not in anchor_residues]
```

### 4.4 Literature Support

- **Rajamani et al. (2004):** Introduced the concept of "anchor residues" in protein-protein interactions. Showed that anchor side chains assume specific conformations that are evolutionarily conserved.
- **Clackson & Wells (1995):** Demonstrated through alanine scanning that binding energy is concentrated in a small number of "hot spot" residues.
- **Bogan & Thorn (1998):** Statistical analysis showing hot spot residues are enriched in Trp, Tyr, and Arg, and are often surrounded by a ring of energetically less important residues.

---

## 5. Proposed Enhancement 3: Hot Spot Preservation via Computational Solvent Mapping

### 5.1 Physical Basis

Druggable binding sites contain **hot spots**—regions with high propensity for binding small organic molecules. These hot spots are characterized by:
- Concave surface topology
- Favorable balance of hydrophobic and polar interactions
- Ability to bind multiple probe types (promiscuity)

The FTMap algorithm (Brenke et al., 2009) identifies hot spots by computationally "mapping" the protein surface with small organic probes (similar to experimental fragment screening). Regions where multiple probe types cluster are designated as **consensus sites** or hot spots.

**Key insight from Zerbe et al. (2012):** Hot spot residues (defined by alanine scanning) almost always protrude into consensus sites identified by FTMap. Therefore, **rotamer transitions that move residues away from consensus sites will disrupt druggability.**

### 5.2 Mathematical Formulation

For each residue `r` in the binding site, define the **hot spot disruption penalty**:

$$P_{hotspot}(r, \chi_{new}) = \begin{cases}
-k \cdot \Delta d & \text{if } \Delta d > 0 \\
0 & \text{otherwise}
\end{cases}$$

Where:
- `Δd = d_new - d_orig`
- `d_orig` = distance from sidechain centroid to nearest hot spot center (original rotamer)
- `d_new` = distance from sidechain centroid to nearest hot spot center (new rotamer)
- `k` = penalty coefficient (recommended: 0.5–1.0 kcal/mol/Å)

This penalizes rotamer states that move sidechain atoms **away** from hot spots.

### 5.3 Hot Spot Identification Protocol

#### Using FTMap (Recommended)

1. Submit the reference structure to [FTMap server](https://ftmap.bu.edu/)
2. Download results containing probe cluster positions
3. Identify **consensus sites** (clusters where ≥16 probe molecules bind)
4. For each binding site residue, compute distance to nearest consensus site

#### Alternative: Local Hydrophobicity Mapping

If FTMap is not available, approximate hot spots using local hydrophobicity analysis:

```python
def compute_local_hydrophobicity(struct, position, radius=5.0):
    """
    Compute local hydrophobicity index around a position.
    
    Uses Kyte-Doolittle hydrophobicity scale.
    """
    kyte_doolittle = {
        'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
        'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
        'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLU': -3.5,
        'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5
    }
    
    total_hydrophobicity = 0
    count = 0
    
    for residue in struct.residue:
        centroid = get_residue_centroid(residue)
        if distance(centroid, position) < radius:
            resname = residue.pdbres.strip()
            total_hydrophobicity += kyte_doolittle.get(resname, 0)
            count += 1
    
    return total_hydrophobicity / max(count, 1)
```

### 5.4 Implementation

```python
import numpy as np

class HotSpotPreserver:
    """
    Manages hot spot identification and preservation during rotamer sampling.
    """
    
    def __init__(self, ftmap_results=None, probe_threshold=16):
        """
        Parameters
        ----------
        ftmap_results : dict, optional
            FTMap output containing probe cluster positions.
            Format: {'consensus_sites': [{'center': (x,y,z), 'probes': n}, ...]}
        probe_threshold : int
            Minimum probes to define a consensus site. Default 16 per Vajda et al.
        """
        self.consensus_sites = []
        
        if ftmap_results is not None:
            for site in ftmap_results['consensus_sites']:
                if site['probes'] >= probe_threshold:
                    self.consensus_sites.append(np.array(site['center']))
    
    def get_sidechain_centroid(self, struct, chain, resnum):
        """
        Compute centroid of sidechain heavy atoms.
        """
        residue = get_residue(struct, chain, resnum)
        backbone_atoms = {'N', 'CA', 'C', 'O', 'H', 'HA'}
        
        coords = []
        for atom in residue.atom:
            if atom.pdbname.strip() not in backbone_atoms:
                coords.append([atom.x, atom.y, atom.z])
        
        if not coords:
            return None
        
        return np.mean(coords, axis=0)
    
    def distance_to_nearest_hotspot(self, centroid):
        """
        Compute distance from centroid to nearest consensus site.
        """
        if centroid is None or not self.consensus_sites:
            return None
        
        distances = [np.linalg.norm(centroid - hs) for hs in self.consensus_sites]
        return min(distances)
    
    def compute_hotspot_penalty(
        self,
        struct_orig,
        struct_new,
        chain,
        resnum,
        k=0.5
    ):
        """
        Compute penalty for moving sidechain away from hot spots.
        
        Parameters
        ----------
        struct_orig : Structure
            Structure with original rotamer state.
        struct_new : Structure  
            Structure with proposed new rotamer state.
        chain : str
            Chain identifier.
        resnum : int
            Residue number.
        k : float
            Penalty coefficient in kcal/mol/Å. Default 0.5.
            Based on typical van der Waals interaction energies.
        
        Returns
        -------
        float
            Penalty (0 or negative). Units: kcal/mol.
        """
        if not self.consensus_sites:
            return 0.0
        
        centroid_orig = self.get_sidechain_centroid(struct_orig, chain, resnum)
        centroid_new = self.get_sidechain_centroid(struct_new, chain, resnum)
        
        if centroid_orig is None or centroid_new is None:
            return 0.0
        
        d_orig = self.distance_to_nearest_hotspot(centroid_orig)
        d_new = self.distance_to_nearest_hotspot(centroid_new)
        
        delta_d = d_new - d_orig
        
        if delta_d > 0:  # Moving away from hot spot
            return -k * delta_d
        
        return 0.0
```

### 5.5 Literature Support

- **Brenke et al. (2009):** Original FTMap paper demonstrating computational solvent mapping reproduces experimental fragment screening results.
- **Zerbe et al. (2012):** "The residues protruding into hot spot regions identified by computational mapping are almost always themselves hot spot residues as defined by alanine scanning experiments."
- **Vajda et al. (2018):** "A druggable site can be defined as a region that comprises a main hot spot binding at least 16 probe clusters."

---

## 6. Proposed Enhancement 4: MD-Derived Conformational Energy Penalties

### 6.1 Physical Basis

The current ΔG calculation uses backbone-dependent rotamer library frequencies:

$$\Delta G_{library} = -k_BT \ln\left(\frac{p_{suggested}}{p_{native}}\right)$$

These frequencies represent **average** preferences across the PDB and do not account for the specific conformational landscape of the target protein.

Fischer et al. (2014) and Merski et al. (2021) demonstrated that using **experimentally-derived or MD-derived conformational populations** dramatically improves docking selectivity by preventing domination of high-energy states.

From the Merski et al. study on T4 lysozyme L99A:

> "Protein flexibility remains a major challenge in library docking because of difficulties in sampling conformational ensembles with accurate probabilities. Here, we use the model cavity site of T4 lysozyme L99A to test flexible receptor docking with energy penalties from molecular dynamics (MD) simulations."

### 6.2 Mathematical Formulation

Run MD simulations on the sample structure and compute **empirical rotamer populations** for each binding site residue:

$$p_i^{MD} = \frac{n_i}{\sum_j n_j}$$

Where:
- `n_i` = number of MD frames where residue adopts rotamer state `i`

The MD-derived conformational penalty is:

$$\Delta G_{MD} = -k_BT \ln\left(\frac{p_{suggested}^{MD}}{p_{native}^{MD}}\right)$$

Where `k_BT ≈ 0.6 kcal/mol` at 300K.

### 6.3 Rotamer Assignment from MD Trajectories

To assign rotamer states from MD, compare χ dihedral angles to library rotamers:

$$\text{rotamer}_i = \arg\min_j \sum_k (1 - \cos(\chi_k^{frame} - \chi_k^{lib,j}))$$

This handles angle periodicity correctly.

### 6.4 Implementation

```python
import math
import numpy as np
from collections import Counter

# kT at 300K
KT_300K = 0.001987 * 300  # kcal/mol

def assign_rotamer_from_chi_angles(chi_angles, rotamer_library):
    """
    Assign frame to closest rotamer state in library.
    
    Parameters
    ----------
    chi_angles : list of float
        Chi dihedral angles from MD frame (degrees).
    rotamer_library : list of dict
        Library rotamers with 'chi_angles' key.
    
    Returns
    -------
    int
        Index of closest library rotamer.
    """
    def angular_distance(angles1, angles2):
        """Compute angular distance accounting for periodicity."""
        total = 0
        for a1, a2 in zip(angles1, angles2):
            diff = abs(a1 - a2)
            diff = min(diff, 360 - diff)  # Periodicity
            total += diff ** 2
        return math.sqrt(total)
    
    min_dist = float('inf')
    best_idx = 0
    
    for idx, rot in enumerate(rotamer_library):
        dist = angular_distance(chi_angles, rot['chi_angles'])
        if dist < min_dist:
            min_dist = dist
            best_idx = idx
    
    return best_idx


def compute_md_rotamer_populations(
    trajectory,
    topology,
    chain,
    resnum,
    rotamer_library,
    chi_atom_names
):
    """
    Compute rotamer state populations from MD trajectory.
    
    Parameters
    ----------
    trajectory : str
        Path to MD trajectory file.
    topology : str
        Path to topology file.
    chain : str
        Chain identifier.
    resnum : int
        Residue number.
    rotamer_library : list of dict
        Backbone-dependent rotamer library for this residue.
    chi_atom_names : list of list
        Atom names defining each chi dihedral.
    
    Returns
    -------
    dict
        Mapping of rotamer_idx -> population fraction.
    
    Notes
    -----
    Requires MDAnalysis or MDTraj for trajectory analysis.
    """
    import mdtraj as md
    
    traj = md.load(trajectory, top=topology)
    
    # Get atom indices for chi dihedrals
    # (implementation depends on trajectory format)
    
    rotamer_counts = Counter()
    
    for frame in range(traj.n_frames):
        # Compute chi angles for this frame
        chi_angles = compute_chi_angles(traj, frame, chain, resnum, chi_atom_names)
        
        # Assign to closest library rotamer
        rot_idx = assign_rotamer_from_chi_angles(chi_angles, rotamer_library)
        rotamer_counts[rot_idx] += 1
    
    # Convert to fractions
    total = sum(rotamer_counts.values())
    populations = {idx: count / total for idx, count in rotamer_counts.items()}
    
    return populations


def compute_md_delta_g(
    target_rotamer_idx,
    native_rotamer_idx,
    md_populations,
    min_population=0.001
):
    """
    Compute conformational free energy change using MD-derived populations.
    
    Parameters
    ----------
    target_rotamer_idx : int
        Index of proposed rotamer state.
    native_rotamer_idx : int
        Index of native rotamer state.
    md_populations : dict
        Mapping of rotamer_idx -> population fraction from MD.
    min_population : float
        Minimum population to avoid log(0). Default 0.001 (0.1%).
        States with lower population are penalized as if at this level.
    
    Returns
    -------
    float
        ΔG in kcal/mol. Positive = less favorable than native.
    
    Notes
    -----
    The minimum population floor prevents infinite penalties for states
    never observed in MD. A value of 0.001 corresponds to ~4.1 kcal/mol
    penalty at 300K, which is substantial but not prohibitive.
    """
    p_native = max(md_populations.get(native_rotamer_idx, min_population), min_population)
    p_target = max(md_populations.get(target_rotamer_idx, min_population), min_population)
    
    delta_g = -KT_300K * math.log(p_target / p_native)
    
    return delta_g
```

### 6.5 Practical Considerations

#### MD Simulation Protocol

| Parameter | Recommended Value | Justification |
|-----------|-------------------|---------------|
| Simulation length | 50–100 ns | Sufficient to sample common rotamer transitions |
| Force field | ff19SB or CHARMM36m | Modern force fields with accurate sidechain dynamics |
| Temperature | 300K | Physiological |
| Sampling interval | 10 ps | Sufficient for rotamer state assignment |
| Equilibration | 10 ns | Remove initial structure bias |

#### Computational Cost

- Running 50 ns MD for each sample structure is computationally expensive
- **Practical compromise:** Run MD on a single representative structure per protein system; assume populations are transferable to other samples with similar backbone conformations

### 6.6 Literature Support

- **Fischer et al. (2014):** "The crystallographically refined occupancies of these conformations, which are observable in an apo receptor structure, define energy penalties for docking."
- **Merski et al. (2021):** Demonstrated that MD-derived conformational penalties improve docking selectivity in T4 lysozyme L99A model system.
- **Bowman et al. (2012):** Showed that MD simulations can capture rare conformational states relevant to ligand binding.

---

## 7. Proposed Enhancement 5: Shape Complementarity Scoring

### 7.1 Physical Basis

Binding affinity correlates strongly with **geometric complementarity** between ligand and binding site. Lawrence & Colman (1993) introduced the **shape complementarity (SC) score** to quantify this:

> "The shape complementarity statistic, Sc, is defined at each point on one surface as the degree to which the surface normals of the two molecules point in opposite directions along the line joining them."

An SC score ranges from 0 (no complementarity) to 1 (perfect complementarity). High-affinity protein-ligand complexes typically have SC > 0.7.

**A large, featureless cavity has poor shape complementarity with any specific ligand, leading to non-selective binding.**

### 7.2 Mathematical Formulation

The SC score at surface point `i` on molecule A is:

$$S_i = \vec{n}_i \cdot \vec{n}_j \cdot e^{-w \cdot d_{ij}^2}$$

Where:
- `n_i` = surface normal at point `i` on molecule A
- `n_j` = surface normal at nearest point `j` on molecule B
- `d_ij` = distance between points
- `w` = distance decay parameter

The overall SC score is the median of all `S_i` values.

### 7.3 Implementation

```python
import numpy as np
from scipy.spatial import cKDTree

def compute_molecular_surface(struct, probe_radius=1.4, grid_spacing=0.5):
    """
    Compute molecular surface and normals.
    
    Parameters
    ----------
    struct : Structure
        Protein or ligand structure.
    probe_radius : float
        Probe radius in Å. Default 1.4 (water).
    grid_spacing : float
        Surface point spacing in Å.
    
    Returns
    -------
    tuple
        (surface_points, surface_normals) as numpy arrays.
    
    Notes
    -----
    Can use Schrodinger's surface calculation or external tools like MSMS.
    """
    # Placeholder - use Schrodinger or MSMS implementation
    from schrodinger.structutils.analyze import generate_molecular_surface
    surface = generate_molecular_surface(struct, probe_radius)
    return surface.points, surface.normals


def compute_shape_complementarity(
    pocket_struct,
    ligand_struct,
    distance_cutoff=3.5,
    w=0.5
):
    """
    Compute Lawrence-Colman shape complementarity score.
    
    Parameters
    ----------
    pocket_struct : Structure
        Binding pocket structure (residues only).
    ligand_struct : Structure
        Ligand structure.
    distance_cutoff : float
        Maximum distance for surface point pairing (Å).
    w : float
        Distance decay parameter.
    
    Returns
    -------
    float
        Shape complementarity score [0, 1].
    
    References
    ----------
    Lawrence & Colman (1993) J Mol Biol 234:946-950
    """
    # Compute molecular surfaces
    pocket_points, pocket_normals = compute_molecular_surface(pocket_struct)
    ligand_points, ligand_normals = compute_molecular_surface(ligand_struct)
    
    # Build KD-tree for efficient nearest neighbor search
    ligand_tree = cKDTree(ligand_points)
    
    sc_values = []
    
    for i, (point, normal) in enumerate(zip(pocket_points, pocket_normals)):
        # Find nearest ligand surface point
        dist, j = ligand_tree.query(point)
        
        if dist > distance_cutoff:
            continue
        
        # Compute SC at this point
        # Normals should point in opposite directions for good complementarity
        normal_dot = -np.dot(normal, ligand_normals[j])  # Negative because opposing
        distance_weight = np.exp(-w * dist ** 2)
        
        sc_i = normal_dot * distance_weight
        sc_values.append(sc_i)
    
    if not sc_values:
        return 0.0
    
    # Return median (robust to outliers)
    return np.median(sc_values)


def shape_complementarity_penalty(sc_score, sc_target=0.7, k=2.0):
    """
    Penalty for poor shape complementarity.
    
    Parameters
    ----------
    sc_score : float
        Computed SC score [0, 1].
    sc_target : float
        Target SC score. Default 0.7 based on typical high-affinity complexes.
    k : float
        Penalty scaling factor.
    
    Returns
    -------
    float
        Penalty (negative if SC < target, zero otherwise).
    """
    if sc_score >= sc_target:
        return 0.0
    
    return -k * (sc_target - sc_score)
```

### 7.4 Literature Support

- **Lawrence & Colman (1993):** Original SC statistic paper. High-affinity antibody-antigen interfaces have SC = 0.64–0.68.
- **Nicholls et al. (1991):** "Protein folding and association: insights from the interfacial and thermodynamic properties of hydrocarbons." Established importance of shape complementarity in molecular recognition.

---

## 8. Proposed Enhancement 6: Desolvation Penalty

### 8.1 Physical Basis

When a ligand binds, both the ligand and the binding pocket must **desolvate**—water molecules are displaced from the interface. This has an energetic cost:

$$\Delta G_{binding} = \Delta G_{intrinsic} + \Delta G_{desolvation}^{ligand} + \Delta G_{desolvation}^{pocket}$$

Large, open pockets expose more surface area to solvent. When this surface is hydrophobic, desolvation upon ligand binding is **energetically favorable** (hydrophobic effect). However, **exposing charged or polar groups** without compensating interactions is unfavorable.

For our purposes, we need to consider:
1. Increased hydrophobic SASA = unfavorable (needs to be desolvated)
2. Exposed polar groups without ligand contacts = unfavorable

### 8.2 Mathematical Formulation

The desolvation penalty for changing pocket conformation:

$$P_{desolvation} = \gamma_{hydrophobic} \cdot \Delta SASA_{hydrophobic} + \gamma_{polar} \cdot \Delta SASA_{polar}^{uncompensated}$$

Where:
- `ΔSASA_hydrophobic` = change in hydrophobic solvent-accessible surface area
- `ΔSASA_polar_uncompensated` = change in polar SASA not compensated by new interactions
- `γ_hydrophobic` ≈ 0.025 kcal/mol/Å² (hydrophobic transfer free energy)
- `γ_polar` ≈ -0.017 kcal/mol/Å² (polar groups prefer solvation)

### 8.3 Implementation

```python
from schrodinger.structutils.analyze import calculate_sasa_by_atom

# Surface tension coefficients (kcal/mol/Å²)
# From Sitkoff et al. (1994) and Sharp et al. (1991)
GAMMA_HYDROPHOBIC = 0.025
GAMMA_POLAR = -0.017

# Hydrophobic atom types
HYDROPHOBIC_ELEMENTS = {'C'}
POLAR_ELEMENTS = {'N', 'O', 'S'}

def compute_sasa_by_type(struct, binding_residues):
    """
    Compute hydrophobic and polar SASA for binding site.
    
    Parameters
    ----------
    struct : Structure
        Protein structure.
    binding_residues : list of tuple
        List of (chain, resnum, resname) for binding site.
    
    Returns
    -------
    dict
        {'hydrophobic': float, 'polar': float} in Å².
    """
    # Get atoms in binding site
    binding_atoms = []
    for chain, resnum, resname in binding_residues:
        residue = get_residue(struct, chain, resnum)
        binding_atoms.extend(list(residue.atom))
    
    # Compute per-atom SASA
    atom_sasa = calculate_sasa_by_atom(struct)
    
    hydrophobic_sasa = 0.0
    polar_sasa = 0.0
    
    for atom in binding_atoms:
        sasa = atom_sasa.get(atom.index, 0)
        element = atom.element
        
        if element in HYDROPHOBIC_ELEMENTS:
            hydrophobic_sasa += sasa
        elif element in POLAR_ELEMENTS:
            polar_sasa += sasa
    
    return {'hydrophobic': hydrophobic_sasa, 'polar': polar_sasa}


def compute_desolvation_penalty(
    struct_before,
    struct_after,
    binding_residues,
    gamma_hydrophobic=GAMMA_HYDROPHOBIC,
    gamma_polar=GAMMA_POLAR
):
    """
    Compute desolvation penalty for pocket opening.
    
    Parameters
    ----------
    struct_before : Structure
        Structure before rotamer change.
    struct_after : Structure
        Structure after rotamer change.
    binding_residues : list of tuple
        Binding site residues.
    gamma_hydrophobic : float
        Surface tension for hydrophobic atoms (kcal/mol/Å²).
    gamma_polar : float
        Surface tension for polar atoms (kcal/mol/Å²).
    
    Returns
    -------
    float
        Desolvation penalty in kcal/mol.
        Positive = unfavorable change.
    
    Notes
    -----
    A positive penalty means the rotamer change increases exposed
    hydrophobic surface area (which will need to be desolvated upon
    ligand binding) without compensating favorable interactions.
    """
    sasa_before = compute_sasa_by_type(struct_before, binding_residues)
    sasa_after = compute_sasa_by_type(struct_after, binding_residues)
    
    delta_hydrophobic = sasa_after['hydrophobic'] - sasa_before['hydrophobic']
    delta_polar = sasa_after['polar'] - sasa_before['polar']
    
    # Hydrophobic SASA increase = unfavorable (positive penalty)
    # Polar SASA increase = favorable for solvation (negative penalty)
    penalty = gamma_hydrophobic * delta_hydrophobic + gamma_polar * delta_polar
    
    return penalty
```

### 8.4 Literature Support

- **Sitkoff et al. (1994):** Parameterization of atomic solvation parameters for calculating solvation free energies.
- **Sharp et al. (1991):** "Reconciling the magnitude of the microscopic and macroscopic hydrophobic effects." Established surface tension coefficients.
- **Mysinger & Shoichet (2010):** Demonstrated that ligand desolvation penalties improve docking enrichment.

---

## 9. Integrated Multi-Objective Function

### 9.1 Proposed Objective Function

Combining all enhancements, the new objective function for rotamer selection is:

$$\text{Score} = w_V \cdot P_{volume}(V, V_{target}) + w_{SC} \cdot SC - w_{\Delta G} \cdot \Delta G_{MD} - w_{WCA} \cdot E_{WCA} - w_{HS} \cdot P_{hotspot} - w_{DS} \cdot P_{desolvation}$$

Where:
- `P_volume` = Gaussian volume penalty (from Enhancement 1)
- `SC` = Shape complementarity score (from Enhancement 5)
- `ΔG_MD` = MD-derived conformational penalty (from Enhancement 4)
- `E_WCA` = WCA steric clash energy (existing)
- `P_hotspot` = Hot spot disruption penalty (from Enhancement 3)
- `P_desolvation` = Desolvation penalty (from Enhancement 6)

### 9.2 Recommended Weights

| Weight | Recommended Value | Justification |
|--------|-------------------|---------------|
| `w_V` | 1.0 | Base weight |
| `w_SC` | 0.5 | Secondary importance |
| `w_ΔG` | 1.0 | Energetically equivalent to volume |
| `w_WCA` | 1.0 | Maintains existing clash avoidance |
| `w_HS` | 2.0 | Higher weight to preserve druggability |
| `w_DS` | 0.5 | Moderate desolvation consideration |

### 9.3 Hard Constraints

In addition to the soft penalties in the objective function, apply **hard constraints** that reject rotamer states regardless of score:

1. **Volume upper bound:** `V_pocket ≤ 2.0 × V_ligand`
2. **WCA threshold:** `E_WCA ≤ 10 REU`
3. **Anchor residue exclusion:** Anchors not sampled

### 9.4 Implementation

```python
def compute_rotamer_score(
    struct_orig,
    struct_new,
    chain,
    resnum,
    ligand_struct,
    binding_residues,
    hot_spot_preserver,
    md_populations,
    rotamer_library,
    native_rotamer_idx,
    target_rotamer_idx,
    V_target,
    weights=None
):
    """
    Compute integrated score for a proposed rotamer state.
    
    Parameters
    ----------
    struct_orig : Structure
        Original structure.
    struct_new : Structure
        Structure with proposed rotamer applied.
    chain : str
        Chain identifier.
    resnum : int
        Residue number.
    ligand_struct : Structure
        Reference ligand structure.
    binding_residues : list of tuple
        Binding site residues.
    hot_spot_preserver : HotSpotPreserver
        Hot spot preservation object.
    md_populations : dict
        MD-derived rotamer populations.
    rotamer_library : list
        Backbone-dependent rotamer library.
    native_rotamer_idx : int
        Index of native rotamer.
    target_rotamer_idx : int
        Index of proposed rotamer.
    V_target : float
        Target pocket volume (Å³).
    weights : dict, optional
        Weight parameters.
    
    Returns
    -------
    dict
        Score components and total score.
    """
    if weights is None:
        weights = {
            'w_V': 1.0,
            'w_SC': 0.5,
            'w_dG': 1.0,
            'w_WCA': 1.0,
            'w_HS': 2.0,
            'w_DS': 0.5
        }
    
    # Compute pocket volume
    V_pocket = measure_binding_site_volume(struct_new, binding_residues)
    
    # Hard constraint: volume upper bound
    V_max = 2.0 * compute_ligand_volume(ligand_struct)
    if V_pocket > V_max:
        return {'rejected': True, 'reason': 'volume_exceeded'}
    
    # 1. Volume penalty
    score_volume = volume_penalty(V_pocket, V_target, sigma=50.0)
    
    # 2. Shape complementarity
    sc_score = compute_shape_complementarity(struct_new, ligand_struct)
    
    # 3. MD-derived ΔG
    delta_g = compute_md_delta_g(target_rotamer_idx, native_rotamer_idx, md_populations)
    
    # 4. WCA energy
    wca_energy = get_wca_energy(struct_new, chain, resnum)
    
    # Hard constraint: WCA threshold
    if wca_energy is not None and wca_energy > 10.0:
        return {'rejected': True, 'reason': 'wca_exceeded'}
    
    # 5. Hot spot penalty
    hs_penalty = hot_spot_preserver.compute_hotspot_penalty(
        struct_orig, struct_new, chain, resnum
    )
    
    # 6. Desolvation penalty
    ds_penalty = compute_desolvation_penalty(struct_orig, struct_new, binding_residues)
    
    # Combine scores
    total_score = (
        weights['w_V'] * score_volume +
        weights['w_SC'] * sc_score -
        weights['w_dG'] * delta_g -
        weights['w_WCA'] * (wca_energy or 0) -
        weights['w_HS'] * abs(hs_penalty) -
        weights['w_DS'] * ds_penalty
    )
    
    return {
        'rejected': False,
        'total_score': total_score,
        'volume': V_pocket,
        'volume_penalty': score_volume,
        'shape_complementarity': sc_score,
        'delta_g': delta_g,
        'wca_energy': wca_energy,
        'hotspot_penalty': hs_penalty,
        'desolvation_penalty': ds_penalty
    }
```

---

## 10. Implementation Priority

Based on expected impact and implementation complexity:

| Priority | Enhancement | Complexity | Expected Impact | Rationale |
|----------|-------------|------------|-----------------|-----------|
| **1** | Reference ligand volume constraint | Low | High | Directly prevents over-opening; simple to implement |
| **2** | Anchor residue exclusion | Low | High | Preserves key specificity determinants |
| **3** | Hot spot preservation (FTMap) | Medium | High | Maintains druggability features |
| **4** | MD-derived ΔG penalties | Medium | Medium | More accurate conformational energetics |
| **5** | Shape complementarity | Medium | Medium | Maintains ligand selectivity |
| **6** | Desolvation penalty | Low | Medium | Physical realism for binding thermodynamics |

### 10.1 Suggested Implementation Roadmap

#### Phase 1: Quick Wins (1-2 days)

1. Implement volume constraint with Gaussian penalty
2. Add ligand volume calculation
3. Implement anchor residue identification (contact-based)
4. Integrate anchor exclusion into sampling loop

#### Phase 2: Hot Spot Integration (3-5 days)

1. Set up FTMap submission workflow
2. Parse FTMap results
3. Implement HotSpotPreserver class
4. Integrate hot spot penalty into scoring

#### Phase 3: MD Integration (1-2 weeks)

1. Establish MD simulation protocol
2. Implement trajectory analysis for rotamer populations
3. Replace/augment library-based ΔG with MD-derived values

#### Phase 4: Refinements (1 week)

1. Implement shape complementarity scoring
2. Add desolvation penalty
3. Tune weights based on retrospective validation

---

## 11. References

1. **Bogan AA, Thorn KS.** (1998) Anatomy of hot spots in protein interfaces. *J Mol Biol* 280:1-9.

2. **Bowman GR, Huang X, Pande VS.** (2009) Using generalized ensemble simulations and Markov state models to identify conformational states. *Methods* 49:197-201.

3. **Brenke R, Kozakov D, Chuang GY, et al.** (2009) Fragment-based identification of druggable 'hot spots' of proteins using Fourier domain correlation techniques. *Bioinformatics* 25:621-627.

4. **Clackson T, Wells JA.** (1995) A hot spot of binding energy in a hormone-receptor interface. *Science* 267:383-386.

5. **Fischer M, Coleman RG, Fraser JS, Shoichet BK.** (2014) Incorporation of protein flexibility and conformational energy penalties in docking screens to improve ligand discovery. *Nat Chem* 6:575-583.

6. **Lawrence MC, Colman PM.** (1993) Shape complementarity at protein/protein interfaces. *J Mol Biol* 234:946-950.

7. **Liang J, Edelsbrunner H, Woodward C.** (1998) Anatomy of protein pockets and cavities: measurement of binding site geometry and implications for ligand design. *Protein Sci* 7:1884-1897.

8. **Merski M, Fischer M, Balius TE, et al.** (2021) Energy penalties enhance flexible receptor docking in a model cavity. *Proc Natl Acad Sci USA* 118:e2106195118.

9. **Mysinger MM, Shoichet BK.** (2010) Rapid context-dependent ligand desolvation in molecular docking. *J Chem Inf Model* 50:1561-1573.

10. **Nicholls A, Sharp KA, Honig B.** (1991) Protein folding and association: insights from the interfacial and thermodynamic properties of hydrocarbons. *Proteins* 11:281-296.

11. **Rajamani D, Thiel S, Vajda S, Camacho CJ.** (2004) Anchor residues in protein-protein interactions. *Proc Natl Acad Sci USA* 101:11287-11292.

12. **Sharp KA, Nicholls A, Fine RF, Honig B.** (1991) Reconciling the magnitude of the microscopic and macroscopic hydrophobic effects. *Science* 252:106-109.

13. **Sitkoff D, Sharp KA, Honig B.** (1994) Accurate calculation of hydration free energies using macroscopic solvent models. *J Phys Chem* 98:1978-1988.

14. **Vajda S, Whitty A, Kozakov D.** (2018) Druggability assessment. In: *Encyclopedia of Biophysics*. Springer, Berlin.

15. **Zerbe BS, Hall DR, Vajda S, Whitty A, Bhardwaj N.** (2012) Relationship between hot spot residues and ligand binding hot spots in protein-protein interfaces. *J Chem Inf Model* 52:2236-2244.

---

## Appendix A: Quick Reference for Key Equations

### A.1 Volume Penalty
$$P_{volume}(V) = -\frac{(V - V_{target})^2}{2\sigma^2}$$

### A.2 Boltzmann Conformational Penalty
$$\Delta G = -k_BT \ln\left(\frac{p_{suggested}}{p_{native}}\right)$$

### A.3 Shape Complementarity
$$S_i = \vec{n}_i \cdot \vec{n}_j \cdot e^{-w \cdot d_{ij}^2}$$

### A.4 Desolvation Penalty
$$P_{desolvation} = \gamma_{hydrophobic} \cdot \Delta SASA_{hydrophobic} + \gamma_{polar} \cdot \Delta SASA_{polar}$$

### A.5 Hot Spot Disruption Penalty
$$P_{hotspot} = -k \cdot \max(0, d_{new} - d_{orig})$$

### A.6 Integrated Score
$$\text{Score} = w_V P_V + w_{SC} \cdot SC - w_{\Delta G} \Delta G - w_{WCA} E_{WCA} - w_{HS} P_{HS} - w_{DS} P_{DS}$$

---

## Appendix B: Parameter Summary Table

| Parameter | Symbol | Recommended Value | Units | Source |
|-----------|--------|-------------------|-------|--------|
| Volume expansion factor | α | 1.5 | - | Liang et al. (1998) |
| Volume penalty σ | σ | 50 | Å³ | Empirical |
| Volume hard constraint | V_max | 2.0 × V_ligand | Å³ | Conservative |
| WCA hard threshold | - | 10 | REU | Empirical |
| Contact threshold | - | 3 | contacts | Structural analysis |
| Conservation threshold | - | 0.8 | fraction | Bogan & Thorn (1998) |
| ΔΔG anchor threshold | - | 2.0 | kcal/mol | Clackson & Wells (1995) |
| FTMap probe threshold | - | 16 | probes | Vajda et al. (2018) |
| Hot spot penalty k | k | 0.5 | kcal/mol/Å | van der Waals scale |
| kT at 300K | kT | 0.6 | kcal/mol | Physical constant |
| γ hydrophobic | γ_h | 0.025 | kcal/mol/Å² | Sharp et al. (1991) |
| γ polar | γ_p | -0.017 | kcal/mol/Å² | Sharp et al. (1991) |
| SC target | - | 0.7 | - | Lawrence & Colman (1993) |

---

*Document generated March 24, 2026*
