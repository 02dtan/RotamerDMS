"""
Rotamer sampling phase module for RotamerDMS.

Samples rotamer states from Schrodinger's backbone-dependent rotamer library
and selects states that maximize binding pocket volume.
"""

import copy
import json
import math
import os
from datetime import datetime

from schrodinger import structure
from schrodinger.protein import rotamers
from schrodinger.structutils.measure import measure_dihedral_angle

from .cavity import measure_binding_site_volume
from .wca_potential import minimize_and_get_fa_rep, _check_pyrosetta_available
from .scoring import (
    rank_rotamers_by_joint_score,
    analyze_pareto_tradeoff,
    compute_normalization_params
)

# Residues with symmetric side chains and which chi angle has 180° symmetry
# For these residues, a 180° flip of the specified chi angle produces an equivalent conformation
SYMMETRIC_CHI_ANGLES = {
    'PHE': 2,  # Chi2 has 180° symmetry (benzene ring flip)
    'TYR': 2,  # Chi2 has 180° symmetry (phenol ring flip)
    'ASP': 2,  # Chi2 has 180° symmetry (carboxylate)
    'GLU': 3,  # Chi3 has 180° symmetry (carboxylate)
}


def get_residue_object(struct, chain, resnum):
    """
    Get a residue object from structure by chain and residue number.
    
    Args:
        struct: Schrodinger Structure object
        chain: Chain identifier
        resnum: Residue number
        
    Returns:
        Residue object or None if not found
    """
    for res in struct.residue:
        if res.chain == chain and res.resnum == resnum:
            return res
    return None


def get_current_chi_angles(struct, residue):
    """
    Measure the current chi angles of a residue in the structure.
    
    Uses residue.getDihedralAtoms() to get atom objects for each chi angle,
    then measures the dihedral angle using measure_dihedral_angle.
    
    Args:
        struct: Schrodinger Structure object (unused, kept for API consistency)
        residue: Schrodinger Residue object
        
    Returns:
        List of chi angles (in degrees)
    """
    chi_angles = []
    # Schrodinger supports Chi1 through Chi5 depending on the residue type
    chi_names = ["Chi1", "Chi2", "Chi3", "Chi4", "Chi5"]
    
    for chi_name in chi_names:
        try:
            # getDihedralAtoms returns a list of 4 atom objects
            atoms = residue.getDihedralAtoms(chi_name)
            if atoms and len(atoms) == 4:
                angle = measure_dihedral_angle(*atoms)
                chi_angles.append(angle)
        except (ValueError, AttributeError, IndexError):
            # ValueError raised if dihedral not found in database for this residue
            break
    return chi_angles


def compute_chi_angular_distance(current_chi, library_chi, resname):
    """
    Compute angular distance between two sets of chi angles.
    
    Handles periodicity (angles wrap at ±180°) and symmetry for residues
    with symmetric side chains (PHE, TYR have 180° symmetry on chi2).
    
    Args:
        current_chi: List of current chi angles
        library_chi: List of library rotamer chi angles
        resname: 3-letter residue code (for symmetry handling)
        
    Returns:
        Angular distance (Euclidean in chi-space), or None if incompatible
    """
    if not library_chi or len(library_chi) != len(current_chi):
        return None
    
    # Check if this residue has a symmetric chi angle
    symmetric_chi_idx = SYMMETRIC_CHI_ANGLES.get(resname)  # 1-indexed
    
    dist_sq = 0.0
    for chi_idx, (c_angle, l_angle) in enumerate(zip(current_chi, library_chi), start=1):
        diff = c_angle - l_angle
        
        # Handle angle periodicity: wrap to [-180, 180]
        while diff > 180:
            diff -= 360
        while diff < -180:
            diff += 360
        
        # For symmetric chi angles, also consider 180° flip
        if chi_idx == symmetric_chi_idx:
            # The 180° flipped difference
            diff_flipped = diff + 180 if diff < 0 else diff - 180
            # Use whichever is smaller in absolute value
            diff = diff if abs(diff) <= abs(diff_flipped) else diff_flipped
        
        dist_sq += diff * diff
    
    return math.sqrt(dist_sq)


def find_closest_library_rotamer(struct, chain, resnum):
    """
    Find the library rotamer that best matches the current conformation.
    
    Compares chi angles of the residue in the structure to chi angles
    of each rotamer in the backbone-dependent library. The rotamer with
    the smallest angular deviation is considered the matching state.
    
    Args:
        struct: Schrodinger Structure object
        chain: Chain identifier
        resnum: Residue number
        
    Returns:
        Dictionary containing:
            - 'index': 1-indexed rotamer index (or None if not found)
            - 'probability': Probability of the closest rotamer
            - 'angular_distance': Angular deviation in degrees
        Returns None if residue not found or has no rotamers
    """
    res = get_residue_object(struct, chain, resnum)
    if res is None:
        return None
    
    resname = res.pdbres.strip()
    if resname in ['ALA', 'GLY']:
        return None
    
    try:
        # Get current chi angles from structure
        current_chi = get_current_chi_angles(struct, res)
        if not current_chi:
            return None
        
        # Load rotamer library
        rotamer_lib = rotamers.Rotamers(struct, res.atom[1])
        all_rotamers = rotamer_lib.rotamers
        
        if not all_rotamers:
            return None
        
        best_idx = None
        min_dist = float('inf')
        best_prob = 0.0
        
        # Iterate through library rotamers (1-indexed)
        for i, rot in enumerate(all_rotamers, start=1):
            try:
                # Use chiAngles method (not chi attribute)
                library_chi = rot.chiAngles()
                
                # Compute angular distance with periodicity and symmetry handling
                dist = compute_chi_angular_distance(current_chi, library_chi, resname)
                if dist is None:
                    continue
                
                if dist < min_dist:
                    min_dist = dist
                    best_idx = i
                    # Get percentage from rotamer (percentage is a method)
                    try:
                        best_prob = rot.percentage()
                    except (AttributeError, TypeError):
                        best_prob = 0.0
            except Exception:
                continue
        
        if best_idx is None:
            return None
        
        return {
            'index': best_idx,
            'probability': best_prob,
            'angular_distance': min_dist
        }
        
    except Exception as e:
        return None


def get_rotamer_library(struct, chain, resnum):
    """
    Load the backbone-dependent rotamer library for a residue.
    
    Args:
        struct: Schrodinger Structure object
        chain: Chain identifier
        resnum: Residue number
        
    Returns:
        Tuple of (Rotamers object, list of rotamer states) or (None, None) if failed
    """
    res = get_residue_object(struct, chain, resnum)
    if res is None:
        print(f"Warning: Residue {chain}:{resnum} not found.")
        return None, None
    
    resname = res.pdbres.strip()
    
    # Ala and Gly have no rotamers
    # This should never come up, as Ala/Gly are pre-filtered out in the mutation phase
    # But this is a good cheap sanity check to keep
    if resname in ['ALA', 'GLY']:
        print(f"Skipping {chain}:{resname}{resnum} - no rotamers for Ala/Gly")
        return None, None
    
    try:
        # Load rotamer library using an atom contained in the residue
        rotamer_lib = rotamers.Rotamers(struct, res.atom[1])
        all_rotamers = rotamer_lib.rotamers
        
        if not all_rotamers:
            print(f"Warning: No rotamers found for {chain}:{resname}{resnum}")
            return None, None
        
        return rotamer_lib, all_rotamers
    except Exception as e:
        print(f"Error loading rotamer library for {chain}:{resname}{resnum}: {e}")
        return None, None


def sample_rotamers_for_residue(struct, chain, resnum, resname, binding_site_residues, 
                                 initial_volume, min_contacts=4, compute_wca=True,
                                 w_vol=1.0, w_dg=1.0, w_wca=1.0, wca_threshold=None,
                                 verbose=True):
    """
    Sample all rotamer states for a single residue and rank by joint optimization score.
    
    Args:
        struct: Schrodinger Structure object (working copy)
        chain: Chain identifier
        resnum: Residue number
        resname: Residue name (3-letter code)
        binding_site_residues: List of binding site residue tuples
        initial_volume: Current binding site volume before rotamer changes
        compute_wca: Whether to compute WCA repulsive potential (default True)
        w_vol: Weight for volume in joint optimization (default 1.0)
        w_dg: Weight for deltaG in joint optimization (default 1.0)
        w_wca: Weight for WCA in joint optimization (default 1.0)
        wca_threshold: Hard WCA threshold to reject clashing rotamers (default None)
        verbose: Whether to print progress
        
    Returns:
        Dictionary containing:
            - 'rotamer_states': List of rotamer indices sorted by joint score
            - 'volume_changes': List of volume changes (same order as rotamer_states)
            - 'percentages': List of rotamer percentages (same order)
            - 'delta_g_values': List of deltaG values (same order) in arbitrary units (AU)
            - 'wca_energies': List of WCA repulsive energies (same order) in Rosetta energy units
            - 'joint_scores': List of joint optimization scores (same order)
            - 'best_rotamer_idx': Index of best rotamer by joint score
            - 'best_volume_change': Volume change from best rotamer
            - 'native_rotamer_idx': Index of closest library rotamer to current conformation
            - 'native_probability': Probability of native/current rotamer state
    """
    if verbose:
        print(f"\n  Sampling rotamers for {chain}:{resname}{resnum}...")
    
    # Get rotamer library
    rotamer_lib, all_rotamers = get_rotamer_library(struct, chain, resnum)
    
    if rotamer_lib is None:
        return None
    
    if verbose:
        print(f"    Found {len(all_rotamers)} rotamer states in library")
    
    # Find the current/native rotamer state by comparing chi angles
    native_rotamer_info = find_closest_library_rotamer(struct, chain, resnum)
    
    if native_rotamer_info is None:
        if verbose:
            print(f"    Warning: Could not identify native rotamer state")
        p_curr = None
        native_rotamer_idx = None
    else:
        native_rotamer_idx = native_rotamer_info['index']
        p_curr = native_rotamer_info['probability']
        if verbose:
            print(f"    Native rotamer: idx={native_rotamer_idx}, "
                  f"p_curr={p_curr:.4f}, angular_dist={native_rotamer_info['angular_distance']:.1f}°")
    
    # Store results for each rotamer
    rotamer_results = []
    
    # Rotamer indices are 1-indexed in Schrodinger's rotamer library
    for rot_idx, rotamer in enumerate(all_rotamers, start=1):
        # Create a copy of the structure
        test_struct = struct.copy()
        
        # Load rotamer library for the copy
        test_rotamer_lib, test_rotamers = get_rotamer_library(test_struct, chain, resnum)
        
        if test_rotamer_lib is None or rot_idx > len(test_rotamers):
            continue
        
        # Apply this rotamer state using 1-indexed key
        try:
            test_rotamers[rot_idx].apply()
        except Exception as e:
            if verbose:
                print(f"    Warning: Failed to apply rotamer {rot_idx}: {e}")
            continue
        
        # Get rotamer percentage (frequency proxy) - computed BEFORE minimization
        # This represents the statistical probability of this rotamer conformation
        try:
            pct = rotamer.percentage
            # Handle if percentage is a method vs property
            percentage = pct() if callable(pct) else float(pct)
        except (AttributeError, TypeError):
            percentage = 0.0
        
        # Compute deltaG = -ln(p_sugg / p_curr) in arbitrary units (AU)
        # This represents the relative free energy change from native to this rotamer
        # Note: deltaG = 0 when p_sugg = p_curr (i.e., suggested rotamer is the native state)
        if p_curr is not None and p_curr > 0 and percentage > 0:
            delta_g = -math.log(percentage / p_curr)
        else:
            # Cannot compute deltaG if probabilities are missing/zero
            delta_g = None
        
        # Minimize structure and get fa_rep energy
        # The minimization relaxes surrounding residues while keeping the target
        # rotamer frozen, producing more physically realistic structures
        fa_rep_energy = None
        fa_dun_energy = None
        struct_for_volume = test_struct  # Default: use non-minimized structure
        
        if compute_wca:
            min_result = minimize_and_get_fa_rep(test_struct, chain, resnum, verbose=False)
            if min_result.get('success'):
                fa_rep_energy = min_result['fa_rep_energy']
                fa_dun_energy = min_result.get('fa_dun_energy')
                struct_for_volume = min_result['minimized_struct']
            # If minimization fails, we still proceed with non-minimized structure
        
        # Measure volume on the (possibly minimized) structure
        result = measure_binding_site_volume(
            struct_for_volume, 
            binding_site_residues,
            min_contacts=min_contacts,
            verbose=False
        )
        
        if not result['success']:
            continue
        
        new_volume = result['total_volume']
        delta_volume = new_volume - initial_volume
        
        rotamer_results.append({
            'index': rot_idx,
            'volume_change': delta_volume,
            'percentage': percentage,
            'delta_g': delta_g,
            'fa_dun': fa_dun_energy,
            'wca_energy': fa_rep_energy,
            'new_volume': new_volume
        })
    
    if not rotamer_results:
        if verbose:
            print(f"    No valid rotamer states found")
        return None
    
    # Rank by joint score (weighted combination of volume, deltaG, WCA)
    # Score = w_vol * z(ΔV) - w_dg * z(ΔG) - w_wca * z(WCA)
    joint_scored_results = rank_rotamers_by_joint_score(
        rotamer_results,
        w_vol=w_vol, w_dg=w_dg, w_wca=w_wca,
        wca_threshold=wca_threshold
    )
    
    # Also sort by volume only (for reference)
    volume_sorted_results = sorted(rotamer_results, key=lambda x: x['volume_change'], reverse=True)
    
    # Pareto front analysis (volume vs WCA tradeoff)
    pareto_analysis = analyze_pareto_tradeoff(rotamer_results, verbose=verbose)
    
    # Extract lists from joint-score-sorted results
    rotamer_states = [r['index'] for r in joint_scored_results]
    volume_changes = [r['volume_change'] for r in joint_scored_results]
    percentages = [r['percentage'] for r in joint_scored_results]
    delta_g_values = [r['delta_g'] for r in joint_scored_results]
    fa_dun_values = [r.get('fa_dun') for r in joint_scored_results]
    wca_energies = [r['wca_energy'] for r in joint_scored_results]
    joint_scores = [r['joint_score'] for r in joint_scored_results]
    
    # Best by joint score
    best = joint_scored_results[0] if joint_scored_results else volume_sorted_results[0]
    # Best by volume only (for comparison)
    best_by_volume = volume_sorted_results[0]
    
    if verbose:
        # Format metrics for display
        best_dg_str = f"{best['delta_g']:.2f} AU" if best.get('delta_g') is not None else "N/A"
        best_fadun_str = f"{best['fa_dun']:.2f} REU" if best.get('fa_dun') is not None else "N/A"
        best_wca_str = f"{best['wca_energy']:.2f} REU" if best.get('wca_energy') is not None else "N/A"
        best_score_str = f"{best['joint_score']:.3f}" if best.get('joint_score') is not None else "N/A"
        
        print(f"    Best rotamer (by joint score): index {best['index']}, "
              f"score={best_score_str}, "
              f"ΔV={best['volume_change']:+.2f} Å³, "
              f"ΔG={best_dg_str}, fa_dun={best_fadun_str}, "
              f"WCA={best_wca_str}")
        
        if best['index'] != best_by_volume['index']:
            vol_dg_str = f"{best_by_volume['delta_g']:.2f} AU" if best_by_volume.get('delta_g') is not None else "N/A"
            vol_fadun_str = f"{best_by_volume['fa_dun']:.2f} REU" if best_by_volume.get('fa_dun') is not None else "N/A"
            vol_wca_str = f"{best_by_volume['wca_energy']:.2f} REU" if best_by_volume.get('wca_energy') is not None else "N/A"
            print(f"    Best by volume only: index {best_by_volume['index']}, "
                  f"ΔV={best_by_volume['volume_change']:+.2f} Å³, "
                  f"ΔG={vol_dg_str}, fa_dun={vol_fadun_str}, WCA={vol_wca_str}")
        
        print(f"    Top 3 rotamers (by joint score) [ΔG=manual, fa_dun=Dunbrack]:")
        for i, r in enumerate(joint_scored_results[:3]):
            dg_str = f"{r['delta_g']:.2f}" if r.get('delta_g') is not None else "N/A"
            fadun_str = f"{r['fa_dun']:.2f}" if r.get('fa_dun') is not None else "N/A"
            wca_str = f"{r['wca_energy']:.2f}" if r.get('wca_energy') is not None else "N/A"
            score_str = f"{r['joint_score']:.3f}" if r.get('joint_score') is not None else "N/A"
            print(f"      {i+1}. idx={r['index']}: score={score_str}, "
                  f"ΔV={r['volume_change']:+.2f} Å³, ΔG={dg_str} AU, fa_dun={fadun_str} REU, WCA={wca_str} REU")
    
    return {
        'rotamer_states': rotamer_states,
        'volume_changes': volume_changes,
        'percentages': percentages,
        'delta_g_values': delta_g_values,
        'fa_dun_values': fa_dun_values,
        'wca_energies': wca_energies,
        'joint_scores': joint_scores,
        'best_rotamer_idx': best['index'],
        'best_volume_change': best['volume_change'],
        'best_percentage': best['percentage'],
        'best_delta_g': best.get('delta_g'),
        'best_fa_dun': best.get('fa_dun'),
        'best_wca_energy': best.get('wca_energy'),
        'best_joint_score': best.get('joint_score'),
        'best_by_volume_idx': best_by_volume['index'],
        'pareto_front': pareto_analysis['pareto_2d_vol_wca'],
        'pareto_analysis': pareto_analysis,
        'native_rotamer_idx': native_rotamer_idx,
        'native_probability': p_curr,
        'all_results': joint_scored_results
    }


def run_rotamer_sampling_phase(sample_struct, sorted_residues, binding_site_residues,
                                initial_volume, min_contacts=4, compute_wca=True,
                                w_vol=1.0, w_dg=1.0, w_wca=1.0, wca_threshold=None,
                                checkpoint_dir=None, verbose=True):
    """
    Run the rotamer sampling phase on sorted residues from mutation phase.
    
    Per design.md:
    - Process residues in order (highest volume impact first)
    - For each residue, sample all rotamer states
    - Record volume changes and percentages for each state
    - Rank by joint score combining volume, deltaG, and WCA
    - Compute WCA repulsive potential for each rotamer state
    
    Args:
        sample_struct: Schrodinger Structure object (original sample)
        sorted_residues: List of residue tuples sorted by volume impact from mutation phase
        binding_site_residues: Full list of binding site residues
        initial_volume: Initial binding site volume
        min_contacts: Minimum binding site residues a cavity must contact (default 4)
        compute_wca: Whether to compute WCA repulsive potential (default True)
        w_vol: Weight for volume in joint optimization (default 1.0)
        w_dg: Weight for deltaG in joint optimization (default 1.0)
        w_wca: Weight for WCA in joint optimization (default 1.0)
        wca_threshold: Hard WCA threshold to reject clashing rotamers (default None)
        checkpoint_dir: Optional directory to save checkpoints
        verbose: Whether to print detailed progress
        
    Returns:
        Dictionary containing:
            - 'residue_rotamers': Dict mapping residue to rotamer sampling results
            - 'best_rotamers': Dict mapping residue to best rotamer index
            - 'total_volume_potential': Estimated max volume if all best rotamers applied
    """
    if verbose:
        print("\n" + "="*60)
        print("ROTAMER SAMPLING PHASE: Finding optimal rotamer states")
        print("="*60)
        print(f"\nProcessing {len(sorted_residues)} residues in priority order...")
        if compute_wca:
            print("WCA repulsive potential computation: ENABLED")
            # Check if PyRosetta is available via subprocess
            if _check_pyrosetta_available():
                print("PyRosetta subprocess available")
            else:
                print("WARNING: PyRosetta not available, WCA computation disabled")
                compute_wca = False
    
    residue_rotamers = {}
    best_rotamers = {}
    
    # Working structure that accumulates best rotamer changes
    working_struct = sample_struct.copy()
    current_volume = initial_volume
    
    for rank, (chain, resnum, resname) in enumerate(sorted_residues, 1):
        if verbose:
            print(f"\n[{rank}/{len(sorted_residues)}] Processing {chain}:{resname}{resnum}")
        
        # Sample rotamers for this residue using current working structure
        result = sample_rotamers_for_residue(
            working_struct,
            chain, resnum, resname,
            binding_site_residues,
            current_volume,
            min_contacts=min_contacts,
            compute_wca=compute_wca,
            w_vol=w_vol, w_dg=w_dg, w_wca=w_wca,
            wca_threshold=wca_threshold,
            verbose=verbose
        )
        
        if result is None:
            if verbose:
                print(f"  Skipping - no valid rotamers")
            continue
        
        residue_rotamers[(chain, resnum, resname)] = result
        best_rotamers[(chain, resnum, resname)] = result['best_rotamer_idx']
        
        # Apply the best rotamer to working structure for subsequent residues
        if result['best_volume_change'] > 0:
            rotamer_lib, all_rotamers = get_rotamer_library(working_struct, chain, resnum)
            # best_rotamer_idx is 1-indexed, use directly with Schrodinger's 1-indexed access
            best_idx = result['best_rotamer_idx']
            if rotamer_lib and 1 <= best_idx <= len(all_rotamers):
                try:
                    all_rotamers[best_idx].apply()
                    current_volume += result['best_volume_change']
                    if verbose:
                        print(f"  Applied best rotamer, new volume: {current_volume:.2f} Å^3")
                except Exception as e:
                    if verbose:
                        print(f"  Warning: Failed to apply best rotamer: {e}")
    
    # Summary
    if verbose:
        print("\n" + "="*60)
        print("ROTAMER SAMPLING PHASE COMPLETE")
        print("="*60)
        print(f"\nResidues with valid rotamer data: {len(residue_rotamers)}")
        print(f"Estimated optimized volume: {current_volume:.2f} Å^3")
        print(f"Volume increase from initial: {current_volume - initial_volume:+.2f} Å^3")
    
    result = {
        'residue_rotamers': {str(k): v for k, v in residue_rotamers.items()},
        'best_rotamers': {str(k): v for k, v in best_rotamers.items()},
        'estimated_final_volume': current_volume,
        'volume_increase': current_volume - initial_volume
    }
    
    # Save checkpoint
    if checkpoint_dir:
        save_rotamer_checkpoint(result, sorted_residues, checkpoint_dir)
    
    return result, working_struct


def save_rotamer_checkpoint(result, sorted_residues, checkpoint_dir):
    """
    Save rotamer sampling results to checkpoint directory.
    """
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Prepare checkpoint data (simplify for JSON serialization)
    checkpoint_data = {
        'timestamp': timestamp,
        'phase': 'rotamer_sampling',
        'estimated_final_volume': result['estimated_final_volume'],
        'volume_increase': result['volume_increase'],
        'best_rotamers': [],
        'residue_details': []
    }
    
    for chain, resnum, resname in sorted_residues:
        key = str((chain, resnum, resname))
        if key in result['best_rotamers']:
            checkpoint_data['best_rotamers'].append({
                'chain': chain,
                'resnum': resnum,
                'resname': resname,
                'best_rotamer_idx': result['best_rotamers'][key]
            })
        
        if key in result['residue_rotamers']:
            res_data = result['residue_rotamers'][key]
            # Extract Pareto front indices for checkpoint
            pareto_indices = [r['index'] for r in res_data.get('pareto_front', [])]
            
            # Build all_rotamers list with full data for each rotamer
            all_rotamers = []
            for i, rot_idx in enumerate(res_data['rotamer_states']):
                all_rotamers.append({
                    'index': rot_idx,
                    'volume_change': res_data['volume_changes'][i],
                    'percentage': res_data['percentages'][i],
                    'delta_g': res_data['delta_g_values'][i],
                    'wca_energy': res_data['wca_energies'][i] if res_data.get('wca_energies') else None,
                    'joint_score': res_data['joint_scores'][i] if res_data.get('joint_scores') else None
                })
            
            checkpoint_data['residue_details'].append({
                'chain': chain,
                'resnum': resnum,
                'resname': resname,
                'best_rotamer_idx': res_data['best_rotamer_idx'],
                'best_volume_change': res_data['best_volume_change'],
                'best_percentage': res_data['best_percentage'],
                'best_delta_g': res_data.get('best_delta_g'),
                'best_wca_energy': res_data.get('best_wca_energy'),
                'best_joint_score': res_data.get('best_joint_score'),
                'best_by_volume_idx': res_data.get('best_by_volume_idx'),
                'pareto_front_indices': pareto_indices,
                'num_pareto_optimal': len(pareto_indices),
                'native_rotamer_idx': res_data.get('native_rotamer_idx'),
                'native_probability': res_data.get('native_probability'),
                'num_rotamers_tested': len(res_data['rotamer_states']),
                'all_rotamers': all_rotamers
            })
    
    checkpoint_path = os.path.join(checkpoint_dir, f'rotamer_sampling_{timestamp}.json')
    
    with open(checkpoint_path, 'w') as f:
        json.dump(checkpoint_data, f, indent=2)
    
    print(f"\nCheckpoint saved: {checkpoint_path}")
    
    return checkpoint_path


def apply_best_rotamers(sample_struct, best_rotamers, verbose=True):
    """
    Apply all best rotamer states to a structure.
    
    Args:
        sample_struct: Schrodinger Structure object to modify
        best_rotamers: Dict mapping (chain, resnum, resname) string to best rotamer index
        verbose: Whether to print progress
        
    Returns:
        Modified structure with best rotamers applied
    """
    modified_struct = sample_struct.copy()
    
    if verbose:
        print("\nApplying best rotamer states to structure...")
    
    for key_str, rotamer_idx in best_rotamers.items():
        # Parse the key string back to tuple
        # Key format: "('A', 102, 'MET')"
        try:
            key = eval(key_str)
            chain, resnum, resname = key
        except:
            if verbose:
                print(f"  Warning: Could not parse key {key_str}")
            continue
        
        rotamer_lib, all_rotamers = get_rotamer_library(modified_struct, chain, resnum)
        
        if rotamer_lib is None:
            continue
        
        # rotamer_idx is 1-indexed, use directly with Schrodinger's 1-indexed access
        if rotamer_idx < 1 or rotamer_idx > len(all_rotamers):
            if verbose:
                print(f"  Warning: Rotamer index {rotamer_idx} out of range for {chain}:{resname}{resnum}")
            continue
        
        try:
            all_rotamers[rotamer_idx].apply()
            if verbose:
                print(f"  Applied rotamer {rotamer_idx} to {chain}:{resname}{resnum}")
        except Exception as e:
            if verbose:
                print(f"  Warning: Failed to apply rotamer to {chain}:{resname}{resnum}: {e}")
    
    return modified_struct


def finalize_structure(optimized_struct, binding_site_residues, initial_volume,
                       output_dir, sample_name, min_contacts=4, verbose=True):
    """
    Finalize the optimized structure: verify volume increase and save.
    
    Per design.md:
    - Run pyKVfinder one last time to verify volume increase
    - Check if single cavity or multiple
    - Save with appropriate filename
    
    Args:
        optimized_struct: Structure with best rotamers applied
        binding_site_residues: List of binding site residue tuples
        initial_volume: Initial volume before optimization
        output_dir: Directory to save output structures
        sample_name: Base name for output files
        min_contacts: Minimum binding site residues a cavity must contact (default 4)
        verbose: Whether to print progress
        
    Returns:
        Dictionary with finalization results
    """
    if verbose:
        print("\n" + "="*60)
        print("FINALIZING: Verifying optimized structure")
        print("="*60)
    
    # Measure final volume
    final_result = measure_binding_site_volume(
        optimized_struct,
        binding_site_residues,
        min_contacts=min_contacts,
        verbose=verbose
    )
    
    if not final_result['success']:
        print("ERROR: No cavities found in optimized structure!")
        return {'success': False, 'error': 'No cavities in final structure'}
    
    final_volume = final_result['total_volume']
    num_cavities = final_result['num_cavities']
    volume_change = final_volume - initial_volume
    
    if verbose:
        print(f"\n=== FINAL RESULTS ===")
        print(f"Initial volume: {initial_volume:.2f} Å^3")
        print(f"Final volume: {final_volume:.2f} Å^3")
        print(f"Volume change: {volume_change:+.2f} Å^3")
        print(f"Number of binding site cavities: {num_cavities}")
    
    # Determine output filename
    os.makedirs(output_dir, exist_ok=True)
    
    if num_cavities > 1:
        output_filename = f"{sample_name}_optimized_still_multiple_pockets.pdb"
        if verbose:
            print(f"\nWARNING: Multiple cavities still detected ({num_cavities})")
    else:
        output_filename = f"{sample_name}_optimized.pdb"
        if verbose:
            print(f"\nSUCCESS: Single contiguous cavity achieved")
    
    output_path = os.path.join(output_dir, output_filename)
    
    # Save structure
    with structure.StructureWriter(output_path) as writer:
        writer.append(optimized_struct)
    
    if verbose:
        print(f"\nOptimized structure saved: {output_path}")
    
    # Verify volume increased
    if volume_change <= 0:
        print(f"\nWARNING: Volume did not increase (change: {volume_change:+.2f} Å^3)")
        print("The rotamer optimization may not have been effective for this structure.")
    
    print("\ndone")
    
    return {
        'success': True,
        'initial_volume': initial_volume,
        'final_volume': final_volume,
        'volume_change': volume_change,
        'num_cavities': num_cavities,
        'output_path': output_path,
        'single_cavity': num_cavities == 1
    }
