"""
Rotamer sampling phase module for RotamerDMS.

Samples rotamer states from Schrodinger's backbone-dependent rotamer library
and selects states that maximize binding pocket volume.
"""

import copy
import json
import os
from datetime import datetime

from schrodinger import structure
from schrodinger.protein import rotamers

from .cavity import measure_binding_site_volume


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
                                 initial_volume, min_contacts=4, verbose=True):
    """
    Sample all rotamer states for a single residue and rank by volume impact.
    
    Args:
        struct: Schrodinger Structure object (working copy)
        chain: Chain identifier
        resnum: Residue number
        resname: Residue name (3-letter code)
        binding_site_residues: List of binding site residue tuples
        initial_volume: Current binding site volume before rotamer changes
        verbose: Whether to print progress
        
    Returns:
        Dictionary containing:
            - 'rotamer_states': List of rotamer indices sorted by volume impact
            - 'volume_changes': List of volume changes (same order as rotamer_states)
            - 'percentages': List of rotamer percentages (same order)
            - 'best_rotamer_idx': Index of best rotamer (most volume increase)
            - 'best_volume_change': Volume change from best rotamer
            - 'native_percentage': Percentage of native rotamer state (if available)
    """
    if verbose:
        print(f"\n  Sampling rotamers for {chain}:{resname}{resnum}...")
    
    # Get rotamer library
    rotamer_lib, all_rotamers = get_rotamer_library(struct, chain, resnum)
    
    if rotamer_lib is None:
        return None
    
    if verbose:
        print(f"    Found {len(all_rotamers)} rotamer states in library")
    
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
        
        # Measure new volume
        result = measure_binding_site_volume(
            test_struct, 
            binding_site_residues,
            min_contacts=min_contacts,
            verbose=False
        )
        
        if not result['success']:
            continue
        
        new_volume = result['total_volume']
        delta_volume = new_volume - initial_volume
        
        # Get rotamer percentage (frequency proxy)
        try:
            pct = rotamer.percentage
            # Handle if percentage is a method vs property
            percentage = pct() if callable(pct) else float(pct)
        except (AttributeError, TypeError):
            percentage = 0.0
        
        rotamer_results.append({
            'index': rot_idx,
            'volume_change': delta_volume,
            'percentage': percentage,
            'new_volume': new_volume
        })
    
    if not rotamer_results:
        if verbose:
            print(f"    No valid rotamer states found")
        return None
    
    # Sort by volume change (descending - largest increase first)
    rotamer_results.sort(key=lambda x: x['volume_change'], reverse=True)
    
    # Extract sorted lists
    rotamer_states = [r['index'] for r in rotamer_results]
    volume_changes = [r['volume_change'] for r in rotamer_results]
    percentages = [r['percentage'] for r in rotamer_results]
    
    best = rotamer_results[0]
    
    if verbose:
        print(f"    Best rotamer: index {best['index']}, "
              f"ΔV = {best['volume_change']:+.2f} Å^3, "
              f"percentage = {best['percentage']:.1f}%")
        print(f"    Top 3 rotamers:")
        for i, r in enumerate(rotamer_results[:3]):
            print(f"      {i+1}. idx={r['index']}: ΔV={r['volume_change']:+.2f} Å^3, "
                  f"pct={r['percentage']:.1f}%")
    
    # Try to get native rotamer percentage (rotamer 1 is often closest to native)
    # TODO 3/13/2026: This is an incorrect way to compute the native rotamer percentage, though this is not in use
    # in the current version of the script. We need to compare dihedrals of the rotamer in the unaltered sample structure
    # to the dihedrals of each rotamer in Schrodinger's rotamer library, and select the closest discrete rotamer state
    # for use in population probability computations.
    native_percentage = percentages[rotamer_states.index(1)] if 1 in rotamer_states else None
    
    return {
        'rotamer_states': rotamer_states,
        'volume_changes': volume_changes,
        'percentages': percentages,
        'best_rotamer_idx': best['index'],
        'best_volume_change': best['volume_change'],
        'best_percentage': best['percentage'],
        'native_percentage': native_percentage,
        'all_results': rotamer_results
    }


def run_rotamer_sampling_phase(sample_struct, sorted_residues, binding_site_residues,
                                initial_volume, min_contacts=4, checkpoint_dir=None, verbose=True):
    """
    Run the rotamer sampling phase on sorted residues from mutation phase.
    
    Per design.md:
    - Process residues in order (highest volume impact first)
    - For each residue, sample all rotamer states
    - Record volume changes and percentages for each state
    - Sort states by volume contribution
    
    Args:
        sample_struct: Schrodinger Structure object (original sample)
        sorted_residues: List of residue tuples sorted by volume impact from mutation phase
        binding_site_residues: Full list of binding site residues
        initial_volume: Initial binding site volume
        min_contacts: Minimum binding site residues a cavity must contact (default 4)
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
            checkpoint_data['residue_details'].append({
                'chain': chain,
                'resnum': resnum,
                'resname': resname,
                'best_rotamer_idx': res_data['best_rotamer_idx'],
                'best_volume_change': res_data['best_volume_change'],
                'best_percentage': res_data['best_percentage'],
                'native_percentage': res_data['native_percentage'],
                'num_rotamers_tested': len(res_data['rotamer_states'])
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
