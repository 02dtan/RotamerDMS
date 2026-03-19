"""
Mutation phase module for RotamerDMS.

Mutates binding site residues to Alanine one by one and measures
the impact on cavity volume to identify residues for rotamer sampling.
"""

import copy
import json
import os
from datetime import datetime

from schrodinger import structure
from schrodinger.structutils import build

from .cavity import measure_binding_site_volume


def mutate_residue_to_ala(struct, chain, resnum):
    """
    Mutate a single residue to Alanine.
    
    Args:
        struct: Schrodinger Structure object (will be modified in place)
        chain: Chain identifier
        resnum: Residue number
        
    Returns:
        True if mutation was successful, False otherwise
    """
    target_res = None
    for res in struct.residue:
        if res.chain == chain and res.resnum == resnum:
            target_res = res
            break
    
    if target_res is None:
        print(f"Warning: Residue {chain}:{resnum} not found for mutation.")
        return False
    
    # Skip if already Alanine or Glycine (no sidechain to remove)
    resname = target_res.pdbres.strip()
    if resname in ['ALA', 'GLY']:
        print(f"Skipping {chain}:{resname}{resnum} - already Ala/Gly")
        return False
    
    try:
        # Mutate using the first atom of the residue
        build.mutate(struct, target_res.atom[1], 'ALA')
        return True
    except Exception as e:
        print(f"Error mutating {chain}:{resname}{resnum}: {e}")
        return False


def run_mutation_phase(sample_struct, binding_site_residues, min_contacts=4, checkpoint_dir=None, verbose=True):
    """
    Run the mutation phase: mutate each binding site residue to Ala and measure volume impact.
    
    Per design.md:
    - Mutate each residue to Ala one by one
    - Measure cavity volume change
    - Keep residues that INCREASE volume when mutated
    - Remove residues that have NO effect on volume
    - Sort remaining by volume increase (descending)
    
    Args:
        sample_struct: Schrodinger Structure object (original sample)
        binding_site_residues: List of (chain, resnum, resname) tuples
        min_contacts: Minimum binding site residues a cavity must contact (default 4)
        checkpoint_dir: Optional directory to save checkpoints
        verbose: Whether to print detailed progress
        
    Returns:
        Dictionary containing:
            - 'sorted_residues': List of (chain, resnum, resname) sorted by volume impact
            - 'volume_changes': Dict mapping residue tuple to volume change (Å³)
            - 'initial_volume': Initial binding site volume
            - 'removed_residues': Residues removed due to no volume effect
    """
    if verbose:
        print("\n" + "="*60)
        print("MUTATION PHASE: Identifying volume-affecting residues")
        print("="*60)
    
    # Measure initial volume
    if verbose:
        print("\nMeasuring initial binding site volume...")
    
    initial_result = measure_binding_site_volume(
        sample_struct, 
        binding_site_residues,
        min_contacts=min_contacts,
        verbose=verbose
    )
    
    if not initial_result['success']:
        raise RuntimeError("Failed to measure initial cavity volume. No suitable cavities found.")
    
    initial_volume = initial_result['total_volume']
    
    if verbose:
        print(f"\nInitial binding site volume: {initial_volume:.2f} Å^3")
        print(f"Number of cavities: {initial_result['num_cavities']}")
    
    # Track results
    volume_changes = {}
    removed_residues = []
    
    if verbose:
        print(f"\nTesting {len(binding_site_residues)} residues for volume impact...")
        print("-" * 40)
    
    for i, (chain, resnum, resname) in enumerate(binding_site_residues):
        if verbose:
            print(f"\n[{i+1}/{len(binding_site_residues)}] Testing {chain}:{resname}{resnum}...")
        
        # Skip Ala and Gly
        if resname in ['ALA', 'GLY']:
            if verbose:
                print(f"  Skipping - already Ala/Gly (no sidechain)")
            removed_residues.append((chain, resnum, resname))
            continue
        
        # Create a copy of the structure for this mutation
        mutant_struct = sample_struct.copy()
        
        # Perform mutation
        success = mutate_residue_to_ala(mutant_struct, chain, resnum)
        
        if not success:
            if verbose:
                print(f"  Mutation failed - removing from list")
            removed_residues.append((chain, resnum, resname))
            continue
        
        # Measure new volume
        mutant_result = measure_binding_site_volume(
            mutant_struct, 
            binding_site_residues,
            min_contacts=min_contacts,
            verbose=False
        )
        
        if not mutant_result['success']:
            if verbose:
                print(f"  No cavities found after mutation - removing from list")
            removed_residues.append((chain, resnum, resname))
            continue
        
        mutant_volume = mutant_result['total_volume']
        delta_volume = mutant_volume - initial_volume
        
        if verbose:
            print(f"  Volume after Ala mutation: {mutant_volume:.2f} Å^3")
            print(f"  Volume change: {delta_volume:+.2f} Å^3")
        
        # Per design.md: keep only if mutation INCREASES volume
        if delta_volume > 0:
            volume_changes[(chain, resnum, resname)] = delta_volume
            if verbose:
                print(f"  -> KEPT (volume increased)")
        else:
            removed_residues.append((chain, resnum, resname))
            if verbose:
                print(f"  -> REMOVED (no volume increase)")
    
    # Sort residues by volume change (descending - largest increase first)
    sorted_residues = sorted(
        volume_changes.keys(),
        key=lambda r: volume_changes[r],
        reverse=True
    )
    
    # Print summary
    if verbose:
        print("\n" + "="*60)
        print("MUTATION PHASE COMPLETE")
        print("="*60)
        print(f"\nResidues kept for rotamer sampling: {len(sorted_residues)}")
        print(f"Residues removed (no volume benefit): {len(removed_residues)}")
        
        if removed_residues:
            print("\nRemoved residues:")
            for chain, resnum, resname in removed_residues:
                print(f"  - {chain}:{resname}{resnum}")
        
        print("\nSorted residues (by volume increase upon Ala mutation):")
        for rank, (chain, resnum, resname) in enumerate(sorted_residues, 1):
            delta = volume_changes[(chain, resnum, resname)]
            print(f"  {rank}. {chain}:{resname}{resnum} -> +{delta:.2f} Å^3")
    
    result = {
        'sorted_residues': sorted_residues,
        'volume_changes': {str(k): v for k, v in volume_changes.items()},
        'initial_volume': initial_volume,
        'removed_residues': removed_residues,
        'initial_num_cavities': initial_result['num_cavities']
    }
    
    # Save checkpoint if directory provided
    if checkpoint_dir:
        save_mutation_checkpoint(result, checkpoint_dir)
    
    return result


def save_mutation_checkpoint(result, checkpoint_dir):
    """
    Save mutation phase results to checkpoint directory.
    
    Args:
        result: Result dictionary from run_mutation_phase
        checkpoint_dir: Directory to save checkpoint
    """
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Convert tuples to strings for JSON serialization
    checkpoint_data = {
        'timestamp': timestamp,
        'phase': 'mutation',
        'initial_volume': result['initial_volume'],
        'initial_num_cavities': result['initial_num_cavities'],
        'sorted_residues': [
            {'chain': r[0], 'resnum': r[1], 'resname': r[2]} 
            for r in result['sorted_residues']
        ],
        'volume_changes': [
            {
                'chain': r[0], 
                'resnum': r[1], 
                'resname': r[2],
                'delta_volume': result['volume_changes'][str(r)]
            }
            for r in result['sorted_residues']
        ],
        'removed_residues': [
            {'chain': r[0], 'resnum': r[1], 'resname': r[2]} 
            for r in result['removed_residues']
        ]
    }
    
    checkpoint_path = os.path.join(checkpoint_dir, f'mutation_phase_{timestamp}.json')
    
    with open(checkpoint_path, 'w') as f:
        json.dump(checkpoint_data, f, indent=2)
    
    print(f"\nCheckpoint saved: {checkpoint_path}")
    
    return checkpoint_path


def select_top_residues(sorted_residues, volume_changes, top_percent=100):
    """
    Select top percentage of residues for rotamer sampling.
    
    Per design.md: Give user option to select percentage threshold of
    top-volume-affecting residues.
    
    Args:
        sorted_residues: List of residue tuples sorted by volume impact
        volume_changes: Dict mapping residue tuple to volume change
        top_percent: Percentage of top residues to keep (1-100)
        
    Returns:
        Filtered list of residue tuples
    """
    if top_percent >= 100:
        return sorted_residues
    
    if top_percent <= 0:
        return []
    
    n_keep = max(1, int(len(sorted_residues) * top_percent / 100))
    selected = sorted_residues[:n_keep]
    
    print(f"\nSelected top {top_percent}% of residues ({n_keep} of {len(sorted_residues)}):")
    for chain, resnum, resname in selected:
        key = (chain, resnum, resname)
        delta = volume_changes.get(str(key), volume_changes.get(key, 0))
        print(f"  - {chain}:{resname}{resnum} (+{delta:.2f} Å^3)")
    
    return selected
