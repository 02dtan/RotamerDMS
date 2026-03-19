#!/usr/bin/env python3
"""
RotamerDMS - Main Runner Script

Maximizes binding pocket volume through rotamer sampling.

Usage:
    $SCHRODINGER/run python3 run_rotamer_dms.py \
        --reference <reference_pdb> \
        --sample <sample_pdb> \
        --ligand <ligand_code> \
        [--distance <distance_threshold>] \
        [--top-percent <percentage>]
"""

import argparse
import os
import sys
from datetime import datetime

# Add src to path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from src.preparation import (
    load_structure,
    align_sequences,
    get_binding_site_residues,
    map_residues_to_sample
)
from src.mutation import (
    run_mutation_phase,
    select_top_residues
)
from src.rotamer_sampling import (
    run_rotamer_sampling_phase,
    finalize_structure
)


def parse_args():
    parser = argparse.ArgumentParser(
        description='RotamerDMS: Maximize binding pocket volume through rotamer sampling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage with defaults
    $SCHRODINGER/run python3 run_rotamer_dms.py \\
        --reference data/4W59_n-hexylbenzene.pdb \\
        --sample data/s0067_rec.pdb \\
        --ligand 3GZ

    # With custom distance threshold and top residue percentage
    $SCHRODINGER/run python3 run_rotamer_dms.py \\
        --reference data/4W59_n-hexylbenzene.pdb \\
        --sample data/s0067_rec.pdb \\
        --ligand 3GZ \\
        --distance 5.0 \\
        --top-percent 50
        """
    )
    
    parser.add_argument(
        '--reference', '-r',
        required=True,
        help='Path to experimental reference structure with bound ligand'
    )
    
    parser.add_argument(
        '--sample', '-s',
        required=True,
        help='Path to sample structure to optimize'
    )
    
    parser.add_argument(
        '--ligand', '-l',
        required=True,
        help='3-letter residue code of the ligand in reference structure'
    )
    
    parser.add_argument(
        '--distance', '-d',
        type=float,
        default=7.0,
        help='Distance threshold (Å) for binding site residue selection (default: 7.0)'
    )
    
    parser.add_argument(
        '--top-percent', '-t',
        type=float,
        default=100.0,
        help='Percentage of top volume-affecting residues to sample (default: 100)'
    )
    
    parser.add_argument(
        '--min-contacts', '-m',
        type=int,
        default=4,
        help='Minimum binding site residues a cavity must contact to be counted (default: 4)'
    )
    
    parser.add_argument(
        '--output-dir', '-o',
        default=None,
        help='Output directory (default: ./output)'
    )
    
    parser.add_argument(
        '--checkpoint-dir', '-c',
        default=None,
        help='Checkpoint directory (default: ./checkpoints)'
    )
    
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Reduce output verbosity'
    )
    
    return parser.parse_args()


def main():
    args = parse_args()
    verbose = not args.quiet
    
    # Set up directories
    base_dir = script_dir
    output_dir = args.output_dir or os.path.join(base_dir, 'output')
    checkpoint_dir = args.checkpoint_dir or os.path.join(base_dir, 'checkpoints')
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    # Extract sample name for output files
    sample_name = os.path.splitext(os.path.basename(args.sample))[0]
    
    print("="*70)
    print("RotamerDMS - Binding Pocket Volume Maximization")
    print("="*70)
    print(f"\nReference structure: {args.reference}")
    print(f"Sample structure: {args.sample}")
    print(f"Ligand code: {args.ligand}")
    print(f"Distance threshold: {args.distance} Å")
    print(f"Top residue percentage: {args.top_percent}%")
    print(f"Min cavity contacts: {args.min_contacts}")
    print(f"Output directory: {output_dir}")
    print(f"Checkpoint directory: {checkpoint_dir}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # =========================================================================
    # PHASE 0: PREPARATION
    # =========================================================================
    print("\n" + "="*70)
    print("PHASE 0: PREPARATION")
    print("="*70)
    
    # Load structures
    print("\nLoading structures...")
    try:
        ref_struct = load_structure(args.reference)
        print(f"  Reference structure loaded: {ref_struct.atom_total} atoms")
    except IOError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    try:
        sample_struct = load_structure(args.sample)
        print(f"  Sample structure loaded: {sample_struct.atom_total} atoms")
    except IOError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    # Align sequences between reference and sample
    print("\nAligning sequences...")
    try:
        common_residues, ref_only, sample_only = align_sequences(ref_struct, sample_struct)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    # Identify binding site residues (only from common residues)
    print(f"\nIdentifying binding site residues (within {args.distance} Å of ligand {args.ligand})...")
    try:
        binding_site_residues = get_binding_site_residues(
            ref_struct, 
            args.ligand, 
            distance_threshold=args.distance
        )
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    if not binding_site_residues:
        print("ERROR: No binding site residues found!")
        sys.exit(1)
    
    # Filter binding site to only include common residues
    common_set = set((r[0], r[1]) for r in common_residues)
    filtered_binding_site = [
        r for r in binding_site_residues 
        if (r[0], r[1]) in common_set
    ]
    
    excluded_from_binding = len(binding_site_residues) - len(filtered_binding_site)
    if excluded_from_binding > 0:
        print(f"  Note: {excluded_from_binding} binding site residue(s) excluded (not in both structures)")
    
    binding_site_residues = filtered_binding_site
    
    # Map to sample structure
    print("\nMapping binding site residues to sample structure...")
    mapped_residues = map_residues_to_sample(binding_site_residues, sample_struct)
    print(f"  {len(mapped_residues)} residues mapped successfully")
    
    # =========================================================================
    # PHASE 1: MUTATION
    # =========================================================================
    print("\n" + "="*70)
    print("PHASE 1: MUTATION (Identifying volume-affecting residues)")
    print("="*70)
    
    mutation_result = run_mutation_phase(
        sample_struct,
        mapped_residues,
        min_contacts=args.min_contacts,
        checkpoint_dir=checkpoint_dir,
        verbose=verbose
    )
    
    initial_volume = mutation_result['initial_volume']
    sorted_residues = mutation_result['sorted_residues']
    volume_changes = mutation_result['volume_changes']
    
    if not sorted_residues:
        print("\nWARNING: No residues found that increase volume upon Ala mutation.")
        print("The binding pocket may already be maximally open or sample structure unsuitable.")
        sys.exit(0)
    
    # Apply top percentage filter
    if args.top_percent < 100:
        sorted_residues = select_top_residues(
            sorted_residues, 
            volume_changes, 
            args.top_percent
        )
    
    # =========================================================================
    # PHASE 2: ROTAMER SAMPLING
    # =========================================================================
    print("\n" + "="*70)
    print("PHASE 2: ROTAMER SAMPLING (Finding optimal states)")
    print("="*70)
    
    rotamer_result, optimized_struct = run_rotamer_sampling_phase(
        sample_struct,
        sorted_residues,
        mapped_residues,
        initial_volume,
        min_contacts=args.min_contacts,
        checkpoint_dir=checkpoint_dir,
        verbose=verbose
    )
    
    # =========================================================================
    # PHASE 3: FINALIZATION
    # =========================================================================
    print("\n" + "="*70)
    print("PHASE 3: FINALIZATION")
    print("="*70)
    
    final_result = finalize_structure(
        optimized_struct,
        mapped_residues,
        initial_volume,
        output_dir,
        sample_name,
        min_contacts=args.min_contacts,
        verbose=verbose
    )
    
    # =========================================================================
    # SUMMARY
    # =========================================================================
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\nInitial binding site volume: {initial_volume:.2f} Å^3")
    print(f"Final binding site volume: {final_result['final_volume']:.2f} Å^3")
    print(f"Volume change: {final_result['volume_change']:+.2f} Å^3")
    print(f"Final number of cavities: {final_result['num_cavities']}")
    print(f"Single contiguous cavity: {'Yes' if final_result['single_cavity'] else 'No'}")
    print(f"\nOutput saved to: {final_result['output_path']}")
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
