#!/usr/bin/env python3
"""
Pareto Front Analysis Script for RotamerDMS

Analyzes rotamer sampling checkpoint files to visualize and export
the Pareto front tradeoff between pocket volume and WCA clash energy.

Usage:
    python analyze_pareto.py --checkpoint <checkpoint.json> [--output <output_dir>]
"""

import argparse
import json
import os
import sys

# Add src to path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

from src.scoring import (
    compute_pareto_front_2d,
    compute_pareto_front_3d,
    rank_rotamers_by_joint_score
)


def load_checkpoint(checkpoint_path):
    """Load rotamer sampling checkpoint file."""
    with open(checkpoint_path, 'r') as f:
        return json.load(f)


def extract_all_rotamers_from_checkpoint(checkpoint_data):
    """
    Extract all rotamer results from checkpoint for per-residue Pareto analysis.
    
    Returns a dict mapping residue identifiers to lists of all rotamer data.
    """
    residue_details = checkpoint_data.get('residue_details', [])
    
    residue_rotamers = {}
    residue_summaries = []
    
    for res in residue_details:
        residue_id = f"{res['chain']}:{res['resname']}{res['resnum']}"
        
        # Extract all rotamers if available (new checkpoint format)
        all_rotamers = res.get('all_rotamers', [])
        
        if all_rotamers:
            # New format: full rotamer data available
            rotamers_with_id = []
            for rot in all_rotamers:
                rotamers_with_id.append({
                    'residue': residue_id,
                    'index': rot['index'],
                    'volume_change': rot['volume_change'],
                    'percentage': rot['percentage'],
                    'delta_g': rot['delta_g'],
                    'wca_energy': rot['wca_energy'],
                    'joint_score': rot.get('joint_score')
                })
            residue_rotamers[residue_id] = rotamers_with_id
        
        # Always store summary for backward compatibility
        residue_summaries.append({
            'residue': residue_id,
            'chain': res['chain'],
            'resnum': res['resnum'],
            'resname': res['resname'],
            'best_rotamer_idx': res['best_rotamer_idx'],
            'best_volume_change': res['best_volume_change'],
            'best_delta_g': res.get('best_delta_g'),
            'best_wca_energy': res.get('best_wca_energy'),
            'best_joint_score': res.get('best_joint_score'),
            'num_pareto_optimal': res.get('num_pareto_optimal', 0),
            'pareto_front_indices': res.get('pareto_front_indices', []),
            'num_rotamers': len(all_rotamers) if all_rotamers else res.get('num_rotamers_tested', 0)
        })
    
    return residue_rotamers, residue_summaries


def analyze_per_residue_pareto(residue_rotamers, verbose=True):
    """
    Compute 3D Pareto fronts (volume, deltaG, WCA) for each residue.
    
    Args:
        residue_rotamers: Dict mapping residue IDs to lists of rotamer data
        verbose: Print detailed analysis
        
    Returns:
        Dict mapping residue IDs to their Pareto analysis results
    """
    per_residue_results = {}
    
    for residue_id, rotamers in residue_rotamers.items():
        # Filter to rotamers with valid data for all 3 metrics
        valid_rotamers = [
            r for r in rotamers 
            if r.get('volume_change') is not None 
            and r.get('delta_g') is not None 
            and r.get('wca_energy') is not None
        ]
        
        if not valid_rotamers:
            # Fall back to 2D if no WCA data
            valid_rotamers = [
                r for r in rotamers 
                if r.get('volume_change') is not None 
                and r.get('delta_g') is not None
            ]
            if valid_rotamers:
                pareto_3d = compute_pareto_front_2d(valid_rotamers, 'volume_change', 'delta_g')
                has_wca = False
            else:
                continue
        else:
            # Compute 3D Pareto front: maximize volume, minimize deltaG and WCA
            pareto_3d = compute_pareto_front_3d(
                valid_rotamers,
                maximize_keys=['volume_change'],
                minimize_keys=['delta_g', 'wca_energy']
            )
            has_wca = True
        
        # Sort Pareto front by volume (descending)
        pareto_3d = sorted(pareto_3d, key=lambda x: x['volume_change'], reverse=True)
        
        per_residue_results[residue_id] = {
            'pareto_front': pareto_3d,
            'num_pareto_optimal': len(pareto_3d),
            'num_rotamers_total': len(rotamers),
            'has_wca': has_wca,
            'all_rotamers': valid_rotamers
        }
        
        if verbose:
            print(f"\n{residue_id}: {len(pareto_3d)} Pareto-optimal rotamers out of {len(rotamers)}")
            print(f"  {'Idx':<6}{'ΔVolume':>10}{'ΔG':>10}{'WCA':>12}{'Score':>10}")
            print(f"  {'-'*48}")
            for r in pareto_3d:
                dg_str = f"{r['delta_g']:.2f}" if r.get('delta_g') is not None else "N/A"
                wca_str = f"{r['wca_energy']:.2f}" if r.get('wca_energy') is not None else "N/A"
                score_str = f"{r['joint_score']:.3f}" if r.get('joint_score') is not None else "N/A"
                print(f"  {r['index']:<6}{r['volume_change']:>+10.2f}{dg_str:>10}{wca_str:>12}{score_str:>10}")
    
    return per_residue_results


def analyze_checkpoint_pareto(checkpoint_path, output_dir=None, verbose=True):
    """
    Analyze Pareto front from checkpoint file.
    
    Args:
        checkpoint_path: Path to rotamer sampling checkpoint JSON
        output_dir: Optional directory for output files
        verbose: Print detailed analysis
        
    Returns:
        Analysis results dictionary
    """
    if verbose:
        print(f"\nLoading checkpoint: {checkpoint_path}")
    
    checkpoint = load_checkpoint(checkpoint_path)
    
    if verbose:
        print(f"Timestamp: {checkpoint.get('timestamp', 'N/A')}")
        print(f"Phase: {checkpoint.get('phase', 'N/A')}")
        print(f"Estimated final volume: {checkpoint.get('estimated_final_volume', 'N/A'):.2f} Å³")
        print(f"Volume increase: {checkpoint.get('volume_increase', 'N/A'):.2f} Å³")
    
    # Extract rotamer data
    residue_rotamers, residue_summaries = extract_all_rotamers_from_checkpoint(checkpoint)
    
    if verbose:
        print(f"\nResidues analyzed: {len(residue_summaries)}")
    
    # Check if we have full rotamer data (new checkpoint format)
    has_full_data = len(residue_rotamers) > 0
    
    if has_full_data:
        if verbose:
            print(f"Full rotamer data available: {len(residue_rotamers)} residues")
            print(f"\n{'='*60}")
            print("PER-RESIDUE 3D PARETO ANALYSIS (Volume ↑, ΔG ↓, WCA ↓)")
            print(f"{'='*60}")
        
        # Compute per-residue 3D Pareto fronts
        per_residue_results = analyze_per_residue_pareto(residue_rotamers, verbose=verbose)
    else:
        if verbose:
            print("Note: Checkpoint uses old format (best rotamer only).")
            print("Re-run rotamer sampling to get full Pareto analysis.")
        per_residue_results = {}
    
    # Also show cross-residue analysis using best rotamers
    valid_summaries = [r for r in residue_summaries if r.get('best_wca_energy') is not None]
    
    if valid_summaries and verbose:
        print(f"\n{'='*60}")
        print("CROSS-RESIDUE SUMMARY (Best Rotamer per Residue)")
        print(f"{'='*60}")
        print(f"\n{'Residue':<15}{'Best Idx':>10}{'ΔVolume':>12}{'ΔG':>10}{'WCA':>12}")
        print("-" * 59)
        for r in sorted(valid_summaries, key=lambda x: x['best_volume_change'], reverse=True):
            dg_str = f"{r['best_delta_g']:.2f}" if r.get('best_delta_g') is not None else "N/A"
            print(f"{r['residue']:<15}{r['best_rotamer_idx']:>10}{r['best_volume_change']:>+12.2f}{dg_str:>10}{r['best_wca_energy']:>12.2f}")
    
    # Export results if output_dir specified
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        
        if has_full_data:
            # Export per-residue Pareto fronts as CSV
            pareto_csv_path = os.path.join(output_dir, 'pareto_fronts_per_residue.csv')
            with open(pareto_csv_path, 'w') as f:
                f.write("residue,rotamer_idx,volume_change,delta_g,wca_energy,joint_score,on_pareto_front\n")
                for residue_id, res_data in per_residue_results.items():
                    pareto_indices = set(r['index'] for r in res_data['pareto_front'])
                    for rot in res_data['all_rotamers']:
                        dg = f"{rot['delta_g']:.4f}" if rot.get('delta_g') is not None else ""
                        wca = f"{rot['wca_energy']:.4f}" if rot.get('wca_energy') is not None else ""
                        js = f"{rot['joint_score']:.4f}" if rot.get('joint_score') is not None else ""
                        on_pareto = 'yes' if rot['index'] in pareto_indices else 'no'
                        f.write(f"{residue_id},{rot['index']},{rot['volume_change']:.4f},{dg},{wca},{js},{on_pareto}\n")
            
            if verbose:
                print(f"\nPer-residue Pareto data exported to: {pareto_csv_path}")
        
        # Export residue summaries
        summary_csv_path = os.path.join(output_dir, 'residue_summaries.csv')
        with open(summary_csv_path, 'w') as f:
            f.write("residue,best_rotamer_idx,volume_change,delta_g,wca_energy,joint_score,num_pareto_optimal,num_rotamers\n")
            for r in sorted(residue_summaries, key=lambda x: x['best_volume_change'], reverse=True):
                dg = f"{r['best_delta_g']:.4f}" if r.get('best_delta_g') is not None else ""
                wca = f"{r['best_wca_energy']:.4f}" if r.get('best_wca_energy') is not None else ""
                js = f"{r['best_joint_score']:.4f}" if r.get('best_joint_score') is not None else ""
                f.write(f"{r['residue']},{r['best_rotamer_idx']},{r['best_volume_change']:.4f},{dg},{wca},{js},{r['num_pareto_optimal']},{r['num_rotamers']}\n")
        
        if verbose:
            print(f"Residue summaries exported to: {summary_csv_path}")
        
        # Export summary JSON
        summary = {
            'checkpoint_file': checkpoint_path,
            'num_residues': len(residue_summaries),
            'has_full_rotamer_data': has_full_data,
            'per_residue_pareto': {
                residue_id: {
                    'num_pareto_optimal': res_data['num_pareto_optimal'],
                    'num_rotamers_total': res_data['num_rotamers_total'],
                    'pareto_front': [
                        {
                            'index': r['index'],
                            'volume_change': r['volume_change'],
                            'delta_g': r.get('delta_g'),
                            'wca_energy': r.get('wca_energy'),
                            'joint_score': r.get('joint_score')
                        }
                        for r in res_data['pareto_front']
                    ]
                }
                for residue_id, res_data in per_residue_results.items()
            } if has_full_data else {},
            'residue_summaries': [
                {
                    'residue': r['residue'],
                    'best_rotamer_idx': r['best_rotamer_idx'],
                    'best_volume_change': r['best_volume_change'],
                    'best_delta_g': r.get('best_delta_g'),
                    'best_wca_energy': r.get('best_wca_energy'),
                    'best_joint_score': r.get('best_joint_score')
                }
                for r in residue_summaries
            ]
        }
        
        summary_path = os.path.join(output_dir, 'pareto_analysis_summary.json')
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        
        if verbose:
            print(f"Summary exported to: {summary_path}")
    
    return {
        'success': True,
        'per_residue_results': per_residue_results,
        'residue_summaries': residue_summaries,
        'has_full_data': has_full_data
    }


def main():
    parser = argparse.ArgumentParser(
        description='Analyze Pareto front from RotamerDMS checkpoint',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--checkpoint', '-c',
        required=True,
        help='Path to rotamer sampling checkpoint JSON file'
    )
    
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output directory for analysis results (CSV, JSON)'
    )
    
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Reduce output verbosity'
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.checkpoint):
        print(f"Error: Checkpoint file not found: {args.checkpoint}")
        sys.exit(1)
    
    result = analyze_checkpoint_pareto(
        args.checkpoint,
        output_dir=args.output,
        verbose=not args.quiet
    )
    
    if not result['success']:
        print(f"Analysis failed: {result.get('error', 'Unknown error')}")
        sys.exit(1)
    
    print("\nPareto analysis complete.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
