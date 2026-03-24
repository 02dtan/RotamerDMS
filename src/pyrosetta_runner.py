#!/usr/bin/env python3
"""
Standalone PyRosetta runner script.
Runs in conda's Python (with PyRosetta) and communicates via JSON.

This script performs energy minimization on a structure with a frozen target
residue, then extracts the fa_rep (repulsive) energy for that residue.

The minimization relaxes surrounding residues to relieve steric clashes while
keeping the target rotamer fixed, resulting in more physically realistic
structures with lower clash energies.

Usage:
    echo '{"pdb_file": "path.pdb", "chain": "A", "resnum": 102}' | python pyrosetta_runner.py
"""

import sys
import json
import tempfile
import os

import pyrosetta
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.minimization_packing import MinMover


# Initialize PyRosetta once at module load
pyrosetta.init('-mute all', silent=True)


def find_pose_resnum(pose, chain, resnum):
    """
    Find the pose residue number for a given chain and PDB residue number.
    
    Args:
        pose: PyRosetta Pose object
        chain: Chain identifier (e.g., 'A')
        resnum: PDB residue number
        
    Returns:
        Pose residue number (1-indexed) or None if not found
    """
    pdb_info = pose.pdb_info()
    
    if pdb_info:
        # Try exact match first
        for i in range(1, pose.total_residue() + 1):
            if pdb_info.chain(i) == chain and pdb_info.number(i) == resnum:
                return i
        
        # Fallback: try matching just residue number
        for i in range(1, pose.total_residue() + 1):
            if pdb_info.number(i) == resnum:
                return i
    
    return None


def minimize_with_frozen_residue(pdb_file, chain, resnum, verbose=False):
    """
    Minimize structure with ref2015 while keeping target residue frozen.
    
    This allows surrounding residues to relax and relieve steric clashes
    caused by the rotamer change, while preserving the target rotamer state.
    
    Args:
        pdb_file: Path to PDB file
        chain: Chain identifier for target residue
        resnum: Residue number of target (to freeze)
        verbose: Print details
        
    Returns:
        Dictionary with:
            - success: bool
            - fa_rep_energy: float - repulsive energy for target residue after minimization
            - total_score: float - total ref2015 score after minimization
            - minimized_pdb: str - PDB content of minimized structure
            - error: str (only if success=False)
    """
    try:
        # Load structure
        pose = pyrosetta.pose_from_pdb(pdb_file)
        
        # Find target residue in pose numbering
        target_pose_resnum = find_pose_resnum(pose, chain, resnum)
        
        if target_pose_resnum is None:
            return {
                'success': False,
                'error': f'Could not find residue {chain}:{resnum}'
            }
        
        # Create ref2015 score function
        sfxn = pyrosetta.create_score_function('ref2015')
        
        # Create MoveMap: allow all residues to move EXCEPT the target
        movemap = MoveMap()
        
        # Enable backbone and chi for all residues by default
        movemap.set_bb(True)
        movemap.set_chi(True)
        
        # Freeze the target residue (both backbone and chi)
        movemap.set_bb(target_pose_resnum, False)
        movemap.set_chi(target_pose_resnum, False)
        
        # Create MinMover with ref2015
        minmover = MinMover()
        minmover.movemap(movemap)
        minmover.score_function(sfxn)
        minmover.min_type('lbfgs_armijo_nonmonotone')
        minmover.tolerance(0.01)
        minmover.max_iter(200)
        
        # Score before minimization (optional, for debugging)
        score_before = sfxn(pose)
        
        # Run minimization
        minmover.apply(pose)
        
        # Score after minimization
        total_score = sfxn(pose)
        
        # Extract energy terms for the target residue
        energies = pose.energies()
        residue_energies = energies.residue_total_energies(target_pose_resnum)
        fa_rep_energy = residue_energies[ScoreType.fa_rep]
        fa_dun_energy = residue_energies[ScoreType.fa_dun]
        
        # Write minimized structure to a temp file and read content
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False, mode='w') as f:
            temp_path = f.name
        
        pose.dump_pdb(temp_path)
        
        with open(temp_path, 'r') as f:
            minimized_pdb = f.read()
        
        os.remove(temp_path)
        
        return {
            'success': True,
            'fa_rep_energy': float(fa_rep_energy),
            'fa_dun_energy': float(fa_dun_energy),
            'total_score': float(total_score),
            'score_before': float(score_before),
            'minimized_pdb': minimized_pdb,
        }
        
    except Exception as e:
        return {
            'success': False,
            'error': str(e)
        }


def main():
    """Read JSON from stdin, write JSON to stdout."""
    try:
        input_data = json.load(sys.stdin)
        
        pdb_file = input_data.get('pdb_file')
        chain = input_data.get('chain')
        resnum = input_data.get('resnum')
        verbose = input_data.get('verbose', False)
        
        if not pdb_file:
            result = {'success': False, 'error': 'pdb_file is required'}
        elif chain is None:
            result = {'success': False, 'error': 'chain is required'}
        elif resnum is None:
            result = {'success': False, 'error': 'resnum is required'}
        else:
            result = minimize_with_frozen_residue(pdb_file, chain, resnum, verbose)
        
        json.dump(result, sys.stdout)
        
    except Exception as e:
        json.dump({'success': False, 'error': str(e)}, sys.stdout)
        sys.exit(1)


if __name__ == '__main__':
    main()
