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
from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.numeric import xyzVector_double_t


# Initialize PyRosetta once at module load
pyrosetta.init('-mute all', silent=True)


# Constraint parameters for backbone preservation
CA_CONSTRAINT_SD = 0.5  # Standard deviation in Angstroms - tighter = less movement allowed
CA_CONSTRAINT_WEIGHT = 1.0  # Weight for coordinate_constraint term in score function


def strip_hydrogens_from_pdb(pdb_content):
    """
    Remove hydrogen atoms from PDB content.
    
    PyRosetta adds hydrogens during pose loading, but we want to return
    only heavy atoms to maintain consistency with the input structure.
    
    Args:
        pdb_content: String containing PDB file content
        
    Returns:
        PDB content string with hydrogen atoms removed
    """
    filtered_lines = []
    for line in pdb_content.split('\n'):
        # Keep non-ATOM/HETATM lines (headers, END, etc.)
        if not line.startswith(('ATOM', 'HETATM')):
            filtered_lines.append(line)
            continue
        
        # For ATOM/HETATM lines, check if it's a hydrogen
        # PDB format: columns 77-78 contain element symbol (right-justified)
        # Also check atom name in columns 13-16
        if len(line) >= 78:
            element = line[76:78].strip()
            if element == 'H':
                continue
        
        # Fallback: check atom name (columns 13-16)
        # Hydrogen atom names typically start with H or are like 1H, 2H, etc.
        atom_name = line[12:16].strip()
        if atom_name.startswith('H') or (len(atom_name) >= 2 and atom_name[0].isdigit() and atom_name[1] == 'H'):
            continue
        
        filtered_lines.append(line)
    
    return '\n'.join(filtered_lines)


def add_ca_coordinate_constraints(pose, exclude_resnum=None):
    """
    Add coordinate constraints to all CA atoms to preserve backbone positions.
    
    This applies a harmonic penalty to CA atoms that move from their original
    positions, limiting backbone movement during minimization while still
    allowing local relaxation to relieve clashes.
    
    Args:
        pose: PyRosetta Pose object (will be modified in place)
        exclude_resnum: Pose residue number to exclude from constraints (e.g., target residue)
        
    Returns:
        Number of constraints added
    """
    from pyrosetta.rosetta.core.pose import Pose
    
    constraints_added = 0
    
    # We need a fixed reference atom for CoordinateConstraint
    # Use the first residue's CA as the anchor (it's arbitrary for coordinate constraints)
    # Actually, we need to use a virtual root or the pose's root
    # For simplicity, use residue 1's CA as anchor
    anchor_resnum = 1
    if not pose.residue(anchor_resnum).has("CA"):
        # Find first residue with CA
        for i in range(1, pose.total_residue() + 1):
            if pose.residue(i).has("CA"):
                anchor_resnum = i
                break
    
    anchor_atom_id = AtomID(pose.residue(anchor_resnum).atom_index("CA"), anchor_resnum)
    
    for i in range(1, pose.total_residue() + 1):
        # Skip excluded residue
        if exclude_resnum is not None and i == exclude_resnum:
            continue
        
        residue = pose.residue(i)
        
        # Only constrain protein residues with CA
        if not residue.has("CA"):
            continue
        
        # Get CA atom index and current position
        ca_atom_idx = residue.atom_index("CA")
        ca_xyz = residue.xyz("CA")
        
        # Create AtomID for this CA
        ca_atom_id = AtomID(ca_atom_idx, i)
        
        # Create target position vector
        target_xyz = xyzVector_double_t(ca_xyz.x, ca_xyz.y, ca_xyz.z)
        
        # Create harmonic function: penalty increases as atom moves from target
        # HarmonicFunc(x0, sd) - x0 is the ideal value (0 for distance from target)
        # sd is standard deviation - smaller = tighter constraint
        harmonic_func = HarmonicFunc(0.0, CA_CONSTRAINT_SD)
        
        # Create coordinate constraint
        coord_constraint = CoordinateConstraint(
            ca_atom_id,      # Atom to constrain
            anchor_atom_id,  # Fixed reference atom
            target_xyz,      # Target position
            harmonic_func    # Penalty function
        )
        
        # Add constraint to pose
        pose.add_constraint(coord_constraint)
        constraints_added += 1
    
    return constraints_added


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
        
        # Enable coordinate_constraint term to enforce CA position constraints
        sfxn.set_weight(ScoreType.coordinate_constraint, CA_CONSTRAINT_WEIGHT)
        
        # Add CA coordinate constraints to preserve backbone positions
        # Exclude the target residue from constraints (it's already frozen)
        num_constraints = add_ca_coordinate_constraints(pose, exclude_resnum=target_pose_resnum)
        
        # Create MoveMap: chi-only minimization (no backbone movement)
        # This is appropriate for rotamer optimization where we want to preserve
        # the global backbone structure while allowing sidechains to relax
        movemap = MoveMap()
        
        # Disable backbone movement globally - this is the key change
        # Backbone should already be physically reasonable; only sidechains need to adjust
        movemap.set_bb(False)
        movemap.set_chi(True)
        
        # Also freeze the target residue's chi angles (preserve the applied rotamer)
        movemap.set_chi(target_pose_resnum, False)
        
        # Create MinMover with ref2015
        # Reduced aggressiveness: fewer iterations, higher tolerance
        minmover = MinMover()
        minmover.movemap(movemap)
        minmover.score_function(sfxn)
        minmover.min_type('lbfgs_armijo_nonmonotone')
        minmover.tolerance(0.1)   # Increased from 0.01 - converge earlier
        minmover.max_iter(50)     # Reduced from 200 - fewer iterations
        
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
            minimized_pdb_raw = f.read()
        
        os.remove(temp_path)
        
        # Strip hydrogens to maintain consistency with input structure
        minimized_pdb = strip_hydrogens_from_pdb(minimized_pdb_raw)
        
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


def minimize_full_structure(pdb_file, verbose=False):
    """
    Minimize entire structure with ref2015 (no frozen residues).
    
    Used for final minimization after all rotamer changes have been applied.
    
    Args:
        pdb_file: Path to PDB file
        verbose: Print details
        
    Returns:
        Dictionary with:
            - success: bool
            - total_fa_rep: float - total fa_rep energy for structure
            - total_score: float - total ref2015 score after minimization
            - minimized_pdb: str - PDB content of minimized structure
            - error: str (only if success=False)
    """
    try:
        # Load structure
        pose = pyrosetta.pose_from_pdb(pdb_file)
        
        # Create ref2015 score function
        sfxn = pyrosetta.create_score_function('ref2015')
        
        # Enable coordinate_constraint term to enforce CA position constraints
        sfxn.set_weight(ScoreType.coordinate_constraint, CA_CONSTRAINT_WEIGHT)
        
        # Add CA coordinate constraints to preserve backbone positions
        num_constraints = add_ca_coordinate_constraints(pose, exclude_resnum=None)
        
        # Create MoveMap: chi-only minimization (no backbone movement)
        # This preserves the global backbone structure while allowing all sidechains to relax
        movemap = MoveMap()
        movemap.set_bb(False)  # Disable backbone movement globally
        movemap.set_chi(True)  # Allow all sidechain chi angles to move
        
        # Create MinMover with ref2015
        # Reduced aggressiveness: fewer iterations, higher tolerance
        minmover = MinMover()
        minmover.movemap(movemap)
        minmover.score_function(sfxn)
        minmover.min_type('lbfgs_armijo_nonmonotone')
        minmover.tolerance(0.1)   # Increased from 0.01 - converge earlier
        minmover.max_iter(100)    # Reduced from 500 - fewer iterations
        
        # Score before minimization
        score_before = sfxn(pose)
        
        # Run minimization
        minmover.apply(pose)
        
        # Score after minimization
        total_score = sfxn(pose)
        
        # Sum fa_rep across all residues
        energies = pose.energies()
        total_fa_rep = 0.0
        for i in range(1, pose.total_residue() + 1):
            res_energies = energies.residue_total_energies(i)
            total_fa_rep += res_energies[ScoreType.fa_rep]
        
        # Write minimized structure to a temp file and read content
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False, mode='w') as f:
            temp_path = f.name
        
        pose.dump_pdb(temp_path)
        
        with open(temp_path, 'r') as f:
            minimized_pdb_raw = f.read()
        
        os.remove(temp_path)
        
        # Strip hydrogens to maintain consistency with input structure
        minimized_pdb = strip_hydrogens_from_pdb(minimized_pdb_raw)
        
        return {
            'success': True,
            'total_fa_rep': float(total_fa_rep),
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
        mode = input_data.get('mode', 'frozen')  # 'frozen' or 'full'
        chain = input_data.get('chain')
        resnum = input_data.get('resnum')
        verbose = input_data.get('verbose', False)
        
        if not pdb_file:
            result = {'success': False, 'error': 'pdb_file is required'}
        elif mode == 'full':
            # Full structure minimization (no frozen residues)
            result = minimize_full_structure(pdb_file, verbose)
        elif chain is None:
            result = {'success': False, 'error': 'chain is required for frozen mode'}
        elif resnum is None:
            result = {'success': False, 'error': 'resnum is required for frozen mode'}
        else:
            result = minimize_with_frozen_residue(pdb_file, chain, resnum, verbose)
        
        json.dump(result, sys.stdout)
        
    except Exception as e:
        json.dump({'success': False, 'error': str(e)}, sys.stdout)
        sys.exit(1)


if __name__ == '__main__':
    main()
