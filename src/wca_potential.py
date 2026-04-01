"""
Energy minimization and fa_rep extraction module.

Performs local energy minimization on structures after rotamer changes,
with the target rotamer residue frozen. This allows surrounding residues
to relax and relieve steric clashes while preserving the desired rotamer state.

After minimization, extracts the fa_rep (repulsive Lennard-Jones) energy
for the target residue as a measure of remaining steric strain.

PyRosetta is called via subprocess to avoid conflicts between Schrodinger's
Python environment and the conda environment where PyRosetta is installed.
"""

import os
import subprocess
import json
import tempfile
import sys

# Path to the PyRosetta runner script and conda Python
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PYROSETTA_RUNNER = os.path.join(SCRIPT_DIR, "pyrosetta_runner.py")

# Conda environment configuration
CONDA_ENV_NAME = "rotamer_dms"


def _find_conda_python():
    """Find the Python executable in the rotamer_dms conda environment."""
    conda_bases = []
    
    # Check CONDA_EXE environment variable
    conda_exe = os.environ.get('CONDA_EXE')
    if conda_exe:
        conda_bases.append(os.path.dirname(os.path.dirname(conda_exe)))
    
    # Check common paths
    home = os.path.expanduser('~')
    conda_bases.extend([
        os.path.join(home, 'miniconda3'),
        os.path.join(home, 'anaconda3'),
        os.path.join(home, 'miniforge3'),
        '/opt/conda',
    ])
    
    # Also check the software directory specific to this environment
    conda_bases.append('/lustre/fs6/lyu_lab/scratch/dantan/software/miniconda3')
    
    for base in conda_bases:
        python_path = os.path.join(base, 'envs', CONDA_ENV_NAME, 'bin', 'python')
        if os.path.exists(python_path):
            return python_path
    
    return None


CONDA_PYTHON = _find_conda_python()
_pyrosetta_available = None


def _check_pyrosetta_available():
    """Check if PyRosetta is available via subprocess."""
    global _pyrosetta_available
    
    if _pyrosetta_available is not None:
        return _pyrosetta_available
    
    if CONDA_PYTHON is None:
        print(f"Warning: Could not find conda environment '{CONDA_ENV_NAME}'")
        _pyrosetta_available = False
        return False
    
    if not os.path.exists(PYROSETTA_RUNNER):
        print(f"Warning: PyRosetta runner script not found: {PYROSETTA_RUNNER}")
        _pyrosetta_available = False
        return False
    
    _pyrosetta_available = True
    return True


def _save_structure_temp(schrodinger_struct):
    """Save Schrodinger structure to a temporary PDB file."""
    from schrodinger import structure as schrod_structure
    
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False, mode='w') as f:
        temp_path = f.name
    
    with schrod_structure.StructureWriter(temp_path) as writer:
        writer.append(schrodinger_struct)
    
    return temp_path


def _parse_pdb_coordinates(pdb_content):
    """
    Parse atom coordinates from PDB content.
    
    Returns a dictionary mapping (chain, resnum, atom_name) -> (x, y, z)
    Only parses ATOM/HETATM records for heavy atoms (hydrogens already stripped).
    """
    coords = {}
    for line in pdb_content.split('\n'):
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        
        # PDB format columns (0-indexed):
        # 12-15: atom name
        # 17-19: residue name
        # 21: chain ID
        # 22-25: residue sequence number
        # 30-37: x coordinate
        # 38-45: y coordinate
        # 46-53: z coordinate
        try:
            atom_name = line[12:16].strip()
            chain = line[21:22].strip()
            resnum = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            
            # Use (chain, resnum, atom_name) as key
            # If chain is empty, use ' ' as placeholder
            if not chain:
                chain = ' '
            coords[(chain, resnum, atom_name)] = (x, y, z)
        except (ValueError, IndexError):
            continue
    
    return coords


def _transfer_coordinates(original_struct, minimized_pdb_content, verbose=False):
    """
    Transfer coordinates from minimized PDB back to original Schrodinger structure.
    
    This preserves the original structure's connectivity, atom properties, and
    residue identifiers while updating only the XYZ coordinates.
    
    Args:
        original_struct: Original Schrodinger Structure object
        minimized_pdb_content: PDB content string from PyRosetta minimization
        verbose: Print debug information
        
    Returns:
        Dictionary with:
            - success: bool
            - structure: Updated Schrodinger Structure (copy of original with new coords)
            - atoms_updated: Number of atoms that had coordinates updated
            - atoms_total: Total number of heavy atoms in original structure
            - atoms_missing: List of atoms not found in minimized PDB
            - rmsd: RMSD between original and updated coordinates
            - error: str (if failed)
    """
    import math
    
    # Parse coordinates from minimized PDB
    minimized_coords = _parse_pdb_coordinates(minimized_pdb_content)
    
    if not minimized_coords:
        return {
            'success': False,
            'error': 'No coordinates found in minimized PDB'
        }
    
    # Make a copy of the original structure to modify
    updated_struct = original_struct.copy()
    
    atoms_updated = 0
    atoms_missing = []
    sum_sq_dist = 0.0
    
    for atom in updated_struct.atom:
        # Skip hydrogens (we only care about heavy atoms)
        if atom.atomic_number == 1:
            continue
        
        # Get atom identifiers
        chain = atom.chain.strip() if atom.chain else ' '
        resnum = atom.resnum
        atom_name = atom.pdbname.strip()
        
        # Look up coordinates in minimized structure
        key = (chain, resnum, atom_name)
        
        if key in minimized_coords:
            new_x, new_y, new_z = minimized_coords[key]
            
            # Calculate squared distance for RMSD
            dx = new_x - atom.x
            dy = new_y - atom.y
            dz = new_z - atom.z
            sum_sq_dist += dx*dx + dy*dy + dz*dz
            
            # Update coordinates
            atom.x = new_x
            atom.y = new_y
            atom.z = new_z
            atoms_updated += 1
        else:
            atoms_missing.append(f"{chain}:{resnum}:{atom_name}")
    
    # Count total heavy atoms
    atoms_total = sum(1 for a in original_struct.atom if a.atomic_number != 1)
    
    # Calculate RMSD
    rmsd = math.sqrt(sum_sq_dist / atoms_updated) if atoms_updated > 0 else 0.0
    
    if verbose:
        print(f"    Coordinate transfer: {atoms_updated}/{atoms_total} atoms updated, RMSD={rmsd:.3f} Å")
        if atoms_missing:
            print(f"    Warning: {len(atoms_missing)} atoms not found in minimized structure")
    
    # Validation: check that we updated most atoms
    update_ratio = atoms_updated / atoms_total if atoms_total > 0 else 0
    if update_ratio < 0.9:
        return {
            'success': False,
            'error': f'Only {atoms_updated}/{atoms_total} atoms ({update_ratio:.1%}) could be mapped. '
                     f'Missing atoms: {atoms_missing[:10]}...'
        }
    
    return {
        'success': True,
        'structure': updated_struct,
        'atoms_updated': atoms_updated,
        'atoms_total': atoms_total,
        'atoms_missing': atoms_missing,
        'rmsd': rmsd,
    }


def _load_structure_from_pdb_content(pdb_content):
    """
    DEPRECATED: Load a Schrodinger Structure from PDB content string.
    
    This function is kept for backwards compatibility but should not be used.
    Use _transfer_coordinates() instead to preserve structure consistency.
    """
    from schrodinger import structure as schrod_structure
    
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False, mode='w') as f:
        f.write(pdb_content)
        temp_path = f.name
    
    try:
        structs = list(schrod_structure.StructureReader(temp_path))
        if structs:
            return structs[0]
        return None
    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


def run_minimization_subprocess(pdb_file, chain, resnum, verbose=False):
    """
    Run PyRosetta minimization in a subprocess using conda's Python.
    
    The target residue is frozen during minimization to preserve the rotamer state.
    
    Args:
        pdb_file: Path to PDB file
        chain: Chain identifier
        resnum: Residue number (to freeze)
        verbose: Print details
        
    Returns:
        Dictionary with:
            - success: bool
            - fa_rep_energy: float
            - minimized_pdb: str (PDB content)
            - total_score: float
            - error: str (if failed)
    """
    if not _check_pyrosetta_available():
        return {'success': False, 'error': 'PyRosetta not available'}
    
    input_data = {
        'pdb_file': pdb_file,
        'chain': chain,
        'resnum': resnum,
        'verbose': verbose,
    }
    
    # Run the subprocess with clean environment (remove Schrodinger's PYTHONPATH)
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    
    try:
        result = subprocess.run(
            [CONDA_PYTHON, PYROSETTA_RUNNER],
            input=json.dumps(input_data),
            capture_output=True,
            text=True,
            timeout=300,  # Increased timeout for minimization
            env=env,
        )
        
        if result.returncode != 0:
            error_msg = result.stderr.strip() if result.stderr else ""
            stdout_msg = result.stdout.strip() if result.stdout else ""
            full_error = f"returncode={result.returncode}"
            if error_msg:
                full_error += f", stderr={error_msg}"
            if stdout_msg and len(stdout_msg) < 500:
                full_error += f", stdout={stdout_msg}"
            return {'success': False, 'error': f'PyRosetta subprocess failed: {full_error}'}
        
        output = json.loads(result.stdout)
        return output
        
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'PyRosetta minimization timed out (>300s)'}
    except json.JSONDecodeError as e:
        return {'success': False, 'error': f'Failed to parse PyRosetta output: {e}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}


def minimize_and_get_fa_rep(schrodinger_struct, chain, resnum, verbose=False):
    """
    Minimize structure with frozen target residue and return fa_rep energy.
    
    This function:
    1. Saves the Schrodinger structure to a temp PDB
    2. Calls PyRosetta to minimize with ref2015, freezing the target residue
    3. Returns the minimized structure and fa_rep energy for the target
    
    The minimization allows surrounding residues to relax, relieving steric
    clashes while preserving the target rotamer conformation.
    
    Args:
        schrodinger_struct: Schrodinger Structure object with rotamer applied
        chain: Chain identifier for target residue
        resnum: Residue number of target (will be frozen)
        verbose: Print progress
        
    Returns:
        Dictionary with:
            - success: bool
            - fa_rep_energy: float - repulsive energy for target residue
            - minimized_struct: Schrodinger Structure - minimized structure
            - total_score: float - total ref2015 score
            - error: str (if failed)
    """
    # Save structure to temp PDB
    temp_pdb = _save_structure_temp(schrodinger_struct)
    
    try:
        result = run_minimization_subprocess(temp_pdb, chain, resnum, verbose)
        
        if not result.get('success'):
            print(f"Minimization failed for {chain}:{resnum}: {result.get('error', 'unknown')}", 
                  file=sys.stderr)
            return result
        
        # Transfer coordinates from minimized PDB back to original structure
        # This preserves connectivity, atom properties, and residue identifiers
        transfer_result = _transfer_coordinates(
            schrodinger_struct, 
            result['minimized_pdb'],
            verbose=verbose
        )
        
        if not transfer_result.get('success'):
            return {
                'success': False,
                'error': f"Coordinate transfer failed: {transfer_result.get('error', 'unknown')}"
            }
        
        minimized_struct = transfer_result['structure']
        
        if verbose:
            print(f"    Minimized {chain}:{resnum}: fa_rep={result['fa_rep_energy']:.2f} REU, "
                  f"fa_dun={result.get('fa_dun_energy', 0):.2f} REU "
                  f"(score: {result.get('score_before', 0):.1f} → {result['total_score']:.1f})")
        
        return {
            'success': True,
            'fa_rep_energy': result['fa_rep_energy'],
            'fa_dun_energy': result.get('fa_dun_energy'),
            'minimized_struct': minimized_struct,
            'total_score': result['total_score'],
            'score_before': result.get('score_before'),
            'coord_transfer_rmsd': transfer_result['rmsd'],
            'atoms_updated': transfer_result['atoms_updated'],
            'atoms_total': transfer_result['atoms_total'],
        }
        
    finally:
        # Clean up temp file
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)


def get_fa_rep_energy(schrodinger_struct, chain, resnum, verbose=False):
    """
    Simple interface to get fa_rep energy after minimization.
    
    Note: This discards the minimized structure. Use minimize_and_get_fa_rep()
    if you need both the energy and the minimized structure.
    
    Args:
        schrodinger_struct: Schrodinger Structure with rotamer state applied
        chain: Chain identifier
        resnum: Residue number
        verbose: Print details
        
    Returns:
        fa_rep energy (float) or None if computation failed
    """
    result = minimize_and_get_fa_rep(schrodinger_struct, chain, resnum, verbose=verbose)
    
    if result.get('success'):
        return result['fa_rep_energy']
    else:
        return None


def run_full_minimization_subprocess(pdb_file, verbose=False):
    """
    Run full structure minimization in PyRosetta subprocess (no frozen residues).
    
    Args:
        pdb_file: Path to PDB file
        verbose: Print details
        
    Returns:
        Dictionary with:
            - success: bool
            - total_fa_rep: float - total fa_rep energy
            - total_score: float - total ref2015 score
            - minimized_pdb: str - PDB content of minimized structure
            - error: str (if failed)
    """
    if not _check_pyrosetta_available():
        return {'success': False, 'error': 'PyRosetta not available'}
    
    input_data = {
        'pdb_file': pdb_file,
        'mode': 'full',
        'verbose': verbose,
    }
    
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    
    try:
        result = subprocess.run(
            [CONDA_PYTHON, PYROSETTA_RUNNER],
            input=json.dumps(input_data),
            capture_output=True,
            text=True,
            timeout=600,  # 10 minute timeout for full structure
            env=env,
        )
        
        if result.returncode != 0:
            error_msg = result.stderr or result.stdout or "Unknown error"
            return {'success': False, 'error': f'Subprocess failed: {error_msg}'}
        
        output = json.loads(result.stdout)
        return output
        
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'PyRosetta minimization timed out (>600s)'}
    except json.JSONDecodeError as e:
        return {'success': False, 'error': f'Failed to parse PyRosetta output: {e}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}


def minimize_full_structure(schrodinger_struct, verbose=False):
    """
    Minimize entire structure with PyRosetta (no frozen residues).
    
    Used for final minimization after all rotamer changes have been applied.
    
    Args:
        schrodinger_struct: Schrodinger Structure object
        verbose: Print progress
        
    Returns:
        Dictionary with:
            - success: bool
            - total_fa_rep: float - total fa_rep energy
            - minimized_struct: Schrodinger Structure - minimized structure
            - total_score: float - total ref2015 score
            - error: str (if failed)
    """
    # Save structure to temp PDB
    temp_pdb = _save_structure_temp(schrodinger_struct)
    
    try:
        result = run_full_minimization_subprocess(temp_pdb, verbose)
        
        if not result.get('success'):
            print(f"Full minimization failed: {result.get('error', 'unknown')}", 
                  file=sys.stderr)
            return result
        
        # Transfer coordinates from minimized PDB back to original structure
        # This preserves connectivity, atom properties, and residue identifiers
        transfer_result = _transfer_coordinates(
            schrodinger_struct, 
            result['minimized_pdb'],
            verbose=verbose
        )
        
        if not transfer_result.get('success'):
            return {
                'success': False,
                'error': f"Coordinate transfer failed: {transfer_result.get('error', 'unknown')}"
            }
        
        minimized_struct = transfer_result['structure']
        
        if verbose:
            print(f"    Full minimization complete: fa_rep={result['total_fa_rep']:.2f} REU "
                  f"(score: {result.get('score_before', 0):.1f} → {result['total_score']:.1f})")
            print(f"    Coordinate transfer: RMSD={transfer_result['rmsd']:.3f} Å, "
                  f"{transfer_result['atoms_updated']}/{transfer_result['atoms_total']} atoms")
        
        return {
            'success': True,
            'total_fa_rep': result['total_fa_rep'],
            'minimized_struct': minimized_struct,
            'total_score': result['total_score'],
            'score_before': result.get('score_before'),
            'coord_transfer_rmsd': transfer_result['rmsd'],
            'atoms_updated': transfer_result['atoms_updated'],
            'atoms_total': transfer_result['atoms_total'],
        }
        
    finally:
        # Clean up temp file
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)
