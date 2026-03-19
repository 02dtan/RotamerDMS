"""
Cavity detection module for RotamerDMS.

Uses pyKVFinder for cavity detection and volume measurement.
Parameters are set to match ChimeraX "Find Cavities" defaults.

pyKVFinder is called via subprocess to avoid numpy version conflicts
between Schrodinger's Python and pyKVFinder's requirements.
"""

import subprocess
import json
import tempfile
import os


# ChimeraX Find Cavities default parameters
DEFAULT_PARAMS = {
    'step': 0.6,              # Grid spacing (Å)
    'probe_in': 1.4,          # Inner probe radius (Å)
    'probe_out': 4.0,         # Outer probe radius (Å)
    'removal_distance': 2.4,  # Exterior trim distance (Å)
    'volume_cutoff': 5.0,     # Minimum cavity volume (Å³)
}

# Path to conda Python with pyKVFinder installed
CONDA_PYTHON = "/lustre/fs6/lyu_lab/scratch/dantan/software/miniconda3/envs/rotamer_dms/bin/python"

# Path to the pyKVFinder runner script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PYKVFINDER_RUNNER = os.path.join(SCRIPT_DIR, "pykvfinder_runner.py")

def save_structure_temp(struct):
    """
    Save a Schrodinger Structure to a temporary PDB file.
    
    Args:
        struct: Schrodinger Structure object
        
    Returns:
        Path to the temporary PDB file
    """
    from schrodinger import structure
    
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
    temp_path = temp_file.name
    temp_file.close()
    
    with structure.StructureWriter(temp_path) as writer:
        writer.append(struct)
    
    return temp_path


def run_pykvfinder_subprocess(pdb_file, binding_site_residues, min_contacts=4, params=None):
    """
    Run pyKVFinder in a subprocess using conda's Python.
    
    This avoids numpy version conflicts between Schrodinger and pyKVFinder.
    
    Args:
        pdb_file: Path to PDB file
        binding_site_residues: List of (chain, resnum, resname) tuples
        min_contacts: Minimum binding site contacts for a cavity
        params: Optional pyKVFinder parameters
        
    Returns:
        Dictionary with cavity detection results
    """
    # Prepare input JSON
    input_data = {
        'pdb_file': pdb_file,
        'binding_site_residues': [list(r) for r in binding_site_residues],
        'min_contacts': min_contacts,
        'params': params,
    }
    
    # Run the subprocess with clean environment (remove Schrodinger's PYTHONPATH)
    env = os.environ.copy()
    # Remove any PYTHONPATH that might interfere with conda's Python
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    
    try:
        result = subprocess.run(
            [CONDA_PYTHON, PYKVFINDER_RUNNER],
            input=json.dumps(input_data),
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
            env=env,
        )
        
        if result.returncode != 0:
            error_msg = result.stderr or result.stdout or "Unknown error"
            raise RuntimeError(f"pyKVFinder subprocess failed (exit {result.returncode}): {error_msg}")
        
        output = json.loads(result.stdout)
        
        if 'error' in output:
            raise RuntimeError(f"pyKVFinder error: {output['error']}")
        
        return output
        
    except subprocess.TimeoutExpired:
        raise RuntimeError("pyKVFinder subprocess timed out")
    except json.JSONDecodeError as e:
        raise RuntimeError(f"Failed to parse pyKVFinder output: {e}")


def measure_binding_site_volume(struct, binding_site_residues, min_contacts=4, params=None, verbose=True):
    """
    Measure the total volume of cavities contacting the binding site.
    
    This is the main function for volume measurement in the optimization loop.
    Uses subprocess to call pyKVFinder in conda's Python environment.
    
    Args:
        struct: Schrodinger Structure object
        binding_site_residues: List of (chain, resnum, resname) tuples
        min_contacts: Minimum number of binding site residues a cavity must contact (default 4)
        params: Optional pyKVFinder parameters
        verbose: Whether to print detailed information
        
    Returns:
        Dictionary containing:
            - 'total_volume': sum of binding site cavity volumes (Å³)
            - 'num_cavities': number of cavities contacting binding site
            - 'cavity_details': list with per-cavity information
            - 'success': whether any suitable cavities were found
    """
    # Save structure to temp file
    temp_pdb = save_structure_temp(struct)
    
    try:
        # Run pyKVFinder via subprocess
        result = run_pykvfinder_subprocess(
            temp_pdb,
            binding_site_residues,
            min_contacts=min_contacts,
            params=params
        )
        
        total_volume = result.get('total_volume', 0.0)
        num_cavities = result.get('num_binding_cavities', 0)
        cavity_details = result.get('cavity_details', [])
        
        if num_cavities == 0:
            if verbose:
                print(f"No cavities found contacting at least {min_contacts} binding site residues.")
            return {
                'total_volume': 0.0,
                'num_cavities': 0,
                'cavity_details': [],
                'success': False
            }
        
        if verbose:
            print(f"\n=== Binding Site Cavity Analysis ===")
            print(f"Found {num_cavities} cavity(ies) contacting binding site")
            print(f"Total binding site volume: {total_volume:.2f} Å^3")
            
            for cav in cavity_details:
                print(f"\n  {cav['name']}: {cav['volume']:.2f} Å^3")
                print(f"    Binding site contacts: {cav['contacts']}")
        
        return {
            'total_volume': total_volume,
            'num_cavities': num_cavities,
            'cavity_details': cavity_details,
            'success': True
        }
        
    finally:
        if os.path.exists(temp_pdb):
            os.remove(temp_pdb)


