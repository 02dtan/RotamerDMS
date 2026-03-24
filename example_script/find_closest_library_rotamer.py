import numpy as np
from schrodinger import structure
from schrodinger.protein import rotamers

def get_current_chi_angles(st, residue):
    """Measures the current chi angles of a residue in the structure."""
    chi_angles = []
    # Schrödinger supports up to 5 chi angles depending on the residue type
    for i in range(1, 6):
        try:
            atom_indices = rotamers.get_chi_atom_indices(st, residue, i)
            if atom_indices:
                angle = st.measure(*atom_indices)
                chi_angles.append(angle)
        except (ValueError, AttributeError):
            break
    return chi_angles

def find_closest_rotamer(st_file):
    st = structure.StructureReader.read(st_file)
    
    for residue in st.residue:
        if residue.pdbres.strip() in ['GLY', 'ALA']:
            continue
            
        ca_atom = residue.getAlphaCarbon()
        if not ca_atom:
            continue
            
        try:
            res_rotamers = rotamers.Rotamers(st, ca_atom)
            current_chi = np.array(get_current_chi_angles(st, residue))
            
            best_idx = -1
            min_dist = float('inf')
            best_prob = 0.0
            
            # Iterate through the library rotamers
            for i, rot in enumerate(res_rotamers.rotamers):
                # The 'chi' attribute is a list of angles for that library state
                library_chi = np.array(rot.chi)
                
                # Use periodicity-aware distance if necessary, 
                # but simple Euclidean distance works for basic identification
                dist = np.linalg.norm(current_chi - library_chi)
                
                if dist < min_dist:
                    min_dist = dist
                    best_idx = i
                    best_prob = rot.probability  # Probability from the library
            
            print(f"Residue {residue.pdbres}{residue.resnum}:")
            print(f"  Closest Rotamer Index: {best_idx}")
            print(f"  Library Probability:   {best_prob:.4f}")
            print(f"  Angular Distance:      {min_dist:.2f} degrees")
            
        except Exception as e:
            print(f"Error processing {residue}: {e}")

# Usage
find_closest_rotamer('../data/s0067_rec.pdb')