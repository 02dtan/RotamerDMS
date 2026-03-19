"""
Preparation module for RotamerDMS.

Handles structure loading, sequence validation, and binding site residue identification.
"""

from schrodinger import structure
from schrodinger.structutils import analyze
import math


def load_structure(pdb_path):
    """
    Load a PDB structure using Schrodinger StructureReader.
    
    Args:
        pdb_path: Path to the PDB file
        
    Returns:
        Schrodinger Structure object
        
    Raises:
        IOError: If the file cannot be read
    """
    try:
        st = structure.StructureReader.read(pdb_path)
        return st
    except IOError as e:
        raise IOError(f"Error reading PDB file {pdb_path}: {e}")


def get_sequence(struct):
    """
    Extract the amino acid sequence from a structure.
    
    Args:
        struct: Schrodinger Structure object
        
    Returns:
        List of tuples: (chain, resnum, resname) for each residue
    """
    sequence = []
    for res in struct.residue:
        # Skip non-protein residues (ligands, water, etc.)
        if res.pdbres.strip() in ['HOH', 'EPE', 'WAT', 'NA', 'CL', 'MG', 'CA', 'ZN']:
            continue
        # Check if it's a standard amino acid or modified
        sequence.append((res.chain, res.resnum, res.pdbres.strip()))
    return sequence


def align_sequences(ref_struct, sample_struct):
    """
    Align sequences between reference and sample by matching chain/resnum.
    Trims residues from tails that are missing in the other structure.
    
    Args:
        ref_struct: Reference structure (Schrodinger Structure object)
        sample_struct: Sample structure (Schrodinger Structure object)
        
    Returns:
        Tuple of (common_residues, ref_only, sample_only) where:
        - common_residues: List of (chain, resnum, resname) present in both
        - ref_only: List of (chain, resnum, resname) only in reference
        - sample_only: List of (chain, resnum, resname) only in sample
    """
    ref_seq = get_sequence(ref_struct)
    sample_seq = get_sequence(sample_struct)
    
    # Build lookup sets by (chain, resnum)
    ref_lookup = {(r[0], r[1]): r[2] for r in ref_seq}
    sample_lookup = {(s[0], s[1]): s[2] for s in sample_seq}
    
    ref_keys = set(ref_lookup.keys())
    sample_keys = set(sample_lookup.keys())
    
    # Find common residues
    common_keys = ref_keys & sample_keys
    ref_only_keys = ref_keys - sample_keys
    sample_only_keys = sample_keys - ref_keys
    
    # Build result lists
    common_residues = []
    mismatches = []
    
    for key in sorted(common_keys):
        chain, resnum = key
        ref_resname = ref_lookup[key]
        sample_resname = sample_lookup[key]
        if ref_resname == sample_resname:
            common_residues.append((chain, resnum, ref_resname))
        else:
            mismatches.append((chain, resnum, ref_resname, sample_resname))
    
    ref_only = [(k[0], k[1], ref_lookup[k]) for k in sorted(ref_only_keys)]
    sample_only = [(k[0], k[1], sample_lookup[k]) for k in sorted(sample_only_keys)]
    
    # Print alignment summary
    print(f"\n=== Sequence Alignment Summary ===")
    print(f"Reference residues: {len(ref_seq)}")
    print(f"Sample residues: {len(sample_seq)}")
    print(f"Common residues: {len(common_residues)}")
    
    if ref_only:
        print(f"\nResidues in REFERENCE only (excluded from analysis):")
        for chain, resnum, resname in ref_only:
            print(f"  {chain}:{resname}{resnum}")
    
    if sample_only:
        print(f"\nResidues in SAMPLE only (excluded from analysis):")
        for chain, resnum, resname in sample_only:
            print(f"  {chain}:{resname}{resnum}")
    
    if mismatches:
        print(f"\nWARNING: Residue type mismatches at same positions:")
        for chain, resnum, ref_res, sample_res in mismatches:
            print(f"  {chain}:{resnum} - ref={ref_res}, sample={sample_res}")
    
    if not common_residues:
        raise ValueError("No common residues found between reference and sample structures.")
    
    print(f"\nProceeding with {len(common_residues)} common residues.")
    return common_residues, ref_only, sample_only


def get_ligand_atoms(struct, ligand_code):
    """
    Select all atoms belonging to the specified ligand.
    
    Args:
        struct: Schrodinger Structure object
        ligand_code: 3-letter residue code for the ligand
        
    Returns:
        List of atom objects belonging to the ligand
        
    Raises:
        ValueError: If ligand is not found
    """
    ligand_atoms = []
    for atom in struct.atom:
        if atom.pdbres.strip() == ligand_code:
            ligand_atoms.append(atom)
    
    if not ligand_atoms:
        raise ValueError(f"Ligand '{ligand_code}' not found in structure.")
    
    print(f"Found ligand '{ligand_code}' with {len(ligand_atoms)} atoms.")
    return ligand_atoms


def calculate_distance(atom1, atom2):
    """
    Calculate Euclidean distance between two atoms.
    
    Args:
        atom1, atom2: Schrodinger Atom objects
        
    Returns:
        Distance in Angstroms
    """
    dx = atom1.x - atom2.x
    dy = atom1.y - atom2.y
    dz = atom1.z - atom2.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def get_binding_site_residues(struct, ligand_code, distance_threshold=7.0):
    """
    Find all residues within the specified distance from the ligand.
    
    Args:
        struct: Schrodinger Structure object (reference structure with ligand)
        ligand_code: 3-letter residue code for the ligand
        distance_threshold: Distance in Angstroms (default 7.0)
        
    Returns:
        List of tuples: (chain, resnum, resname) for residues within threshold
    """
    ligand_atoms = get_ligand_atoms(struct, ligand_code)
    
    binding_site_residues = set()
    
    for res in struct.residue:
        # Skip the ligand itself
        if res.pdbres.strip() == ligand_code:
            continue
        # Skip non-protein residues
        if res.pdbres.strip() in ['HOH', 'EPE', 'WAT', 'NA', 'CL', 'MG', 'CA', 'ZN']:
            continue
        
        # Check if any atom in this residue is within threshold of any ligand atom
        for res_atom in res.atom:
            min_dist = float('inf')
            for lig_atom in ligand_atoms:
                dist = calculate_distance(res_atom, lig_atom)
                if dist < min_dist:
                    min_dist = dist
            
            if min_dist <= distance_threshold:
                binding_site_residues.add((res.chain, res.resnum, res.pdbres.strip()))
                break  # Found one atom within threshold, move to next residue
    
    # Convert to sorted list
    residue_list = sorted(list(binding_site_residues), key=lambda x: (x[0], x[1]))
    
    print(f"\n=== Binding Site Residues (within {distance_threshold} Å of ligand) ===")
    print(f"Total: {len(residue_list)} residues")
    for chain, resnum, resname in residue_list:
        print(f"  {chain}:{resname}{resnum}")
    
    return residue_list


def get_residue_by_id(struct, chain, resnum):
    """
    Get a residue object by chain and residue number.
    
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


def map_residues_to_sample(binding_site_residues, sample_struct):
    """
    Map binding site residues from reference to sample structure.
    Since we've validated sequence identity, we map by residue number.
    
    Args:
        binding_site_residues: List of (chain, resnum, resname) from reference
        sample_struct: Sample structure
        
    Returns:
        List of (chain, resnum, resname) tuples validated in sample structure
    """
    mapped_residues = []
    
    for chain, resnum, resname in binding_site_residues:
        res = get_residue_by_id(sample_struct, chain, resnum)
        if res is None:
            print(f"Warning: Residue {chain}:{resname}{resnum} not found in sample structure.")
            continue
        if res.pdbres.strip() != resname:
            print(f"Warning: Residue mismatch at {chain}:{resnum} - "
                  f"expected {resname}, found {res.pdbres.strip()}")
            continue
        mapped_residues.append((chain, resnum, resname))
    
    return mapped_residues
