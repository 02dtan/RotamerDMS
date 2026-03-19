#!/usr/bin/env python3
"""
Standalone pyKVFinder runner script.
Runs in conda's Python (with pyKVFinder/numpy 2.x) and communicates via JSON.

Usage:
    echo '{"pdb_file": "path.pdb", "binding_site_residues": [...]}' | python pykvfinder_runner.py
"""

import sys
import json
import os

import pyKVFinder
import numpy as np


# ChimeraX default parameters for cavity detection
DEFAULT_PARAMS = {
    'step': 0.6,
    'probe_in': 1.4,
    'probe_out': 4.0,
    'removal_distance': 2.4,
    'volume_cutoff': 5.0,
}


def detect_cavities(pdb_file, params=None):
    """Run pyKVFinder cavity detection on a PDB file."""
    if params is None:
        params = DEFAULT_PARAMS
    
    atomic = pyKVFinder.read_pdb(pdb_file)
    step = params.get('step', DEFAULT_PARAMS['step'])
    
    # Get vertices first - required for detect()
    vertices = pyKVFinder.get_vertices(atomic, step=step)
    
    ncav, cavities = pyKVFinder.detect(
        atomic,
        vertices,
        step=step,
        probe_in=params.get('probe_in', DEFAULT_PARAMS['probe_in']),
        probe_out=params.get('probe_out', DEFAULT_PARAMS['probe_out']),
        removal_distance=params.get('removal_distance', DEFAULT_PARAMS['removal_distance']),
        volume_cutoff=params.get('volume_cutoff', DEFAULT_PARAMS['volume_cutoff']),
    )
    
    if ncav == 0:
        return {
            'num_cavities': 0,
            'volumes': {},
            'residues': {},
        }
    
    # pyKVFinder.spatial returns (surface, volume, area)
    surface, volumes, area = pyKVFinder.spatial(cavities, step=step)
    
    # pyKVFinder.constitutional returns residues dict
    residues = pyKVFinder.constitutional(cavities, atomic, vertices, step=step)
    
    return {
        'num_cavities': ncav,
        'volumes': volumes,
        'residues': residues,
    }


def parse_residue_list(res_list):
    """
    Parse pyKVFinder residue list: ['resnum', 'chain', 'resname']
    Returns (chain, resnum, resname) tuple to match our format.
    """
    if len(res_list) >= 3:
        try:
            resnum = int(res_list[0])
            chain = str(res_list[1])
            resname = str(res_list[2])
            return (chain, resnum, resname)
        except (ValueError, IndexError):
            pass
    return None


def measure_binding_site_volume(pdb_file, binding_site_residues, min_contacts=4, params=None):
    """Measure total volume of cavities contacting binding site residues."""
    bs_set = set()
    for res in binding_site_residues:
        if len(res) == 3:
            chain, resnum, resname = res
            bs_set.add((str(chain), int(resnum), str(resname)))
    
    results = detect_cavities(pdb_file, params)
    
    if results['num_cavities'] == 0:
        return {
            'total_volume': 0.0,
            'num_binding_cavities': 0,
            'cavity_details': [],
        }
    
    cavity_details = []
    total_volume = 0.0
    num_binding_cavities = 0
    
    for cav_name, cav_residues in results['residues'].items():
        contacts = 0
        contacting_residues = []
        
        # pyKVFinder.constitutional returns lists like ['resnum', 'chain', 'resname']
        for res_list in cav_residues:
            parsed = parse_residue_list(res_list)
            if parsed and parsed in bs_set:
                contacts += 1
                contacting_residues.append(list(parsed))
        
        cav_volume = results['volumes'].get(cav_name, 0.0)
        
        if contacts >= min_contacts:
            total_volume += cav_volume
            num_binding_cavities += 1
            cavity_details.append({
                'name': cav_name,
                'volume': cav_volume,
                'contacts': contacts,
                'contacting_residues': contacting_residues,
            })
    
    return {
        'total_volume': total_volume,
        'num_binding_cavities': num_binding_cavities,
        'cavity_details': cavity_details,
    }


def main():
    """Read JSON from stdin, write JSON to stdout."""
    try:
        input_data = json.load(sys.stdin)
        
        pdb_file = input_data.get('pdb_file')
        binding_site_residues = input_data.get('binding_site_residues', [])
        min_contacts = input_data.get('min_contacts', 4)
        params = input_data.get('params', None)
        
        if not pdb_file:
            result = {'error': 'pdb_file is required'}
        elif not os.path.exists(pdb_file):
            result = {'error': f'PDB file not found: {pdb_file}'}
        else:
            result = measure_binding_site_volume(pdb_file, binding_site_residues, min_contacts, params)
        
        json.dump(result, sys.stdout)
        
    except Exception as e:
        json.dump({'error': str(e)}, sys.stdout)
        sys.exit(1)


if __name__ == '__main__':
    main()
