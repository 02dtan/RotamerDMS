#!/bin/bash
#
# RotamerDMS Runner Script
# 
# This script sets up the required environment and runs the RotamerDMS algorithm.
#
# Usage:
#   ./run_rotamer_dms.sh --reference <ref.pdb> --sample <sample.pdb> --ligand <code> [options]
#
# Options:
#   --reference, -r    Path to experimental reference structure with bound ligand
#   --sample, -s       Path to sample structure to optimize
#   --ligand, -l       3-letter residue code of the ligand
#   --distance, -d     Distance threshold in Angstroms (default: 7.0)
#   --top-percent, -t  Percentage of top residues to sample (default: 100)
#   --min-contacts, -m Minimum binding site contacts for cavity (default: 4)
#   --output-dir, -o   Output directory (default: ./output)
#   --checkpoint-dir   Checkpoint directory (default: ./checkpoints)
#   --quiet, -q        Reduce output verbosity
#   --help, -h         Show this help message
#

set -e

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Show help if no arguments
if [ $# -eq 0 ]; then
    echo "RotamerDMS - Binding Pocket Volume Maximization"
    echo ""
    echo "Usage: $0 --reference <ref.pdb> --sample <sample.pdb> --ligand <code> [options]"
    echo ""
    echo "Required arguments:"
    echo "  --reference, -r    Path to experimental reference structure with bound ligand"
    echo "  --sample, -s       Path to sample structure to optimize"
    echo "  --ligand, -l       3-letter residue code of the ligand in reference structure"
    echo ""
    echo "Optional arguments:"
    echo "  --distance, -d     Distance threshold (Å) for binding site selection (default: 7.0)"
    echo "  --top-percent, -t  Percentage of top volume-affecting residues to sample (default: 100)"
    echo "  --min-contacts, -m Minimum binding site residues a cavity must contact (default: 4)"
    echo "  --output-dir, -o   Output directory (default: ./output)"
    echo "  --checkpoint-dir   Checkpoint directory (default: ./checkpoints)"
    echo "  --quiet, -q        Reduce output verbosity"
    echo "  --help, -h         Show this help message"
    echo ""
    echo "Example:"
    echo "  $0 --reference data/4W59_n-hexylbenzene.pdb --sample data/s0067_rec.pdb --ligand 3GZ"
    echo ""
    exit 0
fi

# Check for help flag
for arg in "$@"; do
    if [ "$arg" == "--help" ] || [ "$arg" == "-h" ]; then
        echo "RotamerDMS - Binding Pocket Volume Maximization"
        echo ""
        echo "Usage: $0 --reference <ref.pdb> --sample <sample.pdb> --ligand <code> [options]"
        echo ""
        echo "Required arguments:"
        echo "  --reference, -r    Path to experimental reference structure with bound ligand"
        echo "  --sample, -s       Path to sample structure to optimize"
        echo "  --ligand, -l       3-letter residue code of the ligand in reference structure"
        echo ""
        echo "Optional arguments:"
        echo "  --distance, -d     Distance threshold (Å) for binding site selection (default: 7.0)"
        echo "  --top-percent, -t  Percentage of top volume-affecting residues to sample (default: 100)"
        echo "  --output-dir, -o   Output directory (default: ./output)"
        echo "  --checkpoint-dir   Checkpoint directory (default: ./checkpoints)"
        echo "  --quiet, -q        Reduce output verbosity"
        echo "  --help, -h         Show this help message"
        echo ""
        echo "Example:"
        echo "  $0 --reference data/4W59_n-hexylbenzene.pdb --sample data/s0067_rec.pdb --ligand 3GZ"
        echo ""
        exit 0
    fi
done

echo "=============================================="
echo "RotamerDMS Environment Setup"
echo "=============================================="

# Set up Schrodinger environment
echo "Setting up Schrodinger environment..."
export SB_BASE_OVERRIDE=/lustre/fs4/ruit/scratch/rebecca/LYU_SBGRID/INST
# Temporarily disable exit-on-error since sbgrid.shrc may return non-zero
set +e
source /lustre/fs4/ruit/scratch/rebecca/LYU_SBGRID/INST/sbgrid.shrc
set -e

echo "Environment setup complete."
echo ""

# Change to script directory
cd "$SCRIPT_DIR"

# Run the main Python script with Schrodinger's Python
echo "Running RotamerDMS..."
echo ""

$SCHRODINGER/run python3 run_rotamer_dms.py "$@"
