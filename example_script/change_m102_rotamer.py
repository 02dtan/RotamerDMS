from schrodinger import structure
from schrodinger.protein import rotamers

input_file = "/lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS/data/s0067_rec.pdb"
output_file = "/lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS/output/s0067_rec_optimized.pdb"

# Read
try:
    st = structure.StructureReader.read(input_file)
except IOError as e:
    print(f"Error reading PDB file: {e}")
    exit()

# Handle Met102
target_res = None
for res in st.residue:
    if res.chain == 'A' and res.resnum == 102:
        target_res = res
        break

if target_res:
    rotamer_lib = rotamers.Rotamers(st, target_res.atom[1])
    all_rotamers = rotamer_lib.rotamers

    if all_rotamers:
        desired_rotamer = all_rotamers[1]
        desired_rotamer.apply()
        print("Rotamer applied to Met102.")

# Handle Leu118
target_res = None
for res in st.residue:
    if res.chain == 'A' and res.resnum == 118:
        target_res = res
        break

if target_res:
    rotamer_lib = rotamers.Rotamers(st, target_res.atom[1])
    all_rotamers = rotamer_lib.rotamers

    if all_rotamers:
        desired_rotamer = all_rotamers[2]
        desired_rotamer.apply()
        print("Rotamer applied to Leu118.")

# Handle Val87
target_res = None
for res in st.residue:
    if res.chain == 'A' and res.resnum == 87:
        target_res = res
        break

if target_res:
    rotamer_lib = rotamers.Rotamers(st, target_res.atom[1])
    all_rotamers = rotamer_lib.rotamers

    if all_rotamers:
        desired_rotamer = all_rotamers[2]
        desired_rotamer.apply()
        print("Rotamer applied to Val87.")

with structure.StructureWriter(output_file) as writer:
    writer.append(st)
print(f"Mutated structure saved to {output_file}")