from schrodinger import structure
from schrodinger.structutils import build
from schrodinger.protein import rotamers

input_file = "/lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS/data/s0067_rec.pdb"
output_file = "/lustre/fs6/lyu_lab/scratch/dantan/RotamerDMS/output/s0067_rec_optimized.pdb"
chain_id = 'A'
res_num = 102
new_res_code = 'GLY'

try:
    st = structure.StructureReader.read(input_file)
except IOError as e:
    print(f"Error reading PDB file: {e}")
    exit()

target_res = None
for res in st.residue:
    if res.chain == chain_id and res.resnum == res_num:
        target_res = res
        break
    
if target_res is None:
    print(f"Target residue {chain_id}:{res_num} not found")
    exit()

build.mutate(st, target_res.atom[1], new_res_code)
print(f"Mutated residue {chain_id}:{res_num} to {new_res_code}")

with structure.StructureWriter(output_file) as writer:
    writer.append(st)
print(f"Mutated structure saved to {output_file}")