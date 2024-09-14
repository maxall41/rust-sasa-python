import rust_sasa_python
import time
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.MMCIFParser import MMCIFParser

start = time.time()
v = rust_sasa_python.calculate_sasa_at_residue_level("example.cif")
print(v)
end = time.time()
print("Time for RustSasa:")
print(end - start)

start = time.time()

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("TEST", "example.cif")
sr = ShrakeRupley()
sr.compute(structure[0], level="R")

end = time.time()

print("Time for Biopython:")
print(end - start)