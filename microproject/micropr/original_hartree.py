from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from pyscf import gto, scf, tools
from pyscf.tools import cubegen
#import matplotlib.pyplot as plt
from IPython.display import display

def formula_to_3d_geometry(formula):    
    smiles = formula 

    mol = Chem.MolFromSmiles(smiles)
    for bond in mol.GetBonds():
        print(f"Связь: {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()}, Тип: {bond.GetBondType()}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)  
    AllChem.MMFFOptimizeMolecule(mol)  
    #tools.cubegen.orbital(mol, 'orbital.cube', mo_coeff, resolution=0.2)
    display(Draw.MolToImage(mol, size=(400, 400)))
    atoms = []
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        atoms.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
    
    return "\n".join(atoms)

def calculate_electronic_structure(geometry):
    #PySCF.
    mol = gto.M(
        atom=geometry,
        basis="6-31g*",
        spin=0,
        verbose=4
    )
    mf = scf.RHF(mol).run()
    return mf

"""def calculate_electronic_structure(geometry):
    mol = gto.M(
        atom=geometry,
        basis="6-31g*",
        spin=0,
        verbose=4
    )
    mf = my.MyRHF(mol)
    mf.kernel()
    return mf"""

if __name__ == "__main__":
    formula = "CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O"
    
    try:

        geometry = formula_to_3d_geometry(formula)
        print("3D-геометрия:\n", geometry)
        

        mf = calculate_electronic_structure(geometry)
        print(f"Энергия Хартри-Фока: {mf.e_tot:.6f} Hartree")
        
        #tools.molden.from_mo(mf.mol, f"{formula}.molden", mf.mo_coeff)
        for i in range(min(5, mf.mo_coeff.shape[1])):
            cubegen.orbital(
                mf.mol, 
                'originated.ch4.cube', 
                mf.mo_coeff[:,0],  
                nx=60, ny=60, nz=60,
                resolution=0.2,    
                margin=5.0         
            )  
            cubegen.density(mf.mol, 'linol_acid.cube', mf.make_rdm1()) 
            print(f"{formula}_mo_*.cube")
        
    except Exception as e:
        print(e)