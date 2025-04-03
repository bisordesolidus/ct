from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from pyscf import gto, scf, tools
from pyscf.tools import cubegen
#import matplotlib.pyplot as plt
from IPython.display import display

def formula_to_3d_geometry(formula):
    """Конвертирует формулу в 3D-геометрию с помощью RDKit."""
    # Создаем молекулу из формулы (нужна SMILES-строка, поэтому пример для воды)
    smiles = formula 
    """{
        "H2O": "O",
        "C18H32O2": "CC=CCC=CCC=CCCCCCCCCC(=O)O",  # α-линоленовая кислота
    }.get(formula, "")"""
    
    if not smiles:
        raise ValueError(f"Формула {formula} не найдена в базе SMILES.")
    
    mol = Chem.MolFromSmiles(smiles)
    for bond in mol.GetBonds():
        print(f"Связь: {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()}, Тип: {bond.GetBondType()}")
    mol = Chem.AddHs(mol)  # Добавляем водороды
    AllChem.EmbedMolecule(mol)  # 3D-конформация
    AllChem.MMFFOptimizeMolecule(mol)  # Оптимизация геометрии
    #tools.cubegen.orbital(mol, 'orbital.cube', mo_coeff, resolution=0.2)
    display(Draw.MolToImage(mol, size=(400, 400)))
    # Извлекаем координаты атомов
    atoms = []
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        atoms.append(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")
    
    return "\n".join(atoms)

def calculate_electronic_structure(geometry):
    #Расчет электронной структуры в PySCF.
    mol = gto.M(
        atom=geometry,
        basis="6-31g*",
        spin=0,
        verbose=4
    )
    mf = scf.RHF(mol).run()
    return mf

"""def calculate_electronic_structure(geometry):
    #Расчет электронной структуры с собственной реализацией RHF
    mol = gto.M(
        atom=geometry,
        basis="6-31g*",
        spin=0,
        verbose=4
    )
    mf = my.MyRHF(mol)
    mf.kernel()
    return mf"""

# Пример использования
if __name__ == "__main__":
    formula = "CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O"
    
    try:
        # 1. Генерация геометрии
        geometry = formula_to_3d_geometry(formula)
        print("3D-геометрия:\n", geometry)
        
        # 2. Расчет электронной структуры
        mf = calculate_electronic_structure(geometry)
        print(f"Энергия Хартри-Фока: {mf.e_tot:.6f} Hartree")
        
        # 3. Сохранение орбиталей
        #tools.molden.from_mo(mf.mol, f"{formula}.molden", mf.mo_coeff)
        #print(f"Результаты сохранены в {formula}.molden (откройте в Avogadro/VMD)")
        for i in range(min(5, mf.mo_coeff.shape[1])):
    # Явно задаем сетку 60x60x60 точек с шагом 0.2 ангстрема
            cubegen.orbital(
                mf.mol, 
                'originated.ch4.cube', 
                mf.mo_coeff[:,0],  # Первая орбиталь
                nx=60, ny=60, nz=60,
                resolution=0.2,    # Шаг сетки (Å)
                margin=5.0         # Отступ от края молекулы (Å)
            )  
            cubegen.density(mf.mol, 'linol_acid.cube', mf.make_rdm1()) 
            print(f"{formula}_mo_*.cube")
        
    except Exception as e:
        print(e)