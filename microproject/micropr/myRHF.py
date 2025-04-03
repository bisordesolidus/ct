from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import numpy as np
from scipy.linalg import eigh
from pyscf import gto, scf, tools
from pyscf.lib import logger
from pyscf.scf.diis import DIIS
from pyscf.tools import cubegen
from IPython.display import display

class MyRHF:
    def __init__(self, mol):
        self.mol = mol
        self.conv_tol = 1e-8
        self.max_cycle = 50
        self.diis_space = 8
        self.mo_energy = None
        self.mo_coeff = None
        self.e_tot = 0.0
        
    def kernel(self):
        mol = self.mol
        S = mol.intor('int1e_ovlp')
        T = mol.intor('int1e_kin')
        V = mol.intor('int1e_nuc')
        H = T + V
        
        X = np.linalg.inv(scf.hf.get_ovlp(mol))
        
        mo_energy, mo_coeff = eigh(H, S)
        nocc = mol.nelectron // 2
        dm = self.make_rdm1(mo_coeff, nocc)
        
        diis = DIIS()
        diis.space = self.diis_space

        old_dm = dm
        for cycle in range(self.max_cycle):
            J, K = self.get_jk(mol, dm)
            F = H + J - K

            F = diis.update(F, dm, S, mo_coeff)
            
            # Ротаан
            mo_energy, mo_coeff = eigh(F, S)
            dm = self.make_rdm1(mo_coeff, nocc)
            
            delta_dm = np.linalg.norm(dm - old_dm)
            if delta_dm < self.conv_tol:
                logger.note(self, 'SCF converged in %d cycles', cycle+1)
                break
            old_dm = dm
            
        self.mo_energy = mo_energy
        self.mo_coeff = mo_coeff
        self.e_tot = self.energy_tot(dm, H, F)
        return self.e_tot
        
    def make_rdm1(self, mo_coeff=None, nocc=None):
        if mo_coeff is None:
            mo_coeff = self.mo_coeff
        if nocc is None:
            nocc = self.mol.nelectron // 2
    
        if mo_coeff is not None:
            occ_mo = mo_coeff[:, :nocc]
            dm = 2 * np.dot(occ_mo, occ_mo.T)
        else:
            nao = self.mol.nao
            dm = np.zeros((nao, nao))
    
        return dm
        
    def get_jk(self, mol, dm):
        """Вычисление кулоновского и обменного членов"""
        eri = mol.intor('int2e')
        J = np.einsum('pqrs,rs->pq', eri, dm)
        K = np.einsum('prqs,rs->pq', eri, dm)
        return J, K
        
    def energy_tot(self, dm, h1e, fock):
        """Вычисление полной энергии"""
        return np.sum(dm * (h1e + fock)) / 2 + self.mol.energy_nuc()