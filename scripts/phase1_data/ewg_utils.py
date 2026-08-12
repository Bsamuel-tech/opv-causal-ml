"""
ewg_utils.py — Shared EWG feature computation module
=====================================================
Centralises Hammett-weighted EWG score computation so
step1 and step2 use identical logic. Any change to
Hammett weights or SMARTS patterns only needs to be
made here.

Author: Samuel Bizimana | JUNIA ISEN
Supervisor: Dr. Kekeli N'KONOU
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# Hammett sigma-para values from literature
HAMMETT = {
    'cyano':    0.66,
    'carbonyl': 0.50,
    'nitro':    0.78,
    'halogen':  0.23,
    'ester':    0.45,
    'sulfonyl': 0.72
}

# SMARTS patterns
SMARTS = {
    'cyano':    '[C]#[N]',
    'carbonyl': '[C]=[O]',
    'nitro':    '[N+](=O)[O-]',
    'halogen':  '[F,Cl,Br,I]',
    'ester':    '[OX2][CX3](=[OX1])',
    'sulfonyl': '[SX4](=[OX1])(=[OX1])'
}

def compute_ewg_features(smi):
    """
    Compute EWG count, Hammett-weighted EWG score,
    and structural confounders from a SMILES string.

    Returns dict with keys:
        ewg_count, ewg_weighted, mol_weight,
        n_arom_rings, conj_length, logp
    or None if SMILES is invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None

        # Count each EWG type
        counts = {}
        for name, smarts in SMARTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            counts[name] = len(mol.GetSubstructMatches(pattern))

        ewg_count    = sum(counts.values())
        ewg_weighted = sum(counts[k] * HAMMETT[k] for k in counts)

        # Structural confounders
        mol_weight   = Descriptors.MolWt(mol)
        n_arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        logp         = Descriptors.MolLogP(mol)

        # Conjugation length: count conjugated bonds
        conj_bonds = 0
        for bond in mol.GetBonds():
            if bond.GetIsConjugated():
                conj_bonds += 1
        conj_length = conj_bonds

        return {
            'ewg_count':    ewg_count,
            'ewg_weighted': ewg_weighted,
            'mol_weight':   mol_weight,
            'n_arom_rings': n_arom_rings,
            'conj_length':  conj_length,
            'logp':         logp
        }
    except:
        return None