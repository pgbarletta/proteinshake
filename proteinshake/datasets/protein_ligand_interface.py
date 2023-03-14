# -*- coding: utf-8 -*-
import glob
import os
import os.path as osp

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem

import numpy as np

from proteinshake.datasets import Dataset
from proteinshake.utils.pdbbind import parse_pdbbind_PL_index
from proteinshake.utils import extract_tar, download_url

RDLogger.DisableLog('rdApp.*') # disable warnings

class ProteinLigandInterfaceDataset(Dataset):
    """Proteins bound to small molecules from PDBBind with binding site, ligand and affinity information.



    Parameters
    ----------
    root: str
        Root directory where the dataset should be saved.
    name: str
        The name of the dataset.
    version: str
        PDBBind version to use.

    .. list-table:: Dataset stats
       :widths: 100
       :header-rows: 1

       * - # proteins
       * - 4642


   .. list-table:: Annotations
      :widths: 20 55 25
      :header-rows: 1

      * - Attribute
        - Key
        - Sample value
      * - Dissociation constant (kd)
        - :code:`protein['protein']['kd']`
        - :code:`77.0`
      * - Affinity
        - :code:`protein['protein']['neglog_aff']`
        - :code:`4.11000`
      * - Resolution (Angstroms)
        - :code:`protein['protein']['resolution']`
        - :code:`2.20`
      * - Year solved
        - :code:`protein['protein']['year']`
        - :code:`2016`
      * - Ligand identifier (PDB code)
        - :code:`protein['protein']['ligands_id']`
        - :code:`IEE`
      * - Ligand SMILES 
        - :code:`protein['protein']['ligand_smiles']`
        - :code:`'Cc1ccc(CNc2cc(Cl)nc(N)n2)cc1'`
      * - Molecular ingerprints
        - :code:`protein['protein']['fp_maccs']`, :code:`protein['protein']['fp_morgan_r2']`
        - :code:`'[..,0, 0, 1, 0, 1, 0, 0, 0,..]`
      * - Molecular ingerprints
        - :code:`protein['protein']['fp_maccs']`, :code:`protein['protein']['fp_morgan_r2']`
        - :code:`'[..,0, 0, 1, 0, 1, 0, 0, 0,..]`
      * - Binding site (1 if in binding site, 0 else)
        - :code:`protein['residue']['binding_site']`
        - :code:`'[..,0, 0, 1, 0, 1, 0, 0, 0,..]`





    """


    def __init__(self, version='2020', **kwargs):
        self.version = version

        super().__init__(**kwargs)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*/*_protein.pdb')[:self.limit]

    def get_id_from_filename(self, filename):
        return filename[:4]

    def download(self):
        download_url(f'https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v{self.version}_refined.tar.gz', f'{self.root}/raw')
        extract_tar(f'{self.root}/raw/PDBbind_v{self.version}_refined.tar.gz', f'{self.root}/raw')
        os.rename(f'{self.root}/raw/refined-set', f'{self.root}/raw/files')

        self.index_data = parse_pdbbind_PL_index(osp.join(self.root,
                                                          "raw",
                                                          "files",
                                                          "index",
                                                          f"INDEX_refined_data.{self.version}")
                                                )

    def add_protein_attributes(self, protein):
        pocket = self.pdb2df(f'{self.root}/raw/files/{protein["protein"]["ID"]}/{protein["protein"]["ID"]}_pocket.pdb')
        ligand = Chem.MolFromMolFile(osp.join(self.root,
                                              'raw',
                                              'files',
                                              protein['protein']['ID'],
                                              f'{protein["protein"]["ID"]}_ligand.sdf')
                                     )

        if ligand is None:
            return None
        smiles = Chem.MolToSmiles(ligand)
        fp_morgan = list(map(int, AllChem.GetMorganFingerprintAsBitVect(ligand, 2, nBits=1024).ToBitString()))
        fp_maccs = list(map(int, MACCSkeys.GenMACCSKeys(ligand).ToBitString()))

        pocket_res = np.array(pocket.loc[pocket['atom_type'] == 'CA']['residue_number'].tolist())
        protein_res = np.array(protein['residue']['residue_number']).reshape(-1, 1)

        pocket_atom = np.array(pocket['residue_number'].tolist())
        protein_atom = np.array(protein['atom']['residue_number']).reshape(-1, 1)

        is_site_res = np.zeros_like(protein_res)
        is_site_atom  = np.zeros_like(protein_atom)

        is_site_res[(pocket_res == protein_res).sum(axis=1).nonzero()] = 1.
        is_site_atom[(pocket_atom == protein_atom).sum(axis=1).nonzero()] = 1.


        protein['residue']['binding_site'] = list(map(int, is_site_res.squeeze()))
        protein['atom']['binding_site'] = list(map(int, is_site_atom.squeeze()))

        bind_data = self.index_data[protein['protein']['ID']]
        protein['protein']['kd'] = bind_data['kd']['value']
        protein['protein']['neglog_aff'] = bind_data['neglog_aff']
        protein['protein']['resolution'] = bind_data['resolution']
        protein['protein']['year'] = bind_data['date']
        protein['protein']['ligand_id'] = bind_data['ligand_id']
        protein['protein']['ligand_smiles'] = smiles

        protein['protein']['fp_maccs'] = fp_maccs
        protein['protein']['fp_morgan_r2'] = fp_morgan

        return protein

    def describe(self):
        desc = super().describe()
        desc['property'] = "Small Mol. Binding Site (residue-level)"
        desc['values'] = 2
        desc['type'] = 'Binary'
        return desc
