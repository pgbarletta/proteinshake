# -*- coding: utf-8 -*-
import glob
import os
import gzip

from biopandas.pdb import PandasPdb
from unittest.mock import patch

from proteinshake.datasets.dataset import Dataset, AA_THREE_TO_ONE
from proteinshake.utils import download_url, progressbar, extract_tar

EXTENDED_AA_THREE_TO_ONE = {
    **AA_THREE_TO_ONE,
    'CYZ': 'C',
    'CYX': 'C',
    'HIP': 'H',
    'HID': 'H',
    'HIE': 'H',
}

class ProteinLigandDecoysDataset(Dataset):
    """ Proteins (targets) from `DUDE-Z <https://pubs.acs.org/doi/10.1021/acs.jcim.0c00598>`_ with a list of decoys and active molecules for each.
    Each molecule is encoded as a SMILES string, meant to be used in a virtual screen setting.
    In this setting a model is given a protein and a ligand and outputs a score reflecting the likelihood
    that the given molecule is a binder. Then, this score is used to sort the union of all the ligands and
    decoys. A good model places true ligands at the top of this list. This is known as enrichment factor
    analysis.
    `data source <https://dudez.docking.org/>`_.


    .. list-table:: Dataset stats
       :widths: 100
       :header-rows: 1

       * - # proteins
       * - 38


   .. list-table:: Annotations
      :widths: 25 35 45
      :header-rows: 1

      * - Attribute
        - Key
        - Sample value
      * - Non-binders SIMLES
        - :code:`protein['protein']['deocys_smiles']`
        - :code:`['O=C(CSc1nnc(COc2ccccc2)o1)NC1CCCCC1', 'C[N@H+]1CC[C@@](N)(C(=O)NC[C@@H]2CC[C@@H](C[NH3+])CC2)C1',..]`
      * - Non-binders identifiers
        - :code:`protein['protein']['decoys_ids']`
        - :code:`['ZINC000000087599', 'ZINC000648138664',..]`
      * - Binders SIMLES
        - :code:`protein['protein']['ligands_smiles']`
        - :code:`['CC1=CC2=C(NC(=O)[C@H](CC3CC3)C2)C(=O)N1CC(=O)NCC1=CC=C(N)N=C1C', 'ClC1=CC=CC(CC2=NC3=C(NCCC4CCCC[NH2+]4)N=CC=C3O2)=C1',..]`
      * - Binders identifiers 
        - :code:`protein['protein']['ligands_ids']`
        - :code:`['CHEMBL10785', 'CHEMBL439678', 'CHEMBL278985',..]`
      * - Pfam accession code
        - :code:`protein['protein']['Pfam']`
        - ``['PF00102']``


    """

    @patch('proteinshake.datasets.dataset.AA_THREE_TO_ONE', EXTENDED_AA_THREE_TO_ONE)
    def pdb2df(self, path):
        return super().pdb2df(path)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*.pdb')[:self.limit]

    def get_id_from_filename(self, filename):
        return filename.split(".")[0]

    def download(self):
        targets  = ['AA2AR', 'ABL1', 'ACES', 'ADA', 'ADRB2', 'AMPC', 'ANDR', 'CSF1R', 'CXCR4', 'DEF', 'DRD4', 'EGFR', 'FA7', 'FA10', 'FABP4', 'FGFR1', 'FKB1A', 'GLCM', 'HDAC8', 'HIVPR', 'HMDH', 'HS90A', 'ITAL', 'KITH', 'KIT', 'LCK', 'MAPK2', 'MK01', 'MT1', 'NRAM', 'PARP1', 'PLK1', 'PPARA', 'PTN1', 'PUR2', 'RENI', 'ROCK1', 'SRC', 'THRB', 'TRY1', 'TRYB1', 'UROK', 'XIAP']

        for target_id in progressbar(targets, desc='Downloading', verbosity=self.verbosity):
            # grab receptor
            download_url(f"https://dudez.docking.org/DOCKING_GRIDS_AND_POSES/{target_id}/rec.crg.pdb", f"{self.root}/raw/files/", verbosity=self.verbosity)
            os.rename(f'{self.root}/raw/files/rec.crg.pdb', f'{self.root}/raw/files/{target_id}.pdb')
            # grab ligands
            download_url(f"https://dudez.docking.org/property_matched/{target_id}/ligands.tgz", f"{self.root}/raw/files/ligands_{target_id}.tgz", verbosity=self.verbosity)
            extract_tar(f'{self.root}/raw/files/ligands_{target_id}.tgz', f'{self.root}/raw/ligands_{target_id}', verbosity=self.verbosity)
            # grab decoys
            download_url(f"https://dudez.docking.org/property_matched/{target_id}/decoys.tgz", f"{self.root}/raw/files/decoys_{target_id}.tgz", verbosity=self.verbosity)
            extract_tar(f'{self.root}/raw/files/decoys_{target_id}.tgz', f'{self.root}/raw/decoys_{target_id}', verbosity=self.verbosity)



    def db2_parse(self, path, gzipped=True, smiles_only=True):
        """ Parse a DB2 file. Ugly for now.
        """
        data = {}
        if gzipped:
            with gzip.open(path, "rb") as db:
                for i, line in db:
                    if i == 0:
                        data['id'] = line.split()[1]
                    if i == 2:
                        data['smiles'] = line.split()[1]
                    if i > 2:
                        break


        return data

    def add_protein_attributes(self, protein):
        """ We annotate each protein with a list of decoys and a list of active SMILES strings and molecule IDs.
        """
        target = protein['protein']['ID']
        
        for mode in ['decoys', 'ligands']:
            smiles, ids = [], []
            for molfile in os.listdir(f"{self.root}/raw/files/{mode}_{target}/"):
                mol_data = self.db2_parse(f"{self.root}/raw/files/{mode}_{taret}/{molfile}")
                smiles.append(mol_data['smiles'])
                ids.append(mol_data['id'])

            protein['protein'][f'{mode}_smiles'] = smiles
            protein['protein'][f'{mode}_ids'] = ids
            protein['protein'][f'num_{mode}'] = len(ids)

        protein['protein']['num_mols'] = protein['protein']['num_ligands'] + protein['protein']['num_decoys']
        return protein
