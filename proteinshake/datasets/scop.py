import pandas as pd
from joblib import Parallel, delayed

from tqdm import tqdm

from proteinshake.datasets import RCSBDataset
from proteinshake.utils import download_url, progressbar

class SCOPDataset(RCSBDataset):
    """ Proteins with annotated SCOP class.

    .. admonition:: Please cite

      Murzin, Alexey G., et al. "SCOP: a structural classification of proteins database for the investigation of sequences and structures." Journal of molecular biology 247.4 (1995): 536-540.

      Berman, H M et al. “The Protein Data Bank.” Nucleic acids research vol. 28,1 (2000): 235-42. doi:10.1093/nar/28.1.235

    .. admonition:: Source

      Raw data was obtained and modified from `RCSB Protein Data Bank <https://www.rcsb.org/>`_, originally licensed under `CC0 1.0 <https://creativecommons.org/publicdomain/zero/1.0/>`_.


    .. list-table:: Dataset stats
      :widths: 100
      :header-rows: 1

      * - # proteins
      * - 10066


    .. list-table:: Annotations
      :widths: 25 45 35
      :header-rows: 1

      * - Attribute
      - Key
      - Sample value
      * - Protein type
      - :code:`protein['protein']['SCOP-TP']`
      - :code:`'1'`
      * - Protein class 
      - :code:`protein['protein']['SCOP-CL']`
      - ``'1000000'``
      * - Superfamily 
      - :code:`protein['protein']['SCOP-SF']`
      - ``'3000001'``
      * - Family
      - :code:`protein['protein']['SCOP-FA']`
      - ``'4002873'``


    """

    def _parse_scop(self, path):
        names = ['FA-DOMID', 'FA-PDBID', 'FA-PDBREG', 'FA-UNIID', 'FA-UNIREG', 'SF-DOMID', 'SF-PDBID', 'SF-PDBREG', 'SF-UNIID', 'SF-UNIREG', 'SCOPCLA']
        df = pd.read_csv(path, sep=' ', comment='#', names=names, dtype=str)
        return {k: dict([cla.split("=") for cla in v.split(",")]) for k,v in zip(df['FA-PDBID'], df['SCOPCLA'])}

    def download(self):
        # get the annots
        download_url(f'http://scop.mrc-lmb.cam.ac.uk/files/scop-cla-latest.txt', f'{self.root}/raw/scop.txt')
        self.scop = self._parse_scop(f'{self.root}/raw/scop.txt')
        ids = list(self.scop['FA-PDBID'].unique())

        # get the proteins
        if self.n_jobs == 1:
            print('Warning: Downloading an RCSB dataset with use_precompute = False is very slow. Consider increasing n_jobs.')
        ids = ids[:self.limit] # for testing

        failed = Parallel(n_jobs=self.n_jobs)(delayed(self.download_from_rcsb)(id) for id in progressbar(ids, desc='Downloading PDBs', verbosity=verbosity))
        failed = [f for f in failed if not f is True]
        if len(failed)>0:
            print(f'Failed to download {len(failed)} PDB files.')

    def add_protein_attributes(self, protein):
        """ We annotate the protein with the scop classifications at each level.

        SCOPCLA - SCOP domain classification. The abbreviations denote: TP=protein type, CL=protein class, CF=fold, SF=superfamily, FA=family
        """
        protein_id = protein['protein']['ID'].upper()
        if not protein_id in self.scop: return None
        for cla, val in self.scop[protein_id].items():
            protein['protein']['SCOP-' + cla] = val
        return protein
