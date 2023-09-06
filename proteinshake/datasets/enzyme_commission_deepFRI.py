import json

from proteinshake.datasets import RCSBDataset

class EnzymeCommissionDataset(RCSBDataset):
    """ Enzymes with annotated enzyme commission (EC) numbers.

    .. admonition:: Please cite

        Berman, H M et al. “The Protein Data Bank.” Nucleic acids research vol. 28,1 (2000): 235-42. doi:10.1093/nar/28.1.235

    .. admonition:: Source

        Raw data was obtained and modified from `RCSB Protein Data Bank <https://www.rcsb.org/>`_, originally licensed under `CC0 1.0 <https://creativecommons.org/publicdomain/zero/1.0/>`_.


    .. list-table:: Dataset stats
        :widths: 100
        :header-rows: 1

        * - # proteins
        * - 15603 


    .. list-table:: Annotations
        :widths: 25 35 45
        :header-rows: 1

        * - Attribute
        - Key
        - Sample value
        * - Enzyme Commission
        - :code:`protein['protein']['EC']`
        - :code:`'2.7.7.4'`

    """

    description = 'Enzymes'

    def __init__(self,  **kwargs):
        """

        Args:
            query: REST-API query.

        """
        super().__init__(query=query, **kwargs)

    def download(self):
        # get the annots
        download_url(f'https://github.com/flatironinstitute/DeepFRI/blob/master/preprocessing/data/nrPDB-EC_2020.04_annot.tsv', f'{self.root}/raw/annots.tsv')
        self.ec_annots = pd.read_csv(f'{self.root}/raw/annots.tsv')
        ids = list(self.ec_annots['FA-PDBID'].unique())

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

