import numpy as np
from sklearn import metrics

from proteinshake.datasets import ProteinProteinInterfaceDataset
from proteinshake.tasks import Task
from proteinshake.transforms import CenterTransform, RandomRotateTransform, Compose

class ProteinProteinInterfaceTask(Task):
    """ Identify the binding residues of a protein-protein complex. This is a residue-level binary classification.

    NOTE: To make this an interesting task, the loader has to
    split the protein into its chains so that the model only sees
    one chain at a time. You should also pass the following transforms to the dataset
    :meth:`proteinshake.transforms.CenterTransform` and :meth:`proteinshake.transforms.RandomRotateTransform` when using representations that are aware of the atomic coordinates.

    NOTE: This task is currently in beta.


    .. code-block:: python

        >>> from proteinshake.tasks import ProteinProteinInterfaceTask
        >>> from proteinshake.transforms import CenterTransform, RandomRotateTransform
        >>> ta = ProteinProteinInterfaceTask()
        >>> data = ta.dataset.to_voxel(transforms=[CenterTransform(), RandomRotateTransform()).torch()
    """

    DatasetClass = ProteinProteinInterfaceDataset

    @property
    def task_in(self):
        return ('residue', 'residue')

    @property
    def task_type(self):
        return ('residue', 'binary')

    @property
    def task_out(self):
        return ('binary')

    @property
    def out_dim(self):
        return (1)

    def dummy_output(self):
        import random
        return [random.randint(0, 1) for p in self.test_targets]

    def compute_targets(self):
        # compute targets (e.g. for scaling)
        self.train_targets = [p for i in self.train_index for p in self.target(self.proteins[i])]
        self.val_targets = [p for i in self.val_index for p in self.target(self.proteins[i])]
        self.test_targets = [p for i in self.test_index for p in self.target(self.proteins[i])]

    def target(self, protein_1, protein_2):
        chain_1 = protein_1['residue']['chain_id'][0]
        chain_2 = protein_2['residue']['chain_id'][0]

        contacts = np.zeros((len(protein_1['protein']['sequence']), len(protein_2['protein']['sequence'])))
        inds = self._interfaces[pdbid][chain_1][chain_2]
        np.put(contacts, np.ravel_multi_index(np.transpose(inds), contacts.shape), 1)
        return contacts

    @property
    def default_metric(self):
        return 'average_precision'

    def evaluate(self, y_true, y_pred):
        """ Evaluate performance of an interface classifier.
        """
        return {
            'auc_roc': metrics.roc_auc_score(y_true, y_pred),
            'average_precision': metrics.average_precision_score(y_true, y_pred),
        }

    def to_graph(self, *args, **kwargs):
        self.dataset = self.dataset.to_graph(*args, **kwargs, transform=Compose([CenterTransform(), RandomRotateTransform()]))
        return self

    def to_point(self, *args, **kwargs):
        self.dataset = self.dataset.to_point(*args, **kwargs, transform=Compose([CenterTransform(), RandomRotateTransform()]))
        return self

    def to_voxel(self, *args, **kwargs):
        self.dataset = self.dataset.to_voxel(*args, **kwargs, transform=Compose([CenterTransform(), RandomRotateTransform()]))
        return self
