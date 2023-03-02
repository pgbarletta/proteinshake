import numpy as np
from sklearn import metrics

from proteinshake.datasets import ProteinProteinInterfaceDataset
from proteinshake.tasks import Task
from proteinshake.transforms import CenterTransform
from proteinshake.transforms import RandomRotateTransform

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

    def __init__(self, root='data', resolution='residue', *args, **kwargs):
        dataset = ProteinProteinInterfaceDataset(root=root, transforms=[CenterTransform(), RandomRotateTransform()])
        super().__init__(dataset, *args, **kwargs)

    @property
    def task_type(self):
        return "binary-classification"

    def target(self, protein):
        return protein['residue']['is_interface']

    def evaluate(self, pred):
        """ Evaluate performance of an interface classifier.

        pred: list
            One list for each protein in the test set with a probability for each residue in the protein.

        """
        labels = [self.target(p) for p in \
                  self.proteins[self.test_index]]
        labels = np.hstack(labels)
        pred = np.hstack(pred)
        return {
                'auc-roc': metrics.roc_auc_score(labels, pred),
                'average precision': metrics.average_precision_score(labels, pred),
                }

if __name__ == "__main__":
    import random
    ta = ProteinProteinInterfaceTask(use_precomputed=False).to_graph(eps=8).pyg()

    preds = []
    for p in ta.proteins[ta.test_index]:
        preds.append([random.random() for _ in range(len(p['residue']['residue_number']))])

    print(ta.evaluate(preds))
