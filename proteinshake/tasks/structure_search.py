import itertools
from collections import defaultdict
from functools import cached_property

from scipy.stats import spearmanr
import numpy as np
from sklearn import metrics
from sklearn.model_selection import train_test_split

from proteinshake.datasets import TMAlignDataset
from proteinshake.tasks import Task

class StructureSearchTask(Task):
    """ Retrieve similar proteins to a query.
    Evaluation is cast in the setting of recommender systems where we wish
    to retrieve 'relevant' documents from a large pool of documents.
    Here, a protein is a document and the relevant ones are all proteins
    with a minimum similarity to the query protein.

    """

    DatasetClass = TMAlignDataset
    
    type = 'Retrieval'
    input = 'Protein'
    output = 'Similar Proteins'
    default_metric = 'Precision@k'

    def __init__(self, min_sim=0.8, *args, **kwargs):
        self.min_sim = min_sim
        super().__init__(*args, **kwargs)

    @cached_property
    def targets(self):
        """ Precompute the set of similar proteins for each query """
        targets = {}
        for query in self.proteins:
            targets[query['protein']['ID']] = [c['protein']['ID'] for c in self.proteins if self.dataset.lddt(query['protein']['ID'], c['protein']['ID']) >= self.min_sim]
        return targets

    def target(self, protein):
        """ The target for a protein is a list of proteins deemed 'relevant' according to `self.min_sim`.
        """
        return self.targets[protein['protein']['ID']]

    def precision_at_k(self, y_true, y_pred, k):
        return len(set(y_pred[:k]).intersection(set(y_true))) / len(y_pred)

    def recall_at_k(self, y_true, y_pred, k):
        return len(set(y_pred[:k]).intersection(set(y_true))) / len(y_true)

    def evaluate(self, y_true, y_pred, k=5):
        """ Retrieval metrics.

        Arguments
        -----------
        y_pred:
            List of indices of items (hits) in the dataset for a query.
        """
        return {
            'Precision@k': np.mean([self.precision_at_k(yt, yp, k) for yt, yp in zip(y_true, y_pred)]),
            'Recall@k': np.mean([self.recall_at_k(yt, yp, k) for yt, yp in zip(y_true, y_pred)]),
        }

    def dummy(self):
        import random
        pred = []
        ids = [p['protein']['ID'] for p in self.proteins] 
        for query in self.proteins[self.test_index]:
            targets = self.target(query)
            pred.append(random.sample(ids, len(targets)))
        return pred
