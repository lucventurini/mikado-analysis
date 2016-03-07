#!/usr/bin/env python3

import argparse
import csv
import pickle
import collections
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import operator
from Mikado.scales.resultstorer import ResultStorer
from Mikado.loci_objects import Transcript
import re
# from sklearn.cross_validation import cross_val_score
# from sklearn.datasets import make_classification
# from sklearn.ensemble import ExtraTreesClassifier
# import matplotlib.pyplot as plt


class MetricEntry:

    metrics = [_ for _ in Transcript.get_available_metrics() if
               _ not in ["tid", "parent", "score", "snowy_blast_score", "best_bits"]]

    def __init__(self, row):

        self.__tid = row["tid"]
        self.__locus = row["parent"]

        for key in self.metrics:
            if row[key].lower() == "true":
                row[key] = 1.0
            elif row[key].lower() == "false":
                row[key] = 0.0
            else:
                try:
                    row[key] = float(row[key])
                except ValueError as exc:
                    raise ValueError("Invalid value for key {0}: {1}.\n{2}".format(
                        key, row[key],
                        exc))
            setattr(self, key, row[key])

    @property
    def matrix_row(self):
        return [getattr(self, key) for key in self.metrics]

    @property
    def features(self):
        return self.metrics

    @property
    def nb_features(self):
        return len(self.metrics)

    @property
    def tid(self):
        return self.__tid

    @property
    def locus(self):
        return self.__locus


def load_tmap(tmap_file) -> dict:

    """
    Function to serialise the TMAP file into the original ResultStorer objects.
    :param tmap_file:
    :return:
    """

    results = dict()
    with open(tmap_file) as refmap:
        for row in csv.DictReader(refmap, delimiter="\t"):
            vals = []
            for key in ResultStorer.__slots__:
                if key in ("n_prec", "n_recall", "n_f1",
                           "j_prec", "j_recall", "j_f1",
                           "e_prec", "e_recall", "e_f1"):
                    row[key] = tuple(float(_) for _ in row[key].split(","))
                elif key == "distance":
                    if row[key] == "-":
                        row[key] = 5001
                    else:
                        row[key] = int(row[key])
                elif key in ("ccode", "tid", "gid",):
                    row[key] = tuple(row[key].split(","))
                vals.append(row[key])
            result = ResultStorer(*vals)
            results[result.tid[0]] = result

    return results


def load_metrics(metrics_file) -> [MetricEntry]:

    metrics = []
    with open(metrics_file) as metrics_handle:
        reader = csv.DictReader(metrics_handle, delimiter="\t")
        for row in reader:
            del row["score"]
            row = MetricEntry(row)
            metrics.append(row)

    return metrics


def main():
    parser = argparse.ArgumentParser("Script to build a random forest decision tree")
    parser.add_argument("-t", "--tmap", help="The TMAP file with the comparison results.",
                        required=True)
    parser.add_argument("-m", "--metrics", help="The metrics file.",
                        required=True)
    # parser.add_argument("-r", "--reference", required=True, help="The reference BED file to compare against")
    # parser.add_argument("-o", "--output", required=True, help="The output prefix")
    args = parser.parse_args()

    # X should contain a matrix of features derived from the portcullis tab file
    # y should contain the labels (0 not a valid junction, 1 a valid junction).
    # Confirmed with the reference.

    # Load tab file and produce matrix
    # bed, tab = loadtab(args.input)
    tmap_results = load_tmap(args.tmap)

    print("# TMAP results: " + str(len(tmap_results)))

    # Load reference and add labels
    # ref = bed12.loadbed(args.reference, False, False)
    metrics = load_metrics(args.metrics)
    print("# metered transcripts:", len(metrics))

    X = np.zeros((len(metrics), len(MetricEntry.metrics)))
    y = []

    for index, entry in enumerate(metrics):
        X[index] = entry.matrix_row
        if ".orf" in entry.tid:
            tid = re.sub("\.orf[0-9]*$", "", entry.tid)
        else:
            tid = entry.tid

        if entry.tid in tmap_results:
            score = np.mean([tmap_results[entry.tid].j_f1,
                             tmap_results[entry.tid].n_f1])
            y.append(score)
        else:
            y.append(0)

    clf = RandomForestRegressor(n_estimators=int(len(MetricEntry.metrics)/3),
                                max_depth=None,
                                n_jobs=1,
                                random_state=0)

    clf.fit(X, y)
    importances = clf.feature_importances_
    std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)
    # indices = np.argsort(importances)[::-1]

    order = sorted([(MetricEntry.metrics[_], 100 * importances[_]) for _ in range(X.shape[1])],
                   key=operator.itemgetter(1), reverse=True)

    # Print the feature ranking
    print("Feature ranking:")

    for rank, feature in enumerate(order, start=1):
        print("{0}, feature {1} ({2})".format(
            rank, feature[0], feature[1]))

    with open("forest.pickle", "wb") as forest:
        pickle.dump(clf, forest)
    # scores = cross_val_score(clf, X, y)
    # scores.mean()

    with open("Xy.pickle", "wb") as xy:
        pickle.dump([X, y], xy)

    parents = collections.defaultdict(list)
    for m in metrics:
        parents[m.locus].append((m.tid, m.matrix_row))

    for parent, children in parents.items():
        print(parent)
        vals = [_[1] for _ in children]
        proba = clf.predict(vals)
        # scores = clf.predict(vals)
        probs = sorted([(children[_][0], proba[_]) for _ in range(len(children))],
                       key=operator.itemgetter(1), reverse=True)
        print("\t", *probs, sep="\t\n")

    # probs = clf.predict_proba(X)
    # print(*probs)

    # Plot the feature importances of the forest
    # plt.figure()
    # plt.title("Feature importances")
    # plt.bar(range(X.shape[1]), importances[indices],
    #        color="r", yerr=std[indices], align="center")
    # plt.xticks(range(X.shape[1]), TabEntry.sortedFeatures(indices))
    # locs, labels = plt.xticks()
    # plt.setp(labels, rotation=90)
    # plt.xlim([-1, X.shape[1]])
    # plt.tight_layout()
    #
    # plt.savefig(args.output + ".png")

main()
