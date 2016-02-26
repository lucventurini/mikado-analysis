#!/usr/bin/env python3

# import os
import argparse
import csv
# import bed12
# from sklearn.cross_validation import cross_val_score
import pickle
from sklearn.ensemble import RandomForestClassifier
import numpy as np
# import matplotlib.pyplot as plt
from mikado_lib.scales.resultstorer import ResultStorer
from mikado_lib.loci_objects import Transcript
#
# from sklearn.datasets import make_classification
# from sklearn.ensemble import ExtraTreesClassifier


class MetricEntry:

    metrics = [_ for _ in Transcript.get_available_metrics() if
               _ not in ["tid", "parent", "score", "blast_score", "snowy_blast_score"]]

    def __init__(self, row):

        self.__tid = row["tid"]

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
    def tid(self):
        return self.__tid


class TabEntry:
    chrom = ""
    start = 0
    end = 0
    left = 0
    right = 0
    strand = "?"
    # M1 = "N"
    M2 = 0
    M3 = 0
    M4 = 0
    # M5 = 0
    # M6 = 0
    # M7 = 0
    M8 = 0
    M9 = 0
    M10 = 0
    M11 = 0.0
    M12 = 0
    M13 = 0
    M14 = 0
    # M15 = 0.0

    def __init__(self):
        self.data = []

    def __str__(self):

        line = [self.chrom,
                self.start,
                self.end,
                self.left,
                self.right,
                self.strand,
                self.M2,
                self.M3,
                self.M4,
                self.M8,
                self.M9,
                self.M10,
                self.M11,
                self.M12,
                self.M13,
                self.M14]
        return "    ".join([str(_) for _ in line])

    def makeMatrixRow(self):
        return [self.M2,
                self.M3,
                self.M4,
                self.M8,
                self.M9,
                self.M10,
                self.M11,
                self.M12,
                self.M13,
                self.M14]

    @staticmethod
    def features():
        return ["M2-nb-reads",
                "M3-nb_dist_aln",
                "M4-nb_rel_aln",
                "M8-max_min_anc",
                "M9-dif_anc",
                "M10-dist_anc",
                "M11-entropy",
                "M12-maxmmes",
                "M13-hammping5p",
                "M14-hamming3p"]

    @staticmethod
    def featureAt(index):
        return TabEntry.features()[index]

    @staticmethod
    def sortedFeatures(indicies):
        f = TabEntry.features()
        s=[]
        for i in indicies:
            s.append(f[i])
        return s

    @staticmethod
    def nbMetrics():
        return 10

    @staticmethod
    def create_from_tabline(line):

        b = TabEntry()

        parts = line.split("    ")

        b.chrom = parts[2]
        b.start = int(parts[4])
        b.end = int(parts[5])
        b.left = int(parts[6])
        b.right = int(parts[7])
        b.strand = parts[12]

        b.M2 = int(parts[14])
        b.M3 = int(parts[15])
        b.M4 = int(parts[16])
        b.M8 = int(parts[20])
        b.M9 = int(parts[21])
        b.M10 = int(parts[22])
        b.M11 = float(parts[23])
        b.M12 = int(parts[24])
        b.M13 = int(parts[25])
        b.M14 = int(parts[26])

        return b


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


# def loadtab(tabfile):
#
#     bed = list()
#     tab = list()
#     with open(tabfile) as f:
#         # Skip header
#         f.readline()
#
#         for line in f:
#             line.strip()
#             if len(line) > 1:
#                 bed.append(bed12.BedEntry.create_from_tabline(line, False, False))
#                 tab.append(TabEntry.create_from_tabline(line))
#
#     return bed, tab


def main():
    parser = argparse.ArgumentParser("Script to build a random forest decision tree")
    parser.add_argument("-t", "--tmap", help="The TMAP file with the comparison results.")
    parser.add_argument("-m", "--metrics", help="The metrics file.")
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
        score = np.mean([tmap_results[entry.tid].j_f1,
                         tmap_results[entry.tid].n_f1])
        y.append(score)

    clf = RandomForestClassifier(n_estimators=int(len(MetricEntry.metrics)/3),
                                 max_depth=None,
                                 criterion="entropy",
                                 n_jobs=1,
                                 random_state=0)

    clf.fit(X, y)
    importances = clf.feature_importances_
    std = np.std([tree.feature_importances_ for tree in clf.estimators_], axis=0)
    indices = np.argsort(importances)[::-1]

    # Print the feature ranking
    print("Feature ranking:")

    for f in range(X.shape[1]):
        print("{0}, feature {1} ({2})".format(
            f + 1,
            MetricEntry.metrics[f],
            importances[indices[f]]
        ))

    with open("forest.pickle", "wb") as forest:
        pickle.dump(clf, forest)
    # scores = cross_val_score(clf, X, y)
    # scores.mean()

    with open("Xy.pickle", "wb") as xy:
        pickle.dump([X, y], xy)

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
