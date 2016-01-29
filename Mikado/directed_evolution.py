__author__ = 'Luca Venturini'

from mikado_lib.scales.resultstorer import ResultStorer
from mikado_lib.loci_objects.transcript import Transcript
from mikado_lib.scales.assigner import Assigner
import csv
import argparse
import multiprocessing
from collections import defaultdict


def calculate_score(metrics, rescaling, target=None, multiplier=1):
    """
    Private method that calculates a score for each transcript,
    given a target parameter.
    :param param:
    :return:
    """

    assert rescaling in ("max", "min", "target")
    assert multiplier != 0

    metric_vals = [_[1] for _ in metrics]

    if rescaling == "target":
        denominator = max(abs(x - target) for x in metric_vals)
    else:
        try:
            denominator = (max(metric_vals) - min(metric_vals))
        except TypeError as exc:
            raise TypeError("{}\n{}".format(exc, metric_vals))
    if denominator == 0:
        denominator = 1

    scores = dict()
    for tid, tid_metric in metrics:
        score = 0
        # if ("filter" in self.json_conf["scoring"][param] and
        #         self.json_conf["scoring"][param]["filter"] != {}):
        #     check = self.evaluate(tid_metric, self.json_conf["scoring"][param]["filter"])

        if rescaling == "target":
            score = 1 - abs(tid_metric - target) / denominator
        else:
            if min(metrics) == max(metrics):
                score = 1
            elif rescaling == "max":
                score = abs((tid_metric - min(metric_vals)) / denominator)
            elif rescaling == "min":
                score = abs(1 - (tid_metric - min(metric_vals)) / denominator)

        score *= multiplier
        # self.scores[tid][param] = round(score, 2)
        scores[tid] = score
    return scores


def calculate_scores(locus, conf):

    scores_dict = defaultdict(list)
    metrics = dict()

    for param in conf:
        rescaling = conf[param]["rescaling"]
        target = conf[param]["target"]
        multiplier = conf[param]["multiplier"]
        metrics[param] = [(tid, locus[tid][param]) for tid in locus]
        scores = calculate_score(metrics[param],
                                 rescaling,
                                 target,
                                 multiplier)
        for tid in locus:
            assert tid in scores, (tid, scores)
            scores_dict[tid].append(scores[tid])

    return scores_dict


def evaluate(locus_name, locus, assignments, conf):

    scores_dict = calculate_scores(locus, conf)
    print(locus_name, len(locus))
    for tid in scores_dict:
        print(tid, sum(scores_dict[tid]))
    print()

    

def main():

    parser = argparse.ArgumentParser("""Script to perform a genetic algo analysis
    of Mikado metrics.""")
    parser.add_argument("tmap", type=argparse.FileType(),
                        help="The TMAP file to read for scoring.")
    parser.add_argument("metrics", type=argparse.FileType(),
                        help="The Mikado metrics file.")
    args = parser.parse_args()

    results = dict()
    metrics = [_ for _ in Transcript.get_available_metrics() if
               _ not in ("parent", "tid", "score")]

    for num, row in enumerate(csv.reader(args.tmap, delimiter="\t")):
        if num == 1:
            continue  # header
        result = ResultStorer(*row)
        assert result.tid not in results
        results[result.tid] = result

    loci = dict()

    for row in csv.DictReader(args.metrics, delimiter="\t"):
        if row["parent"] not in loci:
            loci[row["parent"]] = dict()
        loci[row["parent"]][row["tid"]] = dict()
        for key in metrics:
            # Yes I know that eval is bad practice.
            loci[row["parent"]][row["tid"]][key] = eval(row[key])

    conf = dict()
    conf["blast_score"] = dict()
    conf["blast_score"]["rescaling"] = "max"
    conf["blast_score"]["target"] = None
    conf["blast_score"]["multiplier"] = 1

    for locus in loci:
        evaluate(locus, loci[locus], results, conf)

if __name__ == "__main__":
    main()