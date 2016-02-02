import csv
import argparse
import numpy
import sys
from mikado_lib.scales.resultstorer import ResultStorer
from mikado_lib.loci_objects.transcript import Transcript
from mikado_lib.scales.assigner import Assigner
# import multiprocessing
from collections import defaultdict

__author__ = 'Luca Venturini'

def calculate_score(metrics,
                    rescaling,
                    target=None,
                    multiplier=1):
    """
    Private method that calculates a score for each transcript,
    given a target parameter.
    """

    assert rescaling in ("max", "min", "target"), rescaling
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
        scores[tid] = round(score, 2)
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

    if set(_.ccode for _ in assignments.values()) == ("u",):
        return None

    new_assignments = defaultdict(list)
    if len(locus) == 1:
        raise AssertionError(locus)

    for res in assignments.values():
        for gid in res.ref_gene:
            new_assignments[gid].append(res)

    best_res = set()

    for gid in new_assignments:
        ordered = sorted(new_assignments[gid], key=Assigner.result_sorter,
                         reverse=True)
        best = ordered[0]
        best_res.add(best.tid)
        if len(ordered) > 1:
            for num in range(1, len(ordered)):
                if ordered[num].ccode == best.ccode:
                    best_res.add(ordered[num].tid)
                else:
                    break
        if any([_.ccode == ("=",) for _ in new_assignments[gid]]):
            equals = set([_.tid for _ in new_assignments[gid] if
                          _.ccode == ("=",)])
            assert set.intersection(equals, best_res) != set(), (equals, best_res)

    scores_dict = calculate_scores(locus, conf)
    final_dict = dict((tid, round(sum(scores_dict[tid]),2)) for tid in scores_dict)

    best_score = sorted(set(final_dict.values()), reverse=True)[0]

    best = set(_ for _ in final_dict if final_dict[_] == best_score)
    assert len(best) > 0
    # print(locus_name, *sorted(final_dict.items()), sep="\n\t")
    inters = tuple(sorted(list(set.intersection(best, best_res))))
    if len(set(final_dict.values())) == 1:
        # If the score is ininfluent, discard it
        inters = tuple()

    # if set.intersection(best, best_res) != set():
    return [locus_name,
            tuple(sorted(list(best_res))),
            inters]
    # else:
    #     return [locus_name]

    # for tid in scores_dict:
    #     print(tid, sum(scores_dict[tid]))


def main():

    parser = argparse.ArgumentParser("""Script to perform a genetic algo analysis
    of Mikado metrics.""")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("tmap", type=argparse.FileType(),
                        help="The TMAP file to read for scoring.")
    parser.add_argument("metrics", type=argparse.FileType(),
                        help="The Mikado metrics file.")
    args = parser.parse_args()

    results = dict()
    metrics = [_ for _ in Transcript.get_available_metrics() if
               _ not in ("parent", "tid", "score")]

    for row in csv.DictReader(args.tmap, delimiter="\t"):

        for key in ["n_prec", "n_recall", "n_f1",
                    "j_prec", "j_recall", "j_f1",
                    "e_prec", "e_recall", "e_f1", "distance"]:
            if row[key] != "-":
                row[key] = tuple([float(_) for _ in row[key].split(",")])
        res = []
        for key in ResultStorer.__slots__:
            res.append(row[key])
        result = ResultStorer(*res)
        assert result.tid not in results, result.tid
        results[result.tid] = result

    loci = dict()

    metric_arrays = defaultdict(set)

    for row in csv.DictReader(args.metrics, delimiter="\t"):
        if row["parent"] not in loci:
            loci[row["parent"]] = dict()
        loci[row["parent"]][row["tid"]] = dict()
        for key in metrics:
            # Yes I know that eval is bad practice.
            if str.isdigit(row[key][0]):
                row[key] = float(row[key])
            else:
                assert row[key] in ("True", "False")
                if row[key] == "True":
                    row[key] = True
                else:
                    row[key] = False
            loci[row["parent"]][row["tid"]][key] = row[key]
            metric_arrays[key].add(row[key])

    conf = dict()

    key_results = defaultdict(list)
    # pool = multiprocessing.Pool(args.threads)

    loci_num = len(loci)
    for metric in metrics:
        if metric not in conf:
            conf[metric] = dict()
        for tag in ("max", "min"):
            conf[metric][(metric, tag)] = dict()
            conf[metric][(metric, tag)][metric] = dict()
            conf[metric][(metric, tag)][metric]["target"] = None
            conf[metric][(metric, tag)][metric]["rescaling"] = tag
            conf[metric][(metric, tag)][metric]["multiplier"] = 1
        # If we are dealing with something non-boolean ...
        if metric_arrays[metric] not in ({True, False},
                                        {False},
                                        {True}):
            max_val = max(metric_arrays[metric])
            min_val = min(metric_arrays[metric])
            name = "target"
            for prop in numpy.arange(0, 1.1, 0.1):
                target = min_val + (max_val - min_val) * prop
                tag = (name, prop, target)
                conf[metric][(metric, tag)] = dict()
                conf[metric][(metric, tag)][metric] = dict()
                conf[metric][(metric, tag)][metric]["target"] = target
                conf[metric][(metric, tag)][metric]["rescaling"] = name
                conf[metric][(metric, tag)][metric]["multiplier"] = 1

        final_res = dict()
        for key in conf[metric]:
            # print(conf[metric][key])
            for locus in sorted(loci):
                if len(loci[locus]) == 1:
                    continue
                locus_results = dict((_, results[_]) for _ in loci[locus])
                key_results[key].append(evaluate(locus, loci[locus], locus_results,
                                                      conf[metric][key]))
                # job = pool.apply_async()
                # # job = job.get()
                # key_results[key].append(job)
            total = [_ for _ in key_results[key] if _ is not None]
            good = [_ for _ in total if _[-1] != tuple()]
            # good_tmaps = list()
            # for _ in good:
            #     good_tmaps.extend(list(_[1]))
            print(key, len(good), len(good) * 100 / loci_num, file=sys.stderr)
            final_res[key] = [len(good), len(good) * 100 / loci_num]

        best_res = sorted(final_res.keys(), key=lambda key: final_res[key][0],
                          reverse=True)[0]

        print(best_res, final_res[best_res])
        print("Best res", best_res, final_res[best_res], file=sys.stderr)

if __name__ == "__main__":
    main()

