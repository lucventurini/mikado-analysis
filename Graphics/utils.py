import yaml
from collections import OrderedDict, defaultdict
from operator import itemgetter
import os
import re
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore", 'This pattern has match groups')  # Remove annoying useless warnings


def parse_refmaps(orig_stats, filtered_stats, transcripts=False, return_sets=False, return_base=False):

    if orig_stats is None:
        if return_sets is True:
            return [set(), set(), set()]
        else:
            return [-1000] * 3

    orig_refmap = "{}.refmap".format(
        re.sub(".stats$", "", orig_stats))
    orig_tmap = "{}.tmap".format(
        re.sub(".stats$", "", orig_stats))

    filtered_refmap = "{}.refmap".format(
        re.sub(".stats$", "", filtered_stats))

    refmap = pd.read_csv(orig_refmap, delimiter="\t")
    filtered_refmap = pd.read_csv(filtered_refmap, delimiter="\t")
    tmap = pd.read_csv(orig_tmap, delimiter="\t")
    
    if return_base is True:
        if transcripts is True:
            full = filtered_refmap.ref_id.unique().astype(set)
            fused = refmap.ref_id.unique().astype(set)
            missed = filtered_refmap.ref_id.unique().astype(set)
        else:
            full = filtered_refmap.ref_gene.unique().astype(set)
            fused = refmap.ref_gene.unique().astype(set)
            missed = filtered_refmap.ref_gene.unique().astype(set)
    else:
        if transcripts is True:
            full = refmap[refmap.ccode.str.contains("(=|_)", regex=True, na=False)].ref_id.unique().astype(set)
            missed = filtered_refmap[filtered_refmap.ccode.isin((np.nan, "p", "P", "i", "I", "ri", "rI", "X", "x"))].ref_id.unique().astype(set)
            fused = tmap[(tmap.ccode.str.contains("f", na=False)) & (~tmap.ccode.str.contains("(=|_)", regex=True, na=False))].ref_id.unique().astype(set)
        else:
            full = set(refmap[refmap.ccode.str.contains("(=|_)", regex=True, na=False)].ref_gene.unique())
            missed = set(filtered_refmap[filtered_refmap.best_ccode.isin((np.nan, "p", "P", "i", "I", "ri", "rI", "X", "x"))].ref_gene.unique())
            fused = set(tmap[(tmap.ccode.str.contains("f", na=False)) & (~tmap.ccode.str.contains("(=|_)", regex=True, na=False))].ref_gene.unique())
        
        assert isinstance(full, set), type(full)
        assert isinstance(missed, set), type(full)
        if len(set.intersection(full, missed)) > 0:
            # Remove from the missed those that are fully reconstructed. This is due to strange annotations in human
            missed = set.difference(missed, set.intersection(full, missed))

    if return_sets is False:
        res = [len(full), len(missed), len(fused)]
    else:
        res = [full, missed, fused]

    return res


def parse_configuration(configuration, exclude_mikado=False, prefix=None):

    options = yaml.load(configuration)
    configuration.close()

    if "divisions" not in options:
        options["divisions"] = defaultdict(dict)
        options["divisions"]["STAR"]["marker"] = "o"
        options["divisions"]["TopHat"]["marker"] = "^"

    for division in options["divisions"]:
        assert "marker" in options["divisions"][division]

    for method in options["methods"]:
        for key in ("colour", "index", *options["divisions"].keys()):
            if key not in options["methods"][method]:
                print("WARNING: {} not found for {}".format(key.capitalize(), method))
                options["methods"][method][key] = None
        for aligner in options["divisions"]:
            if options["methods"][method][aligner] is None:
                continue
            elif not isinstance(options["methods"][method][aligner], list):
                raise TypeError("Invalid type for aligner {}: {}".format(
                    aligner, type(options["methods"][method][aligner])))
            elif len(options["methods"][method][aligner]) != 2:
                raise ValueError("Invalid number of files specified for {} / {}".format(
                    method, aligner))
            if prefix is not None:
                for index, fname in enumerate(options["methods"][method][aligner]):
                    options["methods"][method][aligner][index] = os.path.join(prefix, fname)

            if any(not os.path.exists(_) for _ in options["methods"][method][aligner]):
                raise OSError("Files not found: {}".format(", ".join(
                    options["methods"][method][aligner])))

    new_methods = OrderedDict()

    for index, method in sorted([(options["methods"][method]["index"], method)
                                 for method in options["methods"]], key=itemgetter(0)):
        if exclude_mikado is True and "mikado" in method.lower():
            continue
        new_methods[method] = options["methods"][method]

    options["methods"] = new_methods

    if len(set(options["methods"][method]["colour"]
               for method in options["methods"])) != len(options["methods"]):
        raise ValueError("Invalid unique number of colours specified!")

    if options["format"] not in ("svg", "png", "tiff"):
        raise ValueError("Invalid output format specified: {}".format(
            options["format"]))

    return options

