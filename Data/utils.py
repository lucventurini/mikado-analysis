import yaml
from collections import OrderedDict, defaultdict
from operator import itemgetter
import os


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

