"""
Copyright (c) 2019 10x Genomics, Inc. All rights reserved.

parse a csv of parameters and determine inputs tweak parameters to analysis stages
"""

# Imports
import martian
from collections import namedtuple
import os
import csv

from utils import get_cell_barcodes

# MRO docstring
__MRO__ = """
stage PARSE_ANALYSIS_INPUTS(
    in  csv      parameters,
    in  string   reference_path,
    in  csv      cell_barcodes,
    out csv      cell_barcodes,
    out string[] factorization,
    out int      tsne_perplexity,
    out int      tsne_max_dims,
    out int      tsne_input_pcs,
    out int      tsne_max_iter,
    out int      tsne_stop_lying_iter,
    out int      tsne_mom_switch_iter,
    out float    tsne_theta,
    out int      random_seed,
    out int      max_clusters,
    out float    neighbor_a,
    out float    neighbor_b,
    out int      graphclust_neighbors,
    out int      num_comps,
    out int      num_dr_bcs,
    out int      num_dr_features,
    out int      num_analysis_bcs,
    src py       "stages/reanalysis/parse_analysis_inputs",
) using (
    volatile = strict,
)
"""

type_default = namedtuple("type_default", ["type", "default"])

ANALYSIS_PARAMS = {
    "dim_reduce": type_default(str, ["lsa"]),
    "tsne_perplexity": type_default(int, 30),
    "tsne_max_dims": type_default(int, 2),
    "tsne_input_pcs": type_default(int, 15),
    "tsne_max_iter": type_default(int, 1000),
    "tsne_stop_lying_iter": type_default(int, 250),
    "tsne_mom_switch_iter": type_default(int, 250),
    "tsne_theta": type_default(float, 0.5),
    "random_seed": type_default(int, 0),
    "max_clusters": type_default(int, 10),
    "neighbor_a": type_default(float, -230.),
    "neighbor_b": type_default(float, 120.),
    "graphclust_neighbors": type_default(int, 0),
    "num_comps": type_default(int, 15),
    "num_dr_bcs": type_default(int, None),
    "num_dr_features": type_default(int, None),
    "num_analysis_bcs": type_default(int, None),
}

def main(args, outs):
    parsed = parse_parameters(args.parameters)
    for param in ANALYSIS_PARAMS:
        reparam = "factorization" if param == "dim_reduce" else param
        setattr(outs, reparam, parsed[param] if param in parsed else ANALYSIS_PARAMS[param].default)

    # handle cell barcode format conversion
    if args.cell_barcodes is not None:
        cell_barcodes = get_cell_barcodes(args.cell_barcodes, args.reference_path, with_species=True)
        with open(outs.cell_barcodes, 'w') as outfile:
            for species in cell_barcodes:
                outfile.write(species + ",")
                outfile.write(",".join(cell_barcodes[species]) + "\n")
    else:
        outs.cell_barcodes = None

def parse_parameters(filename):
    if filename is None:
        return {}

    if not os.path.exists(filename):
        martian.exit("Parameters file does not exist: %s" % filename)

    if not os.access(filename, os.R_OK):
        martian.exit("Parameters file is not readable, please check file permissions: %s" % filename)

    params = {}
    with open(filename, 'rU') as f:
        # skip comment lines
        ff = filter(lambda row: not row.startswith('#') and not row.isspace(), f)
        reader = csv.reader(ff)
        for (i, row) in enumerate(reader, start=1):
            if len(row) != 2:
                martian.exit("Row %d is incorrectly formatted (must have exactly 2 columns)" % i)
            name = row[0].strip().lower()
            value = row[1].strip()
            if name not in ANALYSIS_PARAMS:
                martian.exit("Unrecognized parameter: %s" % name)
            if name in params:
                martian.exit("Cannot specify the same parameter twice: %s" % name)
            required_type = ANALYSIS_PARAMS[name].type
            try:
                cast_value = required_type(value)
                params[name] = cast_value
                if name == "dim_reduce":
                    params[name] = [cast_value]
            except ValueError:
                martian.exit("Parameter %s could not be cast to the required type: %s" % (name, required_type))

    return params
