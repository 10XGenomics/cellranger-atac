"""
Compute TSS and CTCF scores of an ATAC-seq library.

Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""

from __future__ import division

import json
from metrics import calculate_ctcf_score_and_profile, calculate_tss_score_and_profile

__MRO__ = """
stage REPORT_TSS_CTCF(
    in  csv  tss_relpos,
    in  csv  ctcf_relpos,
    out json summary_metrics,
    src py   "stages/metrics/report_tss_ctcf",
)
"""


def main(args, outs):
    if args.tss_relpos is None or args.ctcf_relpos is None:
        outs.summary_metrics = None
        return

    tss, _, _ = calculate_tss_score_and_profile(relative_positions=args.tss_relpos)
    ctcf, _, _ = calculate_ctcf_score_and_profile(relative_positions=args.ctcf_relpos)

    summary_metrics = {
        'tss_enrichment_score': tss,
        'ctcf_enrichment_score': ctcf,
    }

    with open(outs.summary_metrics, 'w') as outfile:
        json.dump(summary_metrics, outfile)
