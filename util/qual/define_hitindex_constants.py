"""Define HITindex constants used by qual classification."""

HITINDEX_SPLIT_OVERLAP_BP = 10
HITINDEX_MIN_JUNCTION_READS = 2
TE_BOUNDARY_MIN_SIDE_READS = 2
TE_JUNCTION_INTERNAL_GENE_PREFIX = "__TEJUNC__"
HIT_IDENTITY_PARAMETERS = {
    "HITterminal": "1.0",
    "HIThybrid": "0.3",
    "HITpval": "1",
    "HIT_CI": "none",
    "prob_med": "0.5",
    "prob_high": "0.8",
}
