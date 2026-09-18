"""Recompute unique-HMDB exporter coverage from explicit or selected default files.

Usage: python -m cellmesh.prior_coverage [--enzyme FILE] [--interaction FILE]
No expression filtering or score/weight calibration is performed.
"""
from __future__ import annotations

import argparse
import json

from .database import load_cell_mesh_database


def summarize_prior_coverage(enzyme_file=None, interaction_file=None):
    enzyme, sensor = load_cell_mesh_database(enzyme_file, interaction_file)
    enzyme_ids = set(enzyme['hmdb_id'].dropna())
    export_ids = set(enzyme.loc[enzyme['role'].eq('export'), 'hmdb_id'].dropna())
    shared_ids = enzyme_ids & set(sensor['hmdb_id'].dropna())

    def coverage(universe):
        numerator, denominator = len(export_ids & universe), len(universe)
        return {'exporter_hmdb_count': numerator, 'total_hmdb_count': denominator,
                'percent': 100.0 * numerator / denominator if denominator else None}

    return {
        'definition': 'Unique normalized HMDB IDs; export role present; no expression filtering',
        'prior_inputs': {'enzyme': enzyme.attrs['input_provenance'],
                         'interaction': sensor.attrs['input_provenance']},
        'enzyme_all': coverage(enzyme_ids),
        'shared_with_interaction': coverage(shared_ids),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--enzyme', help='Enzyme CSV; omitted selects the default file')
    parser.add_argument('--interaction', help='Interaction CSV; omitted selects the default file')
    args = parser.parse_args()
    print(json.dumps(summarize_prior_coverage(args.enzyme, args.interaction),
                     ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
