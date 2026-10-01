"""Generate a synthetic 216-event report to demonstrate offline pagination.

Run from the repository root:
    python examples/browse_event_report.py --output examples/report_browser_demo_output

Creates selected HTML plus its selected CSV. Scores are computed from synthetic input;
permutations are disabled, so FDR is explicitly NA. This is a UI demonstration,
not biological validation. The standard result save() interface is unchanged.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from cellmesh import read_example_data, run_cell_mesh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("cellmesh-browser-demo"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    adata = read_example_data("tiny")
    adata.obs["sample"] = np.resize(["S1", "S2", "S3"], adata.n_obs)
    names = [(f"HMDB{i:07d}", f"Synthetic M {i:02d}") for i in range(1, 25)]
    enzyme = pd.DataFrame([
        dict(hmdb_id=identifier, metabolite=name, reaction="prod", gene="Gene1;Gene2",
             role="production", source="synthetic browser example") for identifier, name in names
    ])
    interaction = pd.DataFrame([
        dict(hmdb_id=identifier, metabolite=name, sensor_gene="Gene8", sensor_type="Receptor",
             source="synthetic browser example") for identifier, name in names
    ])
    result = run_cell_mesh(adata, enzyme, interaction, sample_key="sample", sample_mode="sample_aware",
                           min_cells=3, n_perms=0, random_state=42)
    focus = result.filter_events(source="all")
    pairs = list(focus.events[["hmdb_id", "sensor_gene"]].drop_duplicates().head(3).itertuples(index=False, name=None))
    selected = result.select_events(metabolite_sensor_pairs=pairs, source="all")
    pages = result.save(args.output / "analysis", filtered_events=focus, selected_events=selected)
    html = pages["navigation"]
    csv = pages["filtered_csv"]
    pd.testing.assert_frame_equal(pd.read_csv(csv, float_precision="round_trip").drop(columns="positive_status"), focus.events,
                                  check_dtype=False, check_exact=True)
    assert len(focus.events) == 216
    assert sum(key.startswith('pair_report:') for key in pages) == 3
    print(f"{len(focus.events)} interest events, 50 per page, 3 molecular-pair reports; FDR unavailable")
    print(html.resolve())
    print(csv.resolve())


if __name__ == "__main__":
    main()
