"""Synthetic example: result HDF5, full-event CSV, HTML, and optional input H5AD.

Run from the repository root:
    python examples/export_explain_report.py --output examples/report_demo_output --include-input

The generated expression and small example priors demonstrate the workflow;
they are not biological validation data. No downloads are required.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from cellmesh import CellMeshResult, read_example_data, run_cell_mesh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("cellmesh-report-demo"))
    parser.add_argument("--include-input", action="store_true",
                        help="Retain and save the selected cell-by-gene matrix as a companion .h5ad")
    parser.add_argument("--selection-demo", action="store_true",
                        help="Select A/C -> B, export the selected CSV, focused report and dotplot")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    adata = read_example_data("tiny")
    adata.obs["sample"] = np.resize(["S1", "S2", "S3"], adata.n_obs)
    enzyme = pd.DataFrame([
        {"metabolite": "Synthetic M", "hmdb_id": "HMDB0000001", "reaction": reaction,
         "gene": gene, "role": role, "source": "synthetic example"}
        for reaction, gene, role in [
            ("prod", "Gene1;Gene2", "production"),
            ("partial", "Gene3;UNMEASURED", "production"),
            ("consume", "Gene4", "degradation"),
        ]
    ])
    interaction = pd.DataFrame([
        {"metabolite": "Synthetic M", "hmdb_id": "HMDB0000001", "sensor_gene": gene,
         "sensor_type": "Cell surface receptor", "source": "synthetic example"}
        for gene in ["Gene8", "UNMEASURED_SENSOR"]
    ])
    options = dict(sample_key="sample", sample_mode="sample_aware", min_cells=3,
                   n_perms=19, random_state=42)
    result = run_cell_mesh(adata, enzyme, interaction, retain_input=args.include_input, **options)
    selector = dict(hmdb_id="HMDB0000001", sensor_gene="Gene8", sender="A", receiver="B")
    selected = result.select_events(metabolite_sensor_pairs=[("HMDB0000001", "Gene8")],
                                    cell_pairs=[("A", "B")], source="all")
    if args.selection_demo:
        focus = result.filter_events(source="all", sender_labels=["A", "C"], receiver_labels=["B"],
                                     cell_pairs=[("A", "B"), ("C", "B")])
        files = result.save(args.output / "analysis", report=True, filtered_events=focus, selected_events=selected,
                            include_input=args.include_input, include_all_events=True)
        pd.testing.assert_frame_equal(pd.read_csv(files["filtered_csv"], float_precision="round_trip").drop(columns="positive_status"),
                                      focus.events, check_dtype=False, check_exact=True)
        from cellmesh import plot_event_dotplot
        import matplotlib.pyplot as plt
        plot = plot_event_dotplot(focus)
        files["selected_plot"] = args.output / "selected_events.png"
        plot["fig"].savefig(files["selected_plot"], dpi=160, bbox_inches="tight")
        plt.close(plot["fig"])
        print(f"Selection: {focus.summary['n_selected']} / {focus.summary['n_total']} events; "
              f"{plot['selection']['n_plotted']} plotted; 1 explicit molecular-pair report")
    else:
        files = result.save(args.output / "analysis", report=True, selected_events=selected,
                            include_input=args.include_input, include_all_events=True)
    restored = CellMeshResult.load(files["result"])
    pd.testing.assert_frame_equal(result.events, restored.events, check_exact=True)
    csv_events = pd.read_csv(files["events_csv"], float_precision="round_trip")
    pd.testing.assert_frame_equal(result.events, csv_events, check_dtype=False, check_exact=True)
    if args.include_input:
        archived = restored.load_input()
        replay = run_cell_mesh(archived, restored.priors["enzyme"], restored.priors["interaction"],
                               **options)
        pd.testing.assert_frame_equal(result.events, replay.events, check_exact=True)
        print(f"Archived input: {archived.n_obs} cells × {archived.n_vars} genes; exact replay passed")
    detail = restored.explain_event(**selector)
    print(detail["sample_events"][["sample", "cell_mesh_score", "evaluation_status", "reasons"]].to_string(index=False))
    print(restored.coverage_report()["summary"].to_string(index=False))
    for name, path in files.items():
        print(f"{name}: {path.resolve()} ({path.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
