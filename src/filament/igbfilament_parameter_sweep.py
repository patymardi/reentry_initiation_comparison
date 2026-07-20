import pandas as pd
import argparse
from run_igbfilament import create_folder_lookup
import os

def parser():
    parser = argparse.ArgumentParser(description="Generate summary table after running sensitivity analysis for igbfilament")
    parser.add_argument("--for_analysis", default='/mnt/smb/patricia/Bitbucket/04_inducibility_comparison/results/patient12/for_analysis.csv', help="Path to list of selected cases")
    parser.add_argument("--aux_mesh", default="igbfilament", help="Auxiliary mesh output directory")

    return parser

def read_pts_file(pts_file,patient,chamber,node,step,threshold,dt):

    """
    Read a phase_sing*.pts_t file.

    Returns
    -------
    DataFrame
        One row per detected phase singularity.
    """

    rows = []

    with open(pts_file, "r") as f:

        # Skip header
        f.readline()      # ## Format x y z # singularity type (0=PS || 1=FilSeg) Global element index
        f.readline()      # number of nodes

        while True:

            line = f.readline()

            if not line:
                break

            line = line.strip()

            # Skip empty lines
            if line == "":
                continue

            # Read only time markers
            if not line.startswith("#"):
                continue

            time = float(line[1:])

            n_ps = int(f.readline().strip())

            for ps_id in range(n_ps):
                coords = f.readline().split()

                rows.append({
                    "patient": patient,
                    "chamber": chamber,
                    "node": node,
                    "step": step,
                    "threshold": threshold,
                    "dt": dt,
                    "time": time,
                    "ps_id": ps_id,
                    "x": float(coords[0]),
                    "y": float(coords[1]),
                    "z": float(coords[2]),
                    "type": int(coords[-2]),  # 0 = PS, 1 = FilSeg
                    "element": int(coords[-1]),
                })

        return pd.DataFrame(rows)


def main ():

    args = parser().parse_args()
    base_folder = os.path.dirname(args.for_analysis)

    selected = pd.read_csv(args.for_analysis)

    thresholds = [-50, -20, -10, 0]
    dts = [2, 4, 6, 8, 10]

    folder_lookup = create_folder_lookup(base_folder)

    all_ps = []

    for _, row in selected.iterrows():

        node = row["Node"]
        chamber = row["chamber"]
        step = row["step"]
        patient = row["patient"]

        # patient_folder  e.g. '2026-06-02_S1S2_RA_stim_points_patient51_stim_size_15000.0_S1_500_S2_step_-24_spacedt_1_tend_10000'
        patient_folder = folder_lookup[(patient, chamber, step)]

        # Define output directory for PSoutput
        ps_output_dir = os.path.join(base_folder,patient_folder,f"node_{node}", args.aux_mesh)

        for threshold in thresholds:
            for dt in dts:
                pts_file = os.path.join (ps_output_dir,f'phase_sing_t_{threshold}_d_{dt}.pts_t')
                if not os.path.exists(pts_file):
                    print(f"Missing: {pts_file}")
                    continue
                df = read_pts_file(pts_file,patient,chamber,node,step,threshold,dt)

                all_ps.append(df)

    all_ps = pd.concat(all_ps, ignore_index=True)

    all_ps.to_csv(os.path.join(base_folder, "phase_singularities.csv"),index=False)

    frame_summary = (all_ps.groupby(["patient","chamber","node","step","threshold","dt","time"]).size().reset_index(name="n_ps"))

    summary = (frame_summary.groupby(["patient","chamber","node","step","threshold","dt"])
        .agg(total_ps=("n_ps", "sum"),
            mean_ps=("n_ps", "mean"),
            median_ps=("n_ps", "median"),
            max_ps=("n_ps", "max"),
            std_ps=("n_ps", "std"),
            frames_with_ps=("n_ps", "count")).reset_index())

    summary.to_csv(os.path.join(base_folder, "phase_singularities_summary.csv"),index=False)
    print(f'Done...')


if __name__ == "__main__":
    main()
