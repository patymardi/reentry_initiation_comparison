#!/usr/bin/env python3
""""

-------------------------------------------------------

  run_igbfilament.py    Python script to run igbfilament. 
                        Of course igbfilament can also be run
                        using the command line. This script is
                        for python lovers only.  

  Ver. 1.0.0

  Created:       Patricia Martínez      (20.01.2026)

  IHU-Liryc  - L'institut Des Maladies Du Rythme Cardiaque

-------------------------------------------------------------------
Input :     path_to_igbfilament     Directory where igbfilament binary
                                    is located. Usually inside:
                                    path_to_openCARP/_build/bin/igbfilament
                                    You can find more information on:
                                    https://git.opencarp.org/openCARP/openCARP/-/tree/master/tools/igbutils

            igb_file                Directory where IGB file is stored

            threshold               Activation threshold (mV)           

            dt                      Time embedding lag (ms)

            aux_mesh                Output auxiliary mesh directory

Output:

Usage:              

    run_igbfilament.py --igb_file ../results/patient51/2026-06-02_cross_field_patient51_RA_stim_size_15000.0_S1_500_S2_step_-24_spacedt_1_tend_10000/node_304000 
                       --threhold -20 --dt 6 --aux_mesh igbfilament

"""



import subprocess
import argparse
import os
import pyvista as pv
import numpy as np
import pandas as pd
from carputils.carpio import igb
import shutil


def parser():
    parser = argparse.ArgumentParser(description="Run igbfilament with specified parameters")
    parser.add_argument("--path_to_igbfilament", default= 'path_to_openCARP/_build/bin/igbfilament', help="Path to igbfilament binary")
    parser.add_argument("--igb_file", help="Path to IGB data file")
    parser.add_argument("--threshold", type=int, default=-20, help="Threshold value")
    parser.add_argument("--dt", type=int, default=6, help="Time step (dt)")
    parser.add_argument("--aux_mesh", default="igbfilament", help="Auxiliary mesh output directory")
    parser.add_argument("--meshdir", default="../data", help="Directory where bilayer meshes are stored")
    parser.add_argument("--tmp_folder", default=None, help="Optional local directory to copy igb files and run igbhead -j")

    return parser

def main():

    args = parser().parse_args()

    # Split the path
    parts = args.igb_file.strip("/").split("/")

    # Extract components from known path structure
    patient = parts[-3]  # e.g. patient51

    meshname= os.path.join(args.meshdir,patient, 'LA_RA_bilayer_with_fiber_um')

    node = parts[-1]  # e.g. node_304000
    patient_folder = parts[-2]  # e.g. '2026-06-02_cross_field_patient51_RA_stim_size_15000.0_S1_500_S2_step_-24_spacedt_1_tend_10000'
    chamber = patient_folder.split('_')[2]  # RA
    step = int(patient_folder.split('_')[-5])  # -24
    ps_output_dir = os.path.join(args.igb_file, args.aux_mesh)

    # Create analysis folder
    try:
        os.makedirs(ps_output_dir)
    except OSError:
        print("Creation of the directory %s failed" % ps_output_dir)
    else:
        print("Successfully created the directory %s " % ps_output_dir)
   
    jive_igb(args.igb_file, args.tmp_folder)

    igbfilament_cmd(args, meshname, args.igb_file, ps_output_dir, args.threshold, args.dt)

    print('Done ...')


def jive_igb(igb_file, tmp_folder=None):

    # Read original header
    igb_f = igb.IGBFile(igb_file)
    header = igb_f.header()
    duration = header["dim_t"]

    print(f"Original: duration={header['dim_t']}, inc_t={header['inc_t']}")

    # Decide which file to jive
    if tmp_folder:
        os.makedirs(tmp_folder, exist_ok=True)
        jive_file = os.path.join(tmp_folder, "vm.igb")

        print(f"Copying\n  {igb_file}\n-> {jive_file}")
        shutil.copy2(igb_file, jive_file)
    else:
        print("Running jive on the original file")
        jive_file = igb_file


    # Run jive locally
    jive = subprocess.run(
        ["igbhead", "-j", jive_file],
        capture_output=True,
        text=True
    )

    if jive.returncode != 0:
        raise RuntimeError(
            f"igbhead failed:\nstdout:\n{jive.stdout}\nstderr:\n{jive.stderr}"
        )

    # Read updated header
    igb_f = igb.IGBFile(jive_file)
    header = igb_f.header()
    duration_after = header["dim_t"]

    print(f"After jive: duration={duration_after}, inc_t={inc_t_after}")

    # Only overwrite if something changed
    if duration_after != duration:
        print("Header changed. Copying corrected file back.")
        shutil.copy2(tmp_igb, igb_file)
    else:
        print("Header already correct. No copy needed.")


def igbfilament_cmd(args,meshname,igb_file,ps_output_dir,threshold, dt):

    # Run igbfilament
    igbfilament_cmd = [
        args.path_to_igbfilament,
        meshname ,
        igb_file,
        "-t", str(threshold),
        "-d", str(dt),
        "-a", os.path.join(ps_output_dir,f"phase_sing_t_{threshold}_d_{dt}")
    ]

    print("Running igbfilament command:", " ".join(igbfilament_cmd))
    subprocess.run(igbfilament_cmd, check=True)


if __name__ == "__main__":
    main()

