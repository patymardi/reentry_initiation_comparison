#!/usr/bin/env python3
""""

-------------------------------------------------------

  cross_field.py        Program to induce reentry using a
                        cross-field stimulation protocol.
                        A set of N beats (S1) is delivered 
                        followed  a second stimulus (S2). The number 
                        of S1 stimulus is defined by args.prebeats. 
                        The S2 stimulus size is defined wiht vtx file
                        with radius set by args.stimsize

  Ver. 1.0.0

  Created:       Patricia Martínez      (20.01.2026)

  IHU-Liryc  - L'institut Des Maladies Du Rythme Cardiaque

-------------------------------------------------------------------
Input :             mesh        name of the mesh folder e.g patient55

                    par_file    parameter file (.par) that contains the
                                scaling factors for the ionic conductances.
                                Based on Loewe A. et al. 2015 Comp Biomed 9126:439–447

                    stim_file   file containing x,y,z coordinate of points
                                in the left atrium (LA) or right atrium (RA)
                                Generated using getStimPoints.py function
                                e.g. LA_stim_points_patient55.txt

Output:

Usage:              cross_field.py   --mesh patient55 --stim_file LA_stim_points_patient55
"""

import os

EXAMPLE_DESCRIPTIVE_NAME = 'S1S2 inducibility study'
EXAMPLE_AUTHOR = 'Patricia Martinez <patricia.martinez@ihu-liryc.fr>'
EXAMPLE_DIR = os.path.dirname(__file__)

import sys
import shutil
import csv

import vtk
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter as dsa

from datetime import date
from carputils import settings
from carputils import tools
from carputils import model
from carputils.carpio import igb
import numpy as np
import math
import re

def parser():
    # Generate the standard command line parser
    parser = tools.standard_parser()
    # Add arguments

    parser.add_argument('--giL', 
                        type=float, 
                        default=0.3796,
                        help='intracellular longitudinal conductivity to fit CV=0.7 m/s with dx=0.4 mm and dt=20')
    parser.add_argument('--geL', 
                        type=float, 
                        default=1.3635,
                        help='extracellular longitudinal conductivity to fit CV=0.7 m/s with dx=0.4 mm and dt=20')
    parser.add_argument('--cv',
                        type=float, 
                        default=0.7,
                        help='conduction velocity in m/s')
    parser.add_argument('--resolution',
                        default = 400., type = float,
                        help = 'Edge length of triangular slab elements in microns.')
    parser.add_argument('--model',
                        type=str,
                        default='Courtemanche',
                        help='input ionic model')
    parser.add_argument('--mesh',
                        type=str, default='patient55',
                        help='meshname')
    parser.add_argument('--geometry',
                        type=str, default='LA_RA_bilayer_with_fiber_um',
                        help='Geometrical instance in micrometers; \
                        LA_RA_bilayer_with_fiber_um;')
    parser.add_argument('--M_lump',
                        type=int,
                        default='1',
                        help='set 1 for mass lumping, 0 otherwise')
    parser.add_argument('--dt',
                        type=float, default=20.0,
                        help='[microsec]')
    ###############################################
    parser.add_argument('--par_file',
                        type=str,
                        default="parameters.par",
                        help='.par file containing the ionic properties of the atria')
    ######################

    # Single cell prepace
    parser.add_argument('--cell_bcl',
                          type=float,
                          default=1000.0,
                          help='Specify the basic cycle length (ms) to initialize cells')
    parser.add_argument('--numstim',
                          type=int,
                          default=100,
                          help='Specify the number of single cell stimuli with the given bcl')

    #Prepace settings
    parser.add_argument('--prepace_strength',
                        type=float, default=30,
                        help='prepacing injected transmembrane current in uA/cm^2 (2Dcurrent) or uA/cm^3 (3Dcurrent)')
    parser.add_argument('--prepace_bcl',
                        type=float, default=500.0, # should be the same as S1
                        help='initial basic cycle lenght in [ms]')
    parser.add_argument('--prebeats',
                        type=int,
                        default=10,
                        help='Number of beats to prepace the tissue, equivalent to the number of S1 stimuli.'
                             ' S1S2 protocol delivers one additional S1')

    # Stimulus parameters
    parser.add_argument('--stim_current',
                        type=float,
                        default=30,
                        help='Transmembrane current in uA/cm^2 (2Dcurrent)')
    parser.add_argument('--stim_duration',
                        type=float, default=3.0,
                        help='stimulation duration in ms')
    parser.add_argument('--stim_strength',
                        type=float, default=30,
                        help='stimulation transmembrane current in uA/cm^2 (2Dcurrent) or uA/cm^3 (3Dcurrent)')
    parser.add_argument('--stim_size',
                        type=float, default=15000,
                        help='stimulation size in um')
    parser.add_argument('--stim_file',
                        type=str,
                        default= 'LA_stim_points_patient55.txt',
                        help='stimulation point filename')

    # S1S2 Cross-field Protocol
    parser.add_argument('--APD_percent',
                        type=float, default=94.0,
                        help='action potential duration percentage to set as first guess to find the end of the effective refractory period')
    parser.add_argument('--S1',
                        type=int,
                        default=500,
                        help='S1 pacing cycle length') # should be the same as prepace_bcl
    parser.add_argument('--S2_step',
                        type=int,
                        default=0,
                        help='Decrement or increment in S2 coupling in ms, if set to 0 will do S2=ACT+APD')
    parser.add_argument('--tend',
                        type=int,
                        default=10000,
                        help='For how long will the reentry be saved in ms')
    parser.add_argument('--spacedt',
                        type=int,
                        default=1,
                        help='Temporal resolution in ms')
    return parser

def jobID(args):
    today = date.today()
    stim_pt_name = args.stim_file.split('.')[0]

    ID = '../../results/{}/{}_S1S2_{}_stim_size_{}_S1_{}_S2_step_{}_spacedt_{}_tend_{}'.format(args.mesh, today.isoformat(), stim_pt_name, args.stim_size, args.S1,args.S2_step, args.spacedt,args.tend)

    return ID

@tools.carpexample(parser, jobID)
def run(args, job):

    meshdir = '../../data'
    meshname = '{}/{}/{}'.format(meshdir,args.mesh,args.geometry)
    simid = job.ID

    basename = str(job.ID)
    tags = '{}/element_tag.csv'.format(meshdir)
    tag = {}

    #File already provided with all the tag notations  
    with open(tags) as f:
        reader = csv.DictReader(f)
        for row in reader:
            tag[row['name']] = int(row['tag'])

    """
    File with ionic regions definitions
    """
    cmd = tools.carp_cmd('{}/{}'.format(meshdir,args.par_file)) # cmd = tools.carp_cmd('stimulation.par')

    try:
        os.makedirs(job.ID)
    except OSError:
        print ("Creation of the directory %s failed" % job.ID)
    else:
        print ("Successfully created the directory %s " % job.ID)

    geometry = args.geometry.split('_um')[0] # meshname without um, Ids remain the same
    f_slow_conductive = '{}/{}/{}_elems_slow_conductive_U1'.format(meshdir,args.mesh,geometry)
    f_not_conductive = '{}/{}/{}_elems_not_conductive_U1'.format(meshdir,args.mesh,geometry)

    # Region dynamic retagging
    dyn_reg = ['-numtagreg', 2]  # 101 slow_conductive, 103 not_conductive

    # Slow conductive LA
    dyn_reg.extend(tagregopt(0, 'type', 4))  # Type 4 elemfile
    dyn_reg.extend(tagregopt(0, 'elemfile', f_slow_conductive))
    dyn_reg.extend(tagregopt(0, 'tag', 101))

    # Not conductive LA
    dyn_reg.extend(tagregopt(1, 'type', 4))  # # Type 4 elemfile
    dyn_reg.extend(tagregopt(1, 'elemfile', f_not_conductive))
    dyn_reg.extend(tagregopt(1, 'tag', 103))

    cmd += dyn_reg
    cmd_ref = list(cmd)

    stim_file = '{}/{}/{}'.format(meshdir,args.mesh, args.stim_file)
    stim_pts =np.atleast_2d(np.loadtxt(stim_file, skiprows=1, usecols=(0,1,2))) #Skip header in first line
    np.savetxt(job.ID + '/stim_points.txt', stim_pts)

    sinus = ['-stim[0].elec.geomID', tag['sinus_node']]

    steady_state_dir = '{}/{}/init_state/{}'.format(meshdir, args.mesh,args.S1)

    if not os.path.isfile(steady_state_dir+'/single_cell/{}_trace_header.txt'.format(args.model)):
        run_sim=1
    else:
        run_sim=0

    tissue_init = single_cell_initialization(args, job, steady_state_dir + '/single_cell', run_sim)
    startstatef = steady_state_dir + '/vm_last_beat.igb'

    if not os.path.isfile(startstatef):
        print('---------------------------------------------')
        print('Running prepace:')
        print('---------------------------------------------')
        model.prepace.prepace(args, job, cmd, meshname, sinus, steady_state_dir, tissue_init)
    else:
        print('---------------------------------------------')
        print('Prepace file found, skipping prepacing protocol')
        print('---------------------------------------------')

    #Now compute the APD
    igbout = '--output-file=' + job.ID + '/apd_{}.dat'.format(str(args.APD_percent))
    cmd_igb  = [ settings.execs.igbapd,
             '--repol='+str(args.APD_percent),
             igbout,
             startstatef]

    #Run simulation
    job.bash(cmd_igb)

    apdfile = job.ID + '/apd_{}.dat'.format(str(args.APD_percent))
    actfile = steady_state_dir + '/init_acts_vm_act-thresh.dat'

    # Write header for reentries.txt
    write_header = not os.path.exists(simid + '/reentries.txt')

    stim_site = args.stim_file.split('_')[0]

    if stim_site == "LA":  # i < 10: # The first 10 points are in the LA
        tags = [*range(1, 21, 1)]  # Tags for LA
    elif stim_site == "RA":
        tags = [*range(60, 71, 1)]  # for the right atria it has to be done with the epi surface
    else:
        print(
            "Please name the stimulation file starting with LA or RA. Example LA_bilayer_10_points.txt or RA_bilayer_6_points.txt")
        sys.exit()

    tag_str = ','.join([str(i) for i in tags])

    for i in range(len(stim_pts)):

        pacing_coord=stim_pts[i]

        stim_file = simid + '/stim_point_{}.txt'.format(i)

        with open(stim_file, 'w') as f:
            f.write("1\n")
            f.write("{} {} {} 5000\n".format(pacing_coord[0], pacing_coord[1], pacing_coord[2]))

        # Get point ID
        os.system('meshtool query idxlist -msh={} -coord={}'.format(meshname, stim_file))
        index_node = np.loadtxt(stim_file +'.out.txt', dtype=int, skiprows=1, usecols=(0,))
        index_node_str = str(index_node)

        # Remove meshtool helper files
        os.remove(stim_file +'.out.txt')
        os.remove(stim_file)

        node_coord_str = ','.join([str(i) for i in pacing_coord])

        if not os.path.isfile('{}_tag_{}.surf.vtx'.format(meshname, stim_site)):
            # Extract the LA or RA depending on the given LA_stim_points.txt
            os.system("meshtool extract surface -msh={} -surf={}_tag_{} -op={} -coord={}".format(meshname, meshname, stim_site,tag_str, node_coord_str))

        # Extract stimulus vtx with radius defined by args.stimsize
        os.system("meshtool extract surface -msh={} -surf={}_node_{}.vtx -size={} -coord={} -op={}_tag_{}".format(
                meshname, meshname, index_node, float(args.stim_size) , node_coord_str,
                meshname, stim_site))

        vtx_file= meshname + '_node_{}.vtx.surf'.format(index_node)
        # Pacing nodes
        nodes_to_check = np.loadtxt(vtx_file + '.vtx', skiprows=2, dtype=int)


        with open(apdfile) as fp:
            for i, line in enumerate(fp):
                if i == index_node:
                    APD = float(line)
                elif i > index_node:
                    break
        with open(actfile) as fp:
            for i, line in enumerate(fp):
                if i == index_node:
                    ACT = float(line)
                elif i > index_node:
                    break

        start_S2 = float(math.ceil(ACT + APD)) + args.S2_step
        start_S1 = args.prepace_bcl * args.prebeats # S1 should be the same as prepace_bcl = 500 x 10 =  5000.0 ms

        print ('---------------------------------------------')
        print("Node =", index_node)
        print("ACT =", round(ACT,1), "ms")
        print("APD =", round(APD,1), "ms")
        print("Num_S1=", int(args.prebeats))
        print("PCL =", int(args.S1), "ms")
        print("S1 =", int(start_S1), "ms")
        print("S2 =", int(start_S2), "ms")

        if ACT == -1 and APD == -1:
            # Go out of the loop and proceed to the next index_node in stim_pts[i]
            print(f"No arrhythmia initiated at node = {index_node}, because ACT and APD is not defined or null")
            continue   # skips to next iteration of the loop

        tsav_state = args.prepace_bcl * args.prebeats
        startstatef = steady_state_dir + '/state.' + str(tsav_state)

        s1 = ['-num_stim', 2,
                '-stim[0].crct.type', 0, # equivalent to stimtype, 0=transmembrane current
                '-stim[0].elec.geomID', tag['sinus_node'],
                '-stim[0].pulse.strength', args.stim_strength,
                '-stim[0].ptcl.start', start_S1,
                '-stim[0].ptcl.duration', args.stim_duration,
                '-stim[0].ptcl.npls', 1,
                '-stim[0].ptcl.bcl', start_S2 - start_S1]

        s2 = ['-stim[1].crct.type', 0,
                '-stim[1].pulse.strength', args.stim_strength,
                '-stim[1].ptcl.start', start_S2,
                '-stim[1].ptcl.duration', args.stim_duration,
                '-stim[1].ptcl.npls', 1,
                '-stim[1].ptcl.bcl', args.prepace_bcl,
                '-stim[1].elec.vtx_file',        vtx_file]

        #Use this definition if a rectangle is needed
        # p0 = np.array(pacing_coord) - float(args.stim_size) / 2.0
        # p1 = np.array(pacing_coord) + float(args.stim_size) / 2.0

        # s1 = ['-num_stim', 2,
        #         '-stim[0].crct.type', 0,
        #         '-stim[0].elec.geomID', tag['sinus_node'],
        #         '-stim[0].pulse.strength', args.stim_strength,
        #         '-stim[0].ptcl.start', start_S1,
        #         '-stim[0].ptcl.duration', args.stim_duration,
        #         '-stim[0].ptcl.npls', 1,
        #         '-stim[0].ptcl.bcl', start_S2 - start_S1]

        # s2 = ['-stim[1].crct.type', 0,
        #       '-stim[1].pulse.strength', args.stim_strength,
        #       '-stim[1].ptcl.start', start_S2,
        #       '-stim[1].ptcl.duration', args.stim_duration,
        #       '-stim[1].ptcl.npls', 1,
        #       '-stim[1].ptcl.bcl', args.prepace_bcl,
        #       '-stim[1].elec.p0[0]', p0[0],
        #       '-stim[1].elec.p0[1]', p0[1],
        #       '-stim[1].elec.p0[2]', p0[2],
        #       '-stim[1].elec.p1[0]', p1[0],
        #       '-stim[1].elec.p1[1]', p1[1],
        #       '-stim[1].elec.p1[2]', p1[2]]

        #Equivalent old stimulus definition
        # s1 = ['-num_stim', 2,
        #       '-stimulus[0].stimtype', 0,
        #       '-stimulus[0].geometry', tag['sinus_node'],
        #       '-stimulus[0].strength', args.stim_strength,
        #       '-stimulus[0].start', start_S1,
        #       '-stimulus[0].duration', args.stim_duration,
        #       '-stimulus[0].npls', 1,
        #       '-stimulus[0].bcl', start_S2 - start_S1]

        # s2 = ['-stimulus[1].stimtype',         0,
        #     '-stimulus[1].strength',         args.stim_strength,
        #     '-stimulus[1].start',            start_S1, # change
        #     '-stimulus[1].duration',         args.stim_duration,
        #     '-stimulus[1].npls',             1,
        #     '-stimulus[1].vtx_file',        vtx_file,
        #     '-stimulus[1].ctr_def',          1]

        lat = ['-num_LATs', 1,
               '-lats[0].all', 0,
               '-lats[0].measurand', 0,
               '-lats[0].mode', 0,
               '-lats[0].threshold', -50,
               '-lats[0].ID', 1]

        sentinel = ['-t_sentinel', 1,
                    '-t_sentinel_start', start_S2,
                    '-sentinel_ID', 0]

        cmd += lat + sentinel+ s1 + s2

        cmd += ['-simID', simid + '/node_{}'.format(index_node_str),
                '-dt', args.dt,
                '-timedt', 10, # terminal
                '-spacedt', args.spacedt, # change to 1
                '-start_statef', startstatef,
                '-num_tsav', 5,
                '-tsav[0]', tsav_state,
                '-tsav[1]', start_S2 - 1, # before the S2 is delivered
                '-tsav[2]', start_S2 + (args.tend/10),
                '-tsav[3]', start_S2 + (args.tend/2),
                '-tsav[4]', tsav_state + args.tend, #start_S2 + 10000,
                '-mass_lumping', args.M_lump,
                '-tend', start_S2 + args.tend + 0.1,
                '-meshname', meshname]

        # Run simulation
        job.carp(cmd)

        # Find exit files to determine if reentry was terminated
        exit_file = next((f for f in os.listdir(simid + '/node_{}'.format(index_node_str))if f.startswith("exit")),None)

        if exit_file is not None:
            duration = extract_time_from_roe(exit_file)

        else:
            vm_file = simid + '/node_{}'.format(index_node_str) + '/vm.igb'
            # Read header to get duration
            igb_f = igb.IGBFile(f'{vm_file}')
            header = igb_f.header()
            duration = header["t"]  # 10000 + S2 coupling

        with open(simid + '/reentries.txt', 'a') as f:
            if write_header:
                f.write("Node ACT APD Start_S1 Start_S2 duration\n")
                write_header = False
            f.write("{} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n".format(index_node_str, ACT, APD, start_S1, start_S2, duration))

        print ('Done ...')

def tagregopt( reg, field, val ):
    return ['-tagreg['+str(reg)+'].'+field, val ]

def single_cell_initialization(args, job, steady_state_dir,run_sim=0):

    try:
        os.makedirs(steady_state_dir)
    except OSError:
        print ("Creation of the directory %s failed" % steady_state_dir)
    else:
        print ("Successfully created the directory %s " % steady_state_dir)

    duration = args.numstim * args.cell_bcl
    
    '''
    changes from Nikola Fitzen, Jan 2023, due to simulating scenarios with different electrophysiological properties
    Therefore, the single_cell_initialization as part of the prepacing needs to be adjusted according to the selected case of
    remodelling of the properties. Three different states were used by calling the according file name.
    
    al_mk : al = Axel Loewe ,mmk = Martin Krueger --> The parameters used come from the thesis and paper published by them

    _H.par : healthy heart, remodelling only due to heterogeneity
    _M.par : a heart mildy prone to AF
    _S.par : a heart severely prone to AF

    Further explanation can be found in the according files as documentation and of course in the publications
    '''

    g_CaL_reg = [0.45, 0.7515, 0.7515, 0.3015, 0.3015, 0.4770, 0.4770, 0.45, 0.3375]
    g_K1_reg = [2, 2, 2, 2, 2, 2, 2, 2, 1.34]
    blf_g_Kur_reg = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    g_to_reg = [0.35, 0.35, 0.35, 0.5355, 0.5355, 0.238, 0.238, 0.35, 0.2625]
    g_Ks_reg = [2, 2, 2, 2, 2, 2, 2, 2, 3.74]
    maxI_pCa_reg = [1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5]
    maxI_NaCa_reg = [1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6]
    g_Kr_reg = [1, 1, 1, 1.53, 2.44, 1, 1.6, 1.6, 2.4]

    n_regions = len(g_CaL_reg)

    for k in range(n_regions):
        init_file = steady_state_dir+'/init_values_stab_bcl_{}_reg_{}_numstim_{}.sv'.format(args.cell_bcl,k,args.numstim)
        cmd = [settings.execs.BENCH,
        '--imp',                args.model,
        '--imp-par',            'GCaL*{},GK1*{},factorGKur*{},Gto*{},GKs*{},maxIpCa*{},maxINaCa*{},GKr*{}'.format(g_CaL_reg[k],g_K1_reg[k],blf_g_Kur_reg[k],g_to_reg[k],g_Ks_reg[k],maxI_pCa_reg[k],maxI_NaCa_reg[k],g_Kr_reg[k]),
        '--bcl',                args.cell_bcl,
        '--dt-out',             1, # Temporal granularity in ms (if  --dt-out= duration it saves only last value)
        '--stim-curr',          9.5,
        '--stim-dur',           args.stim_duration,
        '--numstim',            args.numstim,
        '--duration',           duration, 
        '--stim-start',         0,
        '--dt',                 args.dt/1000,
        '--fout='               + steady_state_dir + '/{}_numstim_{}_bcl_ms_{}'.format(args.numstim,args.prepace_bcl,k),
        '-S',                   duration,
        '-F',                   init_file,
        '--trace-no',           k]

        if run_sim:
            job.bash(cmd) # Here Generates the JOB.ID Folder

    GCaL_fib = [0.225]      # g_CaL_reg     --> [0.1238]
    GNa_fib = [0.6]         # g_Na_fib      --> [7.8]
    factorGKur_fib = [0.5]  # blf_g_Kur_fib --> [1]
    Gto_fib = [0.35]        # g_to_fib      --> [0.1652]
    GKs_fib = [2]           # g_Ks_fib      --> [0.129]
    maxIpCa_fib = [1.5]     # maxI_pCa_fib  --> [0.275]
    maxINaCa_fib = [1.6]   # maxI_NaCa_fib --> [1600]

    n_regions += len(GCaL_fib)

    for kk in range(len(GCaL_fib)):
        init_file = steady_state_dir + '/init_values_stab_bcl_{}_reg_{}_numstim_{}.sv'.format(args.cell_bcl,
                                                                                              k + 1 + kk,
                                                                                              args.numstim)
        cmd = [settings.execs.BENCH,
               '--imp', args.model,
               '--imp-par', 'GCaL*{},GNa*{},factorGKur*{},Gto*{},GKs*{},maxIpCa*{},maxINaCa*{}'.format(
                GCaL_fib[kk], GNa_fib[kk], factorGKur_fib[kk], Gto_fib[kk], GKs_fib[kk],
                maxIpCa_fib[kk], maxINaCa_fib[kk]),
               '--bcl', args.cell_bcl,
               '--dt-out', 1,
               # Temporal granularity in ms (if  --dt-out= duration it saves only last value)
               '--stim-curr', 9.5,
               '--stim-dur', args.stim_duration,
               '--numstim', args.numstim,
               '--duration', duration,
               '--stim-start', 0,
               '--dt', args.dt / 1000,
               '--fout=' + steady_state_dir + '/{}_numstim_{}_bcl_ms_{}'.format(args.numstim, args.cell_bcl,
                                                                                k + 1 + kk),
               '-S', duration,
               '-F', init_file,
               '--trace-no', k + 1 + kk]
        if run_sim:
            job.bash(cmd)

        # Move trace files to steady_state_dir
    for file in os.listdir(os.getcwd()):
        if file.startswith("Trace") or file.startswith("{}_trace".format(args.model)):
            old_file = os.getcwd() + '/' + file
            new_file = steady_state_dir + '/' + file
            os.rename(old_file, new_file)

    tissue_init = []
    for k in range(n_regions):
        init_file = steady_state_dir + '/init_values_stab_bcl_{}_reg_{}_numstim_{}.sv'.format(args.cell_bcl, k,
                                                                                              args.numstim)
        tissue_init += ['-imp_region[{}].im_sv_init'.format(k), init_file]

    return tissue_init

def prepace(args, job, cmd_ref, meshname, sinus, steady_state_dir, tissue_init):
    """
    Tissue prepacing:

    It consists of a series of pulses at a fixed basic cycle length to let the
    tissue reach a stable limit cycle. An activation time map is computed for
    the last beat. This method will generate an intermediate state to be loaded
    in the other protocols. This is not a protocol to induce arrhythmia.

    Input: parser arguments (args), output directory (job.ID), struct containing imp_regions and gregions (cmd_ref),
    meshname, sinus node location (sinus), prepacing directory (steady_state_dir) and tissue initialization

    Args:
    '--M_lump',
        type=int,
        default='1',
        help='set 1 for mass lumping, 0 otherwise. Mass lumping will speed up the simulation. Use with regular meshes.'
    '--stim_size',
        type=str,
        default='2000.0',
        help='stimulation edge square size in micrometers'
    '--prepace_bcl',
        type=float, default=500.0,
        help='initial basic cycle lenght in [ms]'
    '--prebeats',
        type = int,
        default = 4,
        help='Number of beats to prepace the tissue'

    """
    simid = steady_state_dir

    try:
        os.makedirs(simid)
    except OSError:
        print("Creation of the directory %s failed" % simid)
    else:
        print("Successfully created the directory %s " % simid)

    cmd = list(cmd_ref)
    tsav_state = args.prebeats * args.prepace_bcl  # 4*500= 2000.0
    # Setting the stimulus at the sinus node
    prepace = ['-num_stim', 1,
               '-num_tsav', 1,
               '-tsav[0]', tsav_state,
               '-stimulus[0].stimtype', 0,
               '-stimulus[0].strength', args.prepace_strength,
               '-stimulus[0].duration', args.stim_duration,
               '-stimulus[0].npls', args.prebeats,
               '-stimulus[0].bcl', args.prepace_bcl]
    cmd += tissue_init + prepace + sinus

    cmd += ['-simID', simid,
            '-dt', 20,
            '-spacedt', 10,
            '-mass_lumping', args.M_lump,
            '-timedt', 100,
            '-tend', tsav_state + 0.1,
            '-retagfile', 'retagged.dat', # show regions tagged dynamically
            '-meshname', meshname]
    # Run simulation
    job.carp(cmd)

    old_file = os.path.join(str(simid), "vm.igb")
    new_file = os.path.join(str(simid), "vm_prepace.igb")
    os.rename(old_file, new_file)

    # Last beat
    cmd = list(cmd_ref)

    startstatef = simid + '/state.' + str(tsav_state)
    tsav_state += args.prepace_bcl
    last_beat = ['-num_stim', 1,
                 '-start_statef', startstatef,
                 '-num_tsav', 1,
                 '-tsav[0]', tsav_state,
                 '-stimulus[0].start', args.prebeats * args.prepace_bcl,
                 '-stimulus[0].stimtype', 0,
                 '-stimulus[0].strength', args.prepace_strength,
                 '-stimulus[0].duration', args.stim_duration,
                 '-stimulus[0].npls', 1]

    lat = ['-num_LATs', 1,
           '-lats[0].all', 0,
           '-lats[0].measurand', 0,
           '-lats[0].mode', 0,
           '-lats[0].threshold', -50]

    cmd += last_beat + sinus + lat

    cmd += ['-simID', simid,
            '-dt', 20,
            '-spacedt', 1,
            '-timedt', 100,
            '-mass_lumping', args.M_lump,
            '-tend', tsav_state + 0.1,
            '-meshname', meshname]

    # Run simulation
    job.carp(cmd)

    old_file = os.path.join(str(simid), "vm.igb")
    new_file = os.path.join(str(simid), "vm_last_beat.igb")
    os.rename(old_file, new_file)


def smart_reader(path):
    extension = str(path).split(".")[-1]

    if extension == "vtk":
        data_checker = vtk.vtkDataSetReader()
        data_checker.SetFileName(str(path))
        data_checker.Update()

        if data_checker.IsFilePolyData():
            reader = vtk.vtkPolyDataReader()
        elif data_checker.IsFileUnstructuredGrid():
            reader = vtk.vtkUnstructuredGridReader()

    elif extension == "vtp":
        reader = vtk.vtkXMLPolyDataReader()
    elif extension == "vtu":
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif extension == "obj":
        reader = vtk.vtkOBJReader()
    else:
        print("No polydata or unstructured grid")

    reader.SetFileName(str(path))
    reader.Update()
    output = reader.GetOutput()

    return output


def getEarliestAct(args,job):

    #Read mesh with data
    meshdir = '../../data/meshes/{}/{}'.format(args.mesh,args.scenario)    
    meshname = '{}/bilayer/{}'.format(meshdir, args.geometry)
    mesh = smart_reader(meshname + '.vtk')

    geo_filter = vtk.vtkGeometryFilter()
    geo_filter.SetInputData(mesh)
    geo_filter.Update()
    mesh_surf = geo_filter.GetOutput()

    sinus_coord = mesh_surf.GetPoint(int(args.pacing))

    sinus = ['-stimulus[0].x0', sinus_coord[0],
             '-stimulus[0].xd', args.stim_size,
             '-stimulus[0].y0', sinus_coord[1],
             '-stimulus[0].yd', args.stim_size,
             '-stimulus[0].z0', sinus_coord[2],
             '-stimulus[0].zd', args.stim_size,
             '-stimulus[0].ctr_def', 1]

    return sinus


def extract_time_from_roe(filename):
    match = re.search(r'exit\.\w+\.(\d+\.\d+)\.roe', filename)
    return float(match.group(1)) if match else None

if __name__ == '__main__':
    run()
