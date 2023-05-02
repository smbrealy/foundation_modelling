from ansys.mapdl.core import launch_mapdl
from ansys.mapdl.core.convert import convert_apdl_block
import numpy as np
import pyvista
import matplotlib.pyplot as plt
import os


def wtg_foundation(steel_den=7850, steel_poisson=0.3, steel_YM=200e+9, steel_damping=0.01, # material parameters
                   mon_len = 90, mon_base_r2 = 5, mon_base_r1 = 4.945, mon_top_r2 = 5, # monopile parameters
                   mon_top_r1 = 4.956, mon_embedded_n_nodes=5,mon_submerged_n_nodes=3, mon_air_n_nodes = 2, 
                   tower_len = 135, tower_top_r2 = 3.25, tower_top_r1 = 3.2275, # tower parameters
                   tower_base_r2 = 5, tower_base_r1 = 4.956, tower_n_nodes = 10,
                   nac_mass = 820e+3, nac_height = 2, nac_n_nodes = 1, # nacelle parameters
                   water_switch=1, # water modelling switch
                   fixtures_switch=1, # adds density to the tower section to account for fixtures/fittings
                   embedded_len = 45, 
                   soil_model_switch = 0, soil_k_vec = [100000,100000,100000,100000,100000], soil_eta = 0.01,
                   s_u=50e+3, eps_50=0.008, gamma=8e+3, J=0.5, # p-y soft clay Matlock
                   spring_len = 10, tower_start = 25,
                   wind_load=2e+6,
                   wave_load=0,
                   ):

    submerged_len = mon_len - tower_start - embedded_len
    if submerged_len + embedded_len > mon_len:
        raise Exception("monopile cannot be shorter than submerged length + embedded length")

    mon_n_nodes = mon_embedded_n_nodes + mon_submerged_n_nodes + mon_air_n_nodes # useful quantity

    # nacelle calculations
    vol = nac_mass / steel_den
    nac_rad = np.sqrt((vol/nac_height)/np.pi)
    
    # launch mapdl
    mapdl = launch_mapdl()

    ### Modelling
    mapdl.clear()
    mapdl.prep7()

    ## Material properties
    mapdl.mptemp(1,0) # don't know why I need this temperature table but it doesn't work without it   
    # normal steel
    mapdl.mpdata("EX",1,c1=steel_YM)   
    mapdl.mpdata("PRXY",1,c1=steel_poisson)   
    mapdl.mpdata("DENS",1,c1=steel_den)
    mapdl.mpdata("DMPR",1,c1=steel_damping)
    # submerged steel
    mapdl.mpdata("EX",2,c1=steel_YM)   
    mapdl.mpdata("PRXY",2,c1=steel_poisson)   
    mapdl.mpdata("DENS",2,c1=steel_den + 1030) # added mass method (+ density of sea water)
    mapdl.mpdata("DMPR",2,c1=steel_damping)
    # steel + fixture/fittings
    mapdl.mpdata("EX",3,c1=steel_YM)   
    mapdl.mpdata("PRXY",3,c1=steel_poisson)   
    mapdl.mpdata("DENS",3,c1=steel_den + 650) # nominally 8500kg/m3 ~ + 8%
    mapdl.mpdata("DMPR",3,c1=steel_damping)


    # Create element types
    mapdl.et(1,"BEAM188") # beam element
    mapdl.et(2,"COMBIN14", 0, 0, 0) # linear spring element (linear solution, option 3, 3d longitudinal spring damper)
    mapdl.et(3,"COMBIN39", 0, 0, 0, 1, 0, 1) # non-linear spring element

    ## create keypoints for geometry  

    keypoints = [0, embedded_len, embedded_len + submerged_len, mon_len, 
                 mon_len + tower_len, mon_len + tower_len + nac_height]

    kp_num = 1
    for i in range(len(keypoints)):
        mapdl.k(kp_num, keypoints[i], 0, 0)
        kp_num += 1

    # connect keypoints with lines
    mapdl.lstr(1, 2) # monopile line
    mapdl.lstr(2, 3) # tower line
    mapdl.lstr(3, 4) # nacelle line
    mapdl.lstr(4, 5) # tower line
    mapdl.lstr(5, 6) # nacelle line

    ## cross-sections
    # monopile base
    mapdl.sectype(1, "BEAM", "CTUBE", "monopile", 0)  
    mapdl.secoffset("CENT") 
    mapdl.secdata(mon_base_r1,mon_base_r2,0,0,0,0,0,0,0,0,0,0)  

    # monopile top
    mapdl.sectype(2, "BEAM", "CTUBE", "monopile", 0)  
    mapdl.secoffset("CENT") 
    mapdl.secdata(mon_top_r1,mon_top_r2,0,0,0,0,0,0,0,0,0,0)  

    # tower base
    mapdl.sectype(3, "BEAM", "CTUBE", "tower", 0)  
    mapdl.secoffset("CENT") 
    mapdl.secdata(tower_base_r1,tower_base_r2,0,0,0,0,0,0,0,0,0,0)  

    # tower top
    mapdl.sectype(4, "BEAM", "CTUBE", "tower", 0)  
    mapdl.secoffset("CENT") 
    mapdl.secdata(tower_top_r1,tower_top_r2,0,0,0,0,0,0,0,0,0,0)  

    # monopile taper
    mapdl.sectype(5,"TAPER", name="monopile_taper") 
    mapdl.secdata(1, 0) 
    mapdl.secdata(2, mon_len) 

    # tower taper
    mapdl.sectype(6,"TAPER", name="tower_taper") 
    mapdl.secdata(3, mon_len) 
    mapdl.secdata(4, mon_len + tower_len) 

    # nacelle
    mapdl.sectype(7, "BEAM", "CSOLID", "nacelle", 0)  
    mapdl.secoffset("CENT") 
    mapdl.secdata(nac_rad,0,0,0,0,0,0,0,0,0,0,0)  

    ### ------------- meshing ---------------
    # create keypoint for cross-section orientation
    mapdl.k(7, 0, 20, 0)

    
    ## Meshing line 1 (embedded monopile)
    # set appropriate material
    mapdl.mat(1)
    # select first line
    mapdl.lsel("", "", "", 1)
    # select correct cross section and orientation
    mapdl.latt(1,1,1,7,"",5)
    # Set number of nodes
    mapdl.lesize("ALL", "", "", mon_embedded_n_nodes, "", "", "", "", 1)
    # Mesh line 1
    mapdl.lmesh(1)

    # submerged monopile
    if water_switch == 1:
        mapdl.mat(2)
    mapdl.lsel("", "", "", 2)
    mapdl.latt(1,1,1,7,"",5)
    mapdl.lesize("ALL", "", "", mon_submerged_n_nodes, "", "", "", "", 1)
    mapdl.lmesh(2)

    # monopile in air
    mapdl.mat(1)
    mapdl.lsel("", "", "", 3)
    mapdl.latt(1,1,1,7,"",5)
    mapdl.lesize("ALL", "", "", mon_air_n_nodes, "", "", "", "", 1)
    mapdl.lmesh(3)

    # tower
    if fixtures_switch == 1:
        mapdl.mat(3)
    mapdl.lsel("", "", "", 4)
    mapdl.latt(1,1,1,7,"",6)
    mapdl.lesize("ALL", "", "", tower_n_nodes, "", "", "", "", 1)
    mapdl.lmesh(4)

    # nacelle
    mapdl.mat(2)
    mapdl.lsel("", "", "", 5)
    mapdl.latt(1,1,1,7,"",7)
    mapdl.lesize("ALL", "", "", nac_n_nodes, "", "", "", "", 1)
    mapdl.lmesh(5)

    # Workaround for material selection not working during meshing
    if water_switch:
        mapdl.esel("","","",mon_embedded_n_nodes + 1,mon_embedded_n_nodes + mon_submerged_n_nodes)
        mapdl.mpchg(2,"ALL")
    if fixtures_switch:
        mapdl.esel("","","",mon_n_nodes + 1, mon_n_nodes + tower_n_nodes)
        mapdl.mpchg(3,"ALL") 
    mapdl.esel("ALL")

    # # Springs
    n_springs = mon_embedded_n_nodes + 1 # get appropriate num nodes below soil depth
    spring_gap = embedded_len/mon_embedded_n_nodes # gap between springs (i.e. gap between nodes)
    spring_depths = []
    for i in range(n_springs):
        spring_depths.append(0 + i*spring_gap) # shallowest to deepest

    ### ---------- Soil Modelling -------------

    # Initialise natural frequency array
    pop_size = len(soil_k_vec)

    #nat_freqs = np.zeros((6, pop_size)) # first 6 natural frequencies
    nat_freqs = np.zeros((pop_size,1))

    for p in range(pop_size):
        ## Set up constants

        ## put APDL back in the correct mode
        mapdl.prep7()
        
        # variables to change for population member
        soil_k = soil_k_vec[p]
        
        # encastre
        if soil_model_switch == 0:
            mapdl.dk(1,"ALL",0) # constrains base of tower
        # linear springs
        elif soil_model_switch == 1: 
            mapdl.type(2) # Change to COMBIN14
            mapdl.r(1, soil_k, soil_eta)
        # non-linear springs
        elif soil_model_switch == 2: # Matlock soft clay
            mapdl.type(3) # change to COMBIN39
            y = np.linspace(0,0.5,30) # displacements 
            D = mon_base_r2 * 2 # monopile diameter 
            y_50 = 2.5 * eps_50 * D # normalising parameter
            z = np.array(spring_depths)
            z = np.reshape(z, (len(z), 1)) # depths below mudline (at spring depths) n_springs x 1
            Zr = (6 * s_u * D) / ((gamma * D) + (J * s_u)) # transition depth
            
            # ultimate lateral soil resistance per unit length, for each depth
            Pu1 = (3 + (gamma * z / s_u) + (J * z / D)) * s_u * D # lateral bearing capacity factor (n_springs x 1)
            Pu2 = 9 * s_u * D
            Pu1[Pu1>Pu2] = Pu2
            p_ult = Pu1 
            
            # soil resistance at rest
            p_rest = np.zeros(len(spring_depths)).reshape(-1,1)
            p_rest[z >= Zr] = 0.72 * p_ult[z >= Zr]
            p_rest[z < Zr] = 0.72 * p_ult[z < Zr] * (z[z < Zr] / Zr)

            #p = np.zeros((len(y), len(p_ult)))
            p = np.zeros((len(spring_depths), len(y)))
            for i in range(len(y)): # loop through displacements, y
                if y[i] <= (3 * y_50):
                    p[:,i] = np.squeeze(0.5 * p_ult * ((y[i] / y_50) ** (1/3)))
                elif y[i] > (3 * y_50) and y[i] < (15 * y_50):
                    p[:,i] = np.squeeze(p_rest + (((15 * y_50) - y[i]) / ((12 * y_50) * 0.72 * p_ult - p_rest)))
                else: # y[i] >= 15 * y_50
                    p[:,i] = np.squeeze(p_rest)

            # now need to multiply p vector by length of element we want the spring to apply to since it is force/unit length
            # i.e. each p x gap between springs with half gap at base and surface of monopile
            p[1:-1,:] *= spring_gap # intermediate springs
            p[0,:] *= spring_gap/2 # monopile surface spring (underestimate) #TODO could be improved
            p[-1,:] *= spring_gap/2 # monopile base spring (overestimate) #TODO could be improved

            # set up non-linear curve constants
            for i in range(n_springs):
                args = []
                args.append(i+1)
                for j in range(len(y)):
                    args.append(y[j])
                    args.append(p[i,j])
    #            print(args)

                if len(args) <= 7:
                    mapdl.r(*args)
                    args = []
                else:
                    mapdl.r(*args[0:7])
                    args = args[7:]
                while len(args) >= 6:
                    mapdl.rmore(*args[0:6])
                    args = args[6:]
                if len(args) > 0:
                    mapdl.rmore(*args)
                    args = []

        ## Create springs if needed
        # set up query
        if p == 0: # only need to create springs for the first population member
            q = mapdl.queries
            if soil_model_switch != 0: # i.e. if modelling soil (not encastre), create springs
                ## lateral stiffnesses
                # y-direction
                for i in range(n_springs):
                    # set correct stiffness
                    if soil_model_switch == 1:
                        mapdl.real(1) # linear model
                    elif soil_model_switch == 2:  
                        mapdl.real(i+1) # non linear model
                    else:
                        raise Exception("soil model not supported")
                    # create node
                    mapdl.n("", embedded_len - i*spring_gap, spring_len, 0)
                    # create spring
                    mapdl.e(q.node(embedded_len - i*spring_gap, 0, 0), q.node(embedded_len - i*spring_gap, spring_len, 0))

                # z-direction
                for i in range(n_springs):
                    # set correct stiffness
                    if soil_model_switch == 1:
                        mapdl.real(1) # linear model
                    elif soil_model_switch == 2:  
                        mapdl.real(i+1) # non linear model
                    else:
                        raise Exception("soil model not supported")
                    # create node
                    mapdl.n("", embedded_len - i*spring_gap, 0, spring_len)
                    # create spring
                    mapdl.e(q.node(embedded_len - i*spring_gap, 0, 0), q.node(embedded_len - i*spring_gap, 0, spring_len))

                # x-direction
                # vertical stiffness
                # create node
                mapdl.r(7, soil_k*100000, soil_eta) # assumed linear spring for base
                mapdl.real(n_springs + 1) # change to correct constant set
                mapdl.type(2) # Change to COMBIN14
                #mapdl.dk(1,"UX",0)
                mapdl.n("", -spring_len, 0, 0)
                # create spring
                mapdl.e(q.node(0, 0, 0), q.node(-spring_len, 0, 0))

                # Add node constrains to free end of springs
                #mapdl.d(q.node(0, spring_len, 0),"ALL",0,"",q.node(spring_depths[0], 0, spring_len))
                mapdl.d(q.node(embedded_len, spring_len, 0),"ALL",0,"",q.node(-spring_len, 0, 0)) # first spring to last spring

            ## Apply force for wind loading
            mapdl.fk(6,"FY",wind_load)

        ### Analysis ###

        # ## Nodes of interest
        # # total number of nodes of interest
        # node_count = mon_embedded_n_nodes + mon_submerged_n_nodes + mon_air_n_nodes + tower_n_nodes + nac_n_nodes 

        # # get lengths of segments of monopile/tower/nacelle
        # mon_embedded_node_len = embedded_len / mon_embedded_n_nodes
        # mon_submerged_node_len = submerged_len / mon_submerged_n_nodes
        # mon_air_node_len = (mon_len - embedded_len - submerged_len) / mon_air_n_nodes
        # tower_node_len = tower_len / tower_n_nodes
        # nac_node_len = nac_height / nac_n_nodes

        # mon_n_nodes = mon_embedded_n_nodes + mon_submerged_n_nodes + mon_air_n_nodes
        # # get x coordinates of nodes of interest
        # nodal_x_coord = [0] * (node_count + 1)
        # for i in range(node_count):
        #     if i < mon_embedded_n_nodes:
        #         nodal_x_coord[i+1] = nodal_x_coord[i] + mon_embedded_node_len
        #     elif i < mon_embedded_n_nodes + mon_submerged_n_nodes:
        #         nodal_x_coord[i+1] = nodal_x_coord[i] + mon_submerged_node_len
        #     elif i < mon_embedded_n_nodes + mon_submerged_n_nodes + mon_air_n_nodes:
        #         nodal_x_coord[i+1] = nodal_x_coord[i] + mon_air_node_len
        #     elif i < mon_n_nodes + tower_n_nodes:
        #         nodal_x_coord[i+1] = nodal_x_coord[i] + tower_node_len
        #     elif i < node_count:
        #         nodal_x_coord[i+1] = nodal_x_coord[i] + nac_node_len

        ## modal analysis
        mapdl.run("/SOL")
        mapdl.modal_analysis(nmode=20, freqb=0, freqe=10, mxpand=20)

        mapdl.finish()
        mapdl.post1()

        # natural frequencies
        #test = mapdl.post_processing.frequency_values[:6]

        ## original with 6 nat frequencies
        #nat_freqs[:,p] = mapdl.post_processing.frequency_values[:6]
        nat_freqs[p] = mapdl.post_processing.frequency_values[1]

    
    return nat_freqs




