from ansys.mapdl.core import launch_mapdl
from ansys.mapdl.core.convert import convert_apdl_block
import numpy as np
import pyvista
import matplotlib.pyplot as plt

def wtg_foundation(steel_den=7850, steel_poisson=0.3, steel_YM=200e+9, steel_damping=0.01, # material parameters
                   mon_len = 90, mon_base_r2 = 5, mon_base_r1 = 4.945, mon_top_r2 = 5, # monopile parameters
                   mon_top_r1 = 4.956, mon_embedded_n_nodes=5,mon_submerged_n_nodes=3, mon_air_n_nodes = 2, 
                   tower_len = 135, tower_top_r2 = 3.25, tower_top_r1 = 3.2275, # tower parameters
                   tower_base_r2 = 5, tower_base_r1 = 4.956, tower_n_nodes = 10,
                   nac_mass = 820e+3, nac_height = 2, nac_n_nodes = 1, # nacelle parameters
                   water_switch=1, # water modelling switch
                   fixtures_switch=1, # adds density to the tower section to acocunt for fixtures/fittings
                   soil_embedded_depth = 45, soil_poisson = 0.4, soil_shear_mod = 50e+6, # soil parameters
                   soil_model_switch = 0, soil_k = 10000000, soil_eta = 0.01,
                   s_u=50e+3, eps_50=0.008, gamma=8e+3, J=0.5, # p-y soft clay Matlock
                   soil_kpy=34e+3, soil_C1=3, soil_C2=3.5, soil_C3=2.7, soil_gamma = 17, soil_A=0.9, # p-y API sand
                   spring_len = 10, submerged_len = 20,
                   scour_n_nodes = 0 # scour parameters
                   ):
    
    if submerged_len + soil_embedded_depth > mon_len:
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

    keypoints = [0, soil_embedded_depth, soil_embedded_depth + submerged_len, mon_len, 
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
    spring_gap = soil_embedded_depth/mon_embedded_n_nodes # gap between springs (i.e. gap between nodes)
    spring_depths = []
    for i in range(n_springs):
        spring_depths.append(soil_embedded_depth - i*spring_gap) # deepest to shallowest

    ### ---------- Soil Modelling -------------
    ## Set up constants
    # encastre
    if soil_model_switch == 0:
        mapdl.dk(1,"ALL",0) # constrains base of tower
    # linear springs
    elif soil_model_switch == 1: 
        mapdl.type(2) # Change to COMBIN14
        mapdl.r(1, soil_k, soil_eta)
    # non-linear springs
    elif soil_model_switch == 2: # Matlock soft clay - note this is the static model (not cyclic loading)
        mapdl.type(3) # change to COMBIN39
        y = np.linspace(0,1,6) # displacements
        D = mon_base_r2 * 2
        y_50 = 2.5 * eps_50 * D # normalising parameter
        z = np.array(spring_depths)
        z = np.reshape(z, (len(z), 1)) # n_springs x 1
        Pu1 = (3 + (gamma * z / s_u) + (J * z / D)) * s_u * D # lateral bearing capacity factor (n_springs x 1)
        Pu2 = 9 * s_u * D
        Pu1[Pu1>Pu2] = Pu2
        p_ult = Pu1 # ultimate lateral soil resistance (n_springs x 1) for each depth
        
        p = np.zeros((len(y), len(p_ult)))
        for i in range(len(y)): # loop through displacements
            for j in range(len(p_ult)): # loop through p_ult for each spring at displacement i
                if y[i]/y_50 > 8: 
                    p[i,j] = p_ult[j] # rows displacements, columns lateral soil resistance/m for each depth
                else:
                    p[i,j] = 0.5 * p_ult[j] * ((y[i] / y_50) ** (1/3)) # lateral soil resistance (force/unit length) (displacements x n_springs)

        # now need to multiply p vector by length of element we want the spring to apply to since it is force/unit length
        # i.e. each p x gap between springs with half gap at base and surface of monopile
        p[:,1:-1] *= spring_gap # intermediate springs
        p[:,-1] *= spring_gap/2 # monopile surface spring
        p[:,0] *= spring_gap/2 # monopile base spring
        
        # for i in range(n_springs):
        #     plt.plot(y, p[:,i])
        # plt.show()

        # set up non-linear curve constants
        for i in range(n_springs):
            args = []
            args.append(i+1)
            for j in range(len(y)):
                args.append(y[j])
                args.append(p[j,i])
            #mapdl.r(i+1, y[0], p[0,i], y[1], p[1,i], y[2], p[2,i], y[3], p[3,i], y[4], p[4,i])
            #print(args)

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
    # # p-y API sand model
    # elif soil_model_switch == 3:
    #     mapdl.type(3)
    #     y = np.linspace(0,0.1,11) # lateral displacements
    #     D = mon_base_r2 * 2
    #     z = np.array(spring_depths)
    #     if scour_n_nodes > 0:
    #         scour_depth = scour_n_nodes * spring_gap # note number of nodes for monopile are important here to determine correct depth
    #         z -= scour_depth # update effective spring depths to account for scour
    #         z[z < 0] = 0
    #     #z = np.reshape(z, (len(z), 1)) # n_springs x 1
    #     Pus = ((soil_C1 * z) + (soil_C2 * D)) * soil_gamma * z # ultimate resistance at shallow depth
    #     Pud = soil_C3 * D * soil_gamma * z # ultimate resistance at the deep depth
    #     Pu = np.minimum(Pus, Pud)

    #     p = np.zeros((len(y), len(z)))
    #     for i in range(len(y)): # loop through displacements
    #         for j in range(len(z)): # loop through depth for each spring at displacement i
    #             if Pu[j] == 0:
    #                 p[i,j] = 0 # avoid divide by zero error (zero lateral resistance at surface for sand)
    #             else:   
    #                 p[i,j] = soil_A * Pu[j] * np.tanh(y[i] * soil_kpy * z[j] / (soil_A * Pu[j])) # lateral soil resistance (force/unit length) (displacements x n_springs)

    #     # now need to multiply p vector by length of element we want the spring to apply to since it is force/unit length
    #     # i.e. each p x gap between springs with half gap at base and surface of monopile
    #     p[:,1:-1] *= spring_gap # intermediate springs
    #     p[:,-1] *= spring_gap/2 # monopile surface spring (technically not stiff enough)
    #     p[:,0] *= spring_gap/2 # monopile base spring (technically too stiff)
        
    #     # # Plot spring curves
    #     # for i in range(n_springs):
    #     #     plt.plot(y, p[:,i])
    #     # plt.show()

    #     # set up non-linear curve constants
    #     for i in range(n_springs):
    #         args = []
    #         args.append(i+1)
    #         for j in range(len(y)):
    #             args.append(y[j])
    #             args.append(p[j,i])
    #         #mapdl.r(i+1, y[0], p[0,i], y[1], p[1,i], y[2], p[2,i], y[3], p[3,i], y[4], p[4,i])
    #         if len(args) <= 7:
    #             mapdl.r(*args)
    #             args = []
    #         else:
    #             mapdl.r(*args[0:7])
    #             args = args[7:]
    #         while len(args) >= 6:
    #             mapdl.rmore(*args[0:6])
    #             args = args[6:]
    #         if len(args) > 0:
    #             mapdl.rmore(*args)
    #             args = []
    #     # Rigid foundation flexible soil method - Not yet implemented - model used for 15MW NREL turbine        
    # elif soil_model_switch == 4:
    #     lat_mod = 1 + (0.55 * (2 - soil_poisson) * soil_embedded_depth / mon_base_r2)# lateral soil stiffness modifier
    #     soil_k_lat = 32*(1-soil_poisson)*soil_shear_mod*mon_base_r2*lat_mod/(7-(8*soil_poisson))# lateral soil stiffness
    #     vert_mod = 1 + (0.6 * (1 - soil_poisson) * soil_embedded_depth / mon_base_r2)
    #     soil_k_vert = 4*soil_shear_mod*mon_base_r2*vert_mod/(1-soil_poisson)
    #     rock_mod = 1 + (1.2 * (1 - soil_poisson) * soil_embedded_depth / mon_base_r2) + (0.2 * (2 - soil_poisson) * ((soil_embedded_depth / mon_base_r2)**3))
    #     soil_k_rock = 8*soil_shear_mod*(mon_base_r2 ** 3)*rock_mod/(3*(1-soil_poisson))
    #     soil_k_tors = 16*soil_shear_mod*(mon_base_r2 ** 3)/3    

    ## Create springs if needed
    # set up query
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
            mapdl.n("", i*spring_gap, spring_len, 0)
            # create spring
            mapdl.e(q.node(i*spring_gap, 0, 0), q.node(i*spring_gap, spring_len, 0))

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
            mapdl.n("", i*spring_gap, 0, spring_len)
            # create spring
            mapdl.e(q.node(i*spring_gap, 0, 0), q.node(i*spring_gap, 0, spring_len))

        # x-direction
        # vertical stiffness
        # create node
        mapdl.r(7, soil_k, soil_eta) # assumed linear spring for base
        mapdl.real(n_springs + 1) # change to correct constant set
        mapdl.type(2) # Change to COMBIN14
        #mapdl.dk(1,"UX",0)
        mapdl.n("", -spring_len, 0, 0)
        # create spring
        mapdl.e(q.node(0, 0, 0), q.node(-spring_len, 0, 0))

        # Add node constrains to free end of springs
        #mapdl.d(q.node(0, spring_len, 0),"ALL",0,"",q.node(spring_depths[0], 0, spring_len))
        mapdl.d(q.node(0, spring_len, 0),"ALL",0,"",q.node(-spring_len, 0, 0)) # first spring to last spring

    ## Apply force for wind loading
    mapdl.fk(6,"FY",-100)

    # # plotting
    # mapdl.eshape(1)
    # mapdl.eplot(show_edges=True, smooth_shading=True,
    #             show_node_numbering=True, plot_bc=True)
    # mapdl.nplot(nnum=True)
    # mapdl.finish()

    ### Analysis ###

    ## Nodes of interest
    # total number of nodes of interest
    node_count = mon_embedded_n_nodes + mon_submerged_n_nodes + mon_air_n_nodes + tower_n_nodes + nac_n_nodes 

    # get lengths of segments of monopile/tower/nacelle
    mon_embedded_node_len = soil_embedded_depth / mon_embedded_n_nodes
    mon_submerged_node_len = submerged_len / mon_submerged_n_nodes
    mon_air_node_len = (mon_len - soil_embedded_depth - submerged_len) / mon_air_n_nodes
    tower_node_len = tower_len / tower_n_nodes
    nac_node_len = nac_height / nac_n_nodes

    mon_n_nodes = mon_embedded_n_nodes + mon_submerged_n_nodes + mon_air_n_nodes
    # get x coordinates of nodes of interest
    nodal_x_coord = [0] * (node_count + 1)
    for i in range(node_count):
        if i < mon_embedded_n_nodes:
            nodal_x_coord[i+1] = nodal_x_coord[i] + mon_embedded_node_len
        elif i < mon_embedded_n_nodes + mon_submerged_n_nodes:
            nodal_x_coord[i+1] = nodal_x_coord[i] + mon_submerged_node_len
        elif i < mon_embedded_n_nodes + mon_submerged_n_nodes + mon_air_n_nodes:
            nodal_x_coord[i+1] = nodal_x_coord[i] + mon_air_node_len
        elif i < mon_n_nodes + tower_n_nodes:
            nodal_x_coord[i+1] = nodal_x_coord[i] + tower_node_len
        elif i < node_count:
            nodal_x_coord[i+1] = nodal_x_coord[i] + nac_node_len

    ## modal analysis
    mapdl.run("/SOL")
    mapdl.modal_analysis(nmode=20, freqb=0, freqe=10, mxpand=20)

    mapdl.finish()
    mapdl.post1()

    # natural frequencies
    nat_freqs = mapdl.post_processing.frequency_values

    # initialise matrices
    mode_shapes_x = np.zeros((10,len(nodal_x_coord)))
    mode_shapes_y = np.zeros((10,len(nodal_x_coord)))
    mode_shapes_z = np.zeros((10,len(nodal_x_coord)))

    # get nodal displacements
    for k in range(10):
        mapdl.set(sbstep=k+1) # for each solution substep
        for i in range(len(nodal_x_coord)):
            mode_shapes_x[k,i] = q.ux(q.node(nodal_x_coord[i], 0, 0))
            mode_shapes_y[k,i] = q.uy(q.node(nodal_x_coord[i], 0, 0))
            mode_shapes_z[k,i] = q.uz(q.node(nodal_x_coord[i], 0, 0))
    
    ## Harmonic Analysis
    freq_start = 0
    freq_end = 1
    freq_incs = 500 # set the number of increments between start and end

    mapdl.run("/SOLU")
    mapdl.antype(3)  # HARMONIC ANALYSIS
    mapdl.harfrq(freq_start, freq_end)  # Set beginning and ending frequency
    mapdl.nsubst(freq_incs)  # Set the number of frequency increments
    mapdl.solve()

    mapdl.finish()
    mapdl.post1()

    q = mapdl.queries

    # # All nodes
    # dy = np.zeros((freq_incs, len(nodal_x_coord)))

    # for i in range(freq_incs):
    #     mapdl.set(sbstep=i+1)
    #     for j in range(len(nodal_x_coord)):
    #         dy[i,j] = q.uy(q.node(nodal_x_coord[j], 0, 0))

    #End node only
    dy = np.zeros(freq_incs)
    for i in range(freq_incs):
        mapdl.set(sbstep=i+1)
        dy[i] = q.uy(q.node(nodal_x_coord[-1], 0, 0))


    ### playing with results file ###
    # mapdl.post26()
    # mapdl.file('file','rst')   
    # mapdl.numvar(200)  
    # mapdl.solu(191,"NCMIT")  
    # mapdl.store("MERGE") 
    # mapdl.plcplx(0)
    # mapdl.prcplx(1)
    # mapdl.filldata(191,"","","",1,1) 
    # mapdl.realvar(191,191) 
    # # !*  
    # mapdl.nsol(2,42,"U","Y", "UY_4")
    # mapdl.store("MERGE") 
    # mapdl.xvar(1)  
    # mapdl.plvar(2)
    # mapdl.prvar(2)
    # # )/GOP    ! Resume printing after UNDO process   


    print("analysis complete")
    mapdl.exit()
    
    return nat_freqs, nodal_x_coord, mode_shapes_x, mode_shapes_y, mode_shapes_z, dy


#### ------------- Sense checking ------------------
# Encastre
# nat_freqs, nodal_x_coord, mode_shapes_x, mode_shapes_y, mode_shapes_z, dy = wtg_foundation(
#                 steel_den=7850, steel_poisson=0.3, steel_YM=200e+9, steel_damping=0.02, # material parameters
#                    mon_len = 90, mon_base_r2 = 5, mon_base_r1 = 4, mon_top_r2 = 5, # monopile parameters
#                    mon_top_r1 = 4, mon_embedded_n_nodes=5, mon_submerged_n_nodes=3, mon_air_n_nodes = 2, 
#                    tower_len = 135, tower_top_r2 = 3, tower_top_r1 = 2, # tower parameters
#                    tower_base_r2 = 5, tower_base_r1 = 4, tower_n_nodes = 10,
#                    nac_mass = 789.168e+3, nac_height = 2, nac_n_nodes = 1, # nacelle parameters
#                    water_switch=0,
#                    fixtures_switch=0,
#                    soil_model_switch = 0,
#                    soil_embedded_depth = 45, soil_poisson = 0.4, soil_shear_mod = 50e+6, # soil parameters
#                    soil_k = 10000000, soil_eta = 0.01,
#                    s_u=50e+3, eps_50=0.008, gamma=8e+3, J=0.5, # p-y soft clay Matlock
#                    spring_len = 10, submerged_len = 20,
#                    scour_n_nodes = 0 # scour parameters
#                    )

# linear sprung
# nat_freqs, nodal_x_coord, mode_shapes_x, mode_shapes_y, mode_shapes_z, dy = wtg_foundation(
#                 steel_den=7850, steel_poisson=0.3, steel_YM=200e+9, steel_damping=0.02, # material parameters
#                    mon_len = 90, mon_base_r2 = 5, mon_base_r1 = 4, mon_top_r2 = 5, # monopile parameters
#                    mon_top_r1 = 4, mon_embedded_n_nodes=5, mon_submerged_n_nodes=3, mon_air_n_nodes = 2, 
#                    tower_len = 135, tower_top_r2 = 3, tower_top_r1 = 2, # tower parameters
#                    tower_base_r2 = 5, tower_base_r1 = 4, tower_n_nodes = 10,
#                    nac_mass = 789.168e+3, nac_height = 2, nac_n_nodes = 1, # nacelle parameters
#                    water_switch=0,
#                    fixtures_switch=0,
#                    soil_model_switch = 1,
#                    soil_embedded_depth = 45, soil_poisson = 0.4, soil_shear_mod = 50e+6, # soil parameters
#                    soil_k = 10000000, soil_eta = 0.01,
#                    s_u=50e+3, eps_50=0.008, gamma=8e+3, J=0.5, # p-y soft clay Matlock
#                    spring_len = 10, submerged_len = 20,
#                    scour_n_nodes = 0 # scour parameters
#                    )

# non-linear sprung
# nat_freqs, nodal_x_coord, mode_shapes_x, mode_shapes_y, mode_shapes_z, dy = wtg_foundation(
#                 steel_den=7850, steel_poisson=0.3, steel_YM=200e+9, steel_damping=0.02, # material parameters
#                    mon_len = 90, mon_base_r2 = 5, mon_base_r1 = 4, mon_top_r2 = 5, # monopile parameters
#                    mon_top_r1 = 4, mon_embedded_n_nodes=5, mon_submerged_n_nodes=3, mon_air_n_nodes = 2, 
#                    tower_len = 135, tower_top_r2 = 3, tower_top_r1 = 2, # tower parameters
#                    tower_base_r2 = 5, tower_base_r1 = 4, tower_n_nodes = 10,
#                    nac_mass = 789.168e+3, nac_height = 2, nac_n_nodes = 1, # nacelle parameters
#                    water_switch=0,
#                    fixtures_switch=0,
#                    soil_model_switch = 2,
#                    soil_embedded_depth = 45, soil_poisson = 0.4, soil_shear_mod = 50e+6, # soil parameters
#                    soil_k = 10000000, soil_eta = 0.01,
#                    s_u=50e+3, eps_50=0.008, gamma=8e+3, J=0.5, # p-y soft clay Matlock
#                    spring_len = 10, submerged_len = 20,
#                    scour_n_nodes = 0 # scour parameters
#                    )


### -------- Use cases -----------

# # non-linear sprung
nat_freqs, nodal_x_coord, mode_shapes_x, mode_shapes_y, mode_shapes_z, dy = wtg_foundation(
                steel_den=7850, steel_poisson=0.3, steel_YM=200e+9, steel_damping=0.02, # material parameters
                   mon_len = 90, mon_base_r2 = 5, mon_base_r1 = 4, mon_top_r2 = 5, # monopile parameters
                   mon_top_r1 = 4, mon_embedded_n_nodes=5, mon_submerged_n_nodes=3, mon_air_n_nodes = 2, 
                   tower_len = 135, tower_top_r2 = 3, tower_top_r1 = 2, # tower parameters
                   tower_base_r2 = 5, tower_base_r1 = 4, tower_n_nodes = 10,
                   nac_mass = 789.168e+3, nac_height = 2, nac_n_nodes = 1, # nacelle parameters
                   water_switch=1,
                   fixtures_switch=1,
                   soil_model_switch = 2,
                   soil_embedded_depth = 45, soil_poisson = 0.4, soil_shear_mod = 50e+6, # soil parameters
                   soil_k = 10000000, soil_eta = 0.01,
                   s_u=50e+3, eps_50=0.008, gamma=8e+3, J=0.5, # p-y soft clay Matlock
                   spring_len = 10, submerged_len = 20,
                   scour_n_nodes = 0 # scour parameters
                   )


# 15 MW NREL Ref turbine - soil stiffness to be matched - based on dense sand or gravel soils (35 - 200MPa for shear mod)
# nat_freqs, nodal_x_coord, mode_shapes_x, mode_shapes_y, mode_shapes_z = wtg_foundation(steel_den=7850, steel_poisson=0.3, steel_YM=200e+9, steel_damping=0.02, # material parameters
#                    mon_len = 90, mon_base_r2 = 5, mon_base_r1 = 4.945, mon_top_r2 = 5, # monopile parameters
#                    mon_top_r1 = 4.959, mon_n_nodes = 10,
#                    tower_len = 130, tower_top_r2 = 3.25, tower_top_r1 = 3.226, # tower parameters
#                    tower_base_r2 = 5, tower_base_r1 = 4.959, tower_n_nodes = 10,
#                    nac_mass = 1017e+3, nac_height = 2, nac_n_nodes = 1, # nacelle parameters
#                    soil_embedded_depth = 45, soil_poisson = 0.4, soil_shear_mod = 140e+6, # soil parameters
#                    soil_model_switch = 3, soil_kpy=200e+3, # p-y curve soil model
#                    spring_len = 10, submerged_len = 30,
#                    scour_n_nodes=0)

# # 5 MW NREL ref turbine
# nat_freqs, nodal_x_coord, mode_shapes_x, mode_shapes_y, mode_shapes_z, dy = wtg_foundation(steel_den=7850, steel_poisson=0.3, steel_YM=210e+9, steel_damping=0.01, # material parameters
#                    mon_len = 75, mon_base_r2 = 3, mon_base_r1 = 2.973, mon_top_r2 = 3, # monopile parameters
#                    mon_top_r1 = 2.973, mon_embedded_n_nodes=5, mon_submerged_n_nodes=3, mon_air_n_nodes=2,
#                    tower_len = 90, tower_top_r2 = 1.935, tower_top_r1 = 1.916, # tower parameters
#                    tower_base_r2 = 3, tower_base_r1 = 2.973, tower_n_nodes = 10,
#                    nac_mass = 350e+3, nac_height = 2, nac_n_nodes = 1, # nacelle parameters
#                    soil_embedded_depth = 45, soil_model_switch = 1,
#                    s_u=50e+3, eps_50=0.008, gamma=8e+3, J=0.5, # p-y soft clay Matlock
#                    spring_len = 10, submerged_len = 20,
#                    scour_n_nodes = 0 # scour parameters
#                    )


# 5 MW NREL ref turbine - Rayleigh Ritz Validation
# nat_freqs, nodal_x_coord, mode_shapes_x, mode_shapes_y, mode_shapes_z, dy = wtg_foundation(steel_den=7850, steel_poisson=0.3, steel_YM=210e+9, steel_damping=0.01, # material parameters
#                    mon_len = 75, mon_base_r2 = 3, mon_base_r1 = 2.973, mon_top_r2 = 3, # monopile parameters
#                    mon_top_r1 = 2.973, mon_embedded_n_nodes=5, mon_submerged_n_nodes=3, mon_air_n_nodes=2,
#                    tower_len = 90, tower_top_r2 = 1.935, tower_top_r1 = 1.916, # tower parameters
#                    tower_base_r2 = 3, tower_base_r1 = 2.973, tower_n_nodes = 10,
#                    nac_mass = 0.001, nac_height = 0.001, nac_n_nodes = 1, # nacelle parameters
#                    soil_model_switch = 0, soil_embedded_depth = 45, submerged_len=20
#                    )

# # population testing - small: YM, damping, density, big: soil properties
# for i in range(5):
#     nat_freqs, nodal_x_coord, mode_shapes_x, mode_shapes_y, mode_shapes_z, dy = wtg_foundation(steel_den=7850, steel_poisson=0.3, steel_YM=np.random.normal(210e+9,3e+9), steel_damping=0.01, # material parameters
#                    mon_len = 75, mon_base_r2 = 3, mon_base_r1 = 2.973, mon_top_r2 = 3, # monopile parameters
#                    mon_top_r1 = 2.973, mon_embedded_n_nodes=5, mon_submerged_n_nodes=3, mon_air_n_nodes=2,
#                    tower_len = 90, tower_top_r2 = 1.935, tower_top_r1 = 1.916, # tower parameters
#                    tower_base_r2 = 3, tower_base_r1 = 2.973, tower_n_nodes = 10,
#                    nac_mass = 350e+3, nac_height = 2, nac_n_nodes = 1, # nacelle parameters
#                    soil_embedded_depth = 45, soil_model_switch = 1,
#                    s_u=np.random.normal(50e+3,10e+3), eps_50=0.008, gamma=8e+3, J=0.5, # p-y soft clay Matlock
#                    spring_len = 10, submerged_len=20,
#                    scour_n_nodes = 0 # scour parameters
#                    )
#     freq_start = 0
#     freq_end = 0.5
#     freq_incs = 500
#     print("natural frequencies")
#     print(nat_freqs)
#     plt.plot(np.linspace(freq_start, freq_end, freq_incs),abs(dy[:]))
# plt.show()


## Mode shape plotting
# # Initialise the subplot function using number of rows and columns
# figure, axis = plt.subplots(3, sharex=True, sharey=True) # share x and y axes
# for k in range(10):
#     plt.subplot(311)
#     plt.plot(nodal_x_coord, mode_shapes_x[k,:], label="mode " + str(k+1))
#     plt.subplot(312)
#     plt.plot(nodal_x_coord, mode_shapes_y[k,:], label="mode " + str(k+1))
#     plt.subplot(313)
#     plt.plot(nodal_x_coord, mode_shapes_z[k,:], label="mode " + str(k+1))

# plt.subplot(311)
# plt.title("mode shape UX")
# plt.legend()
# plt.subplot(312)
# plt.title("mode shape UY")
# plt.legend()
# plt.subplot(313)
# plt.title("mode shape UZ")
# plt.legend()
# plt.show()

# # FRF plotting
# freq_start = 0
# freq_end = 3
# freq_incs = 300
# # Plot all nodes
# # for k in range(len(nodal_x_coord)):
# #     plt.plot(np.linspace(freq_start, freq_end, freq_incs),abs(dy[:,k]))

# # Plot end node
# plt.plot(np.linspace(freq_start, freq_end, freq_incs),abs(dy[:]))
# plt.show()

freq_start = 0
freq_end = 1
freq_incs = 500
print("natural frequencies")
print(nat_freqs)
plt.plot(np.linspace(freq_start, freq_end, freq_incs),dy[:])
plt.show()