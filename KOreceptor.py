# receptor KO / KD

import smoldyn
import random
import math
from pathlib import Path
from datetime import datetime
from multiprocessing import Pool

scriptdir = Path(__file__).parent
MAX_CPUS = 56


def build_model_smoldyn(exo_rate=0.016, R_exo_rate=0.016, endo_rate=0.06, time_step=0.01, run_time=1.0):
    start_time = datetime.now()
    print(f"smoldyn: {smoldyn.version()}, building simulation: {start_time.strftime('%H:%M:%S')}")
    s = smoldyn.Simulation(low=[0, 0], high=[100, 20])

    # SPECIES

    ########## CREATING VESICLES
    vesicle_F = s.addSpecies("vesicle_F", difc=dict(all=50), color="yellow", display_size=1)
    vesicle_F.addToSolution(1, lowpos=(0, 0), highpos=(100, 5))
    s.setSpeciesMobility("vesicle_F", smoldyn._smoldyn.MolecState.soln, diffConst=50, drift=[0, 50], difmatrix=[])

    vesicle_F_top = s.addSpecies("vesicle_F_top", difc=dict(all=50), color="yellow", display_size=1)
    vesicle_F_top.addToSolution(1, lowpos=(0, 15), highpos=(100, 20))
    s.setSpeciesMobility("vesicle_F_top", smoldyn._smoldyn.MolecState.soln, diffConst=50, drift=[0, -50], difmatrix=[])


    #


    vesicle_R = s.addSpecies("vesicle_R", difc=dict(all=50), color="orange", display_size=1)
    vesicle_R_top = s.addSpecies("vesicle_R_top", difc=dict(all=50), color="orange", display_size=1)
    vesicle_R_top.addToSolution(1, lowpos=(0, 15), highpos=(100, 20))
    vesicle_R.addToSolution(1, lowpos=(0, 0), highpos=(100, 5))
    s.setSpeciesMobility("vesicle_R", smoldyn._smoldyn.MolecState.soln, diffConst=50, drift=[0, 50], difmatrix=[])
    s.setSpeciesMobility("vesicle_R_top", smoldyn._smoldyn.MolecState.soln, diffConst=50, drift=[0, -50], difmatrix=[])

    vesicle_H = s.addSpecies("vesicle_H", difc=dict(all=50), color="coral", display_size=1)
    vesicle_H_top= s.addSpecies("vesicle_H_top", difc=dict(all=50), color="coral", display_size=1)
    vesicle_H_top.addToSolution(1, lowpos=(0, 15), highpos=(100, 20))
    vesicle_H.addToSolution(1, lowpos=(0, 0), highpos=(100, 5))
    s.setSpeciesMobility("vesicle_H", smoldyn._smoldyn.MolecState.soln, diffConst=50, drift=[0, 50], difmatrix=[])
    s.setSpeciesMobility("vesicle_H_top", smoldyn._smoldyn.MolecState.soln, diffConst=50, drift=[0, -50], difmatrix=[])
    ##########

    ########## VESICLE FUSION
    fused_vesicle_F = s.addSpecies("fused_vesicle_F", difc=dict(all=10), color="purple", display_size=1)
    fused_vesicle_R = s.addSpecies("fused_vesicle_R", difc=dict(all=10), color="purple", display_size=1)
    fused_vesicle_H = s.addSpecies("fused_vesicle_H", difc=dict(all=10), color="purple", display_size=1)
    ##########

    ########## INITIALISING MOLECULES
    F = s.addSpecies("F", difc=dict(all=200000), color="red", display_size=1)
    H = s.addSpecies("H", difc=dict(all=316), color=dict(all="green"), display_size=1)
    R = s.addSpecies("R", difc=dict(all=685), color=dict(all="blue"), display_size=1)
    RR = s.addSpecies("RR", difc=dict(all=316), color=dict(all="blue"), display_size=1)
    FH = s.addSpecies("FH", difc=dict(all=316), color=dict(all="pink"), display_size=1)
    FHRR = s.addSpecies("FHRR", difc=dict(all=316), color=dict(all="purple"), display_size=1)
    ##########

    
    ########## BOUNDARIES
    r1 = smoldyn.Rectangle(corner=[0, 0], dimensions=[20], axis="+0", name="r1")
    r2 = smoldyn.Rectangle(corner=[100, 0], dimensions=[20], axis="-0", name="r2")
    r3 = smoldyn.Rectangle(corner=[0, 0], dimensions=[100], axis="+1", name="r3")
    r4 = smoldyn.Rectangle(corner=[0, 20], dimensions=[100], axis="-1", name="r4")
    sides = s.addSurface("sides", panels=[r1, r2])
    top_and_bottom = s.addSurface("top_and_bottom", panels=[r3, r4])
    top_and_bottom.setStyle('both', "black")
    top_and_bottom.setAction('both', [F, H, R, RR, FH, FHRR, vesicle_F, vesicle_F_top, vesicle_R, vesicle_R_top, vesicle_H, vesicle_H_top, fused_vesicle_H, fused_vesicle_R, fused_vesicle_F], 'reflect')

    sides.setAction('front', [F, H, R, RR, FH, FHRR, vesicle_F, vesicle_F_top, vesicle_R,vesicle_R_top, vesicle_H, vesicle_H_top, fused_vesicle_H, fused_vesicle_R, fused_vesicle_F], 'jump')
    sides.setAction('back', [F, H, R, RR, FH, FHRR, vesicle_F, vesicle_F_top, vesicle_R, vesicle_R_top, vesicle_H, vesicle_H_top,  fused_vesicle_H, fused_vesicle_R, fused_vesicle_F], 'reflect')

    sides.setStyle('front', color="red")
    sides.setStyle('back', color="blue")
    sides.setStyle('both', thickness=1)

    r1.jumpTo('front', r2, 'front', True)
    ##########

    ########## bottom cell PM
    rect1 = smoldyn.Rectangle(corner=(0, 5), dimensions=[100], axis="+y")
    bottom_PM = s.addSurface("bottom_PM", panels=[rect1])
    bottom_PM.setStyle('both', color="black", thickness=1)
    bottom_PM.setAction('both', [F, H, R, RR, FH, FHRR, vesicle_F, vesicle_F_top, vesicle_R, vesicle_R_top, vesicle_H, vesicle_H_top, fused_vesicle_H, fused_vesicle_R, fused_vesicle_F], "reflect")
    bottom_cell = s.addCompartment(name="bottom_cell", surface="bottom_PM", point=[1, 1])
    ##########

    ########## top cell PM
    rect2 = smoldyn.Rectangle(corner=(0, 15), dimensions=[100], axis="+y")
    top_PM = s.addSurface("top_PM", panels=[rect2])
    top_PM.setStyle('both', color="black", thickness=1)
    top_PM.setAction('both', [F, H, R, RR, FH, FHRR, vesicle_F, vesicle_F_top, vesicle_R, vesicle_R_top, vesicle_H, vesicle_H_top, fused_vesicle_H, fused_vesicle_R, fused_vesicle_F], "reflect")
    top_cell = s.addCompartment(name="top_cell", surface="top_PM", point=[99, 19])
    ##########

    ######### REACTIONS ##############################################

    ######### setting up vesicle production and release, and species decay 
    F_vesicle_production = s.addReaction("F_vesicle_production", [], [vesicle_F], rate=exo_rate)
    F_vesicle_production_top = s.addReaction("F_vesicle_production_top", [], [vesicle_F_top], rate=exo_rate)
    bottom_PM.setRate(vesicle_F, 'bsoln', 'back', rate=1000, revrate=0)
    bottom_PM.setRate(vesicle_F, 'back', 'fsoln', rate=1000, new_species=fused_vesicle_F)
    top_PM.setRate(vesicle_F_top, 'fsoln', 'front', rate=1000, revrate=0)
    top_PM.setRate(vesicle_F_top, 'front', 'bsoln', rate=1000, new_species=fused_vesicle_F)
    release_F = s.addReaction("release_F", subs=[fused_vesicle_F], prds=[F] * random.randint(10, 20), rate=1000)
    F_decay = s.addReaction("F_decay", subs=[F], prds=[], rate=endo_rate)

    R_vesicle_production = s.addReaction("R_vesicle_production", [], [vesicle_R], rate=R_exo_rate)
    s.addCommand("killmol vesicle_R", cmd_type='i', on=30, off=run_time, step=0.01)
    s.addCommand("killmol vesicle_R_top", cmd_type='i', on=30, off=run_time, step=0.01)
    R_vesicle_production_top = s.addReaction("R_vesicle_production_top", [], [vesicle_R_top], rate=R_exo_rate)
    top_PM.setRate(vesicle_R_top, 'fsoln', 'front', rate=1000, revrate=0)
    top_PM.setRate(vesicle_R_top, 'front', 'down', rate=1000, new_species=fused_vesicle_R)
    bottom_PM.setRate(vesicle_R, 'bsoln', 'back', rate=1000, revrate=0)
    bottom_PM.setRate(vesicle_R, 'back', 'up', rate=1000, new_species=fused_vesicle_R)
    release_R_top = s.addReaction("release_R_top", subs=[(fused_vesicle_R, 'down')], prds=[(R, 'down')] * random.randint(10, 20), rate=1000)
    release_R_bottom = s.addReaction("release_R_bottom", subs=[(fused_vesicle_R, 'up')], prds=[(R, 'up')] * random.randint(10, 20), rate=1000)
    #R_decay = s.addReaction("R_decay", subs=[R], prds=[], rate=0.8)

    H_vesicle_production = s.addReaction("H_vesicle_production", [], [vesicle_H], rate=exo_rate)
    H_vesicle_production_top = s.addReaction("H_vesicle_production_top", [], [vesicle_H_top], rate=exo_rate)
    top_PM.setRate(vesicle_H_top, 'fsoln', 'front', rate=1000, revrate=0)
    top_PM.setRate(vesicle_H_top, 'front', 'down', rate=1000, new_species=fused_vesicle_H)
    bottom_PM.setRate(vesicle_H, 'bsoln', 'back', rate=1000, revrate=0)
    bottom_PM.setRate(vesicle_H, 'back', 'up', rate=1000, new_species=fused_vesicle_H)
    release_H_top = s.addReaction("release_H_top", subs=[(fused_vesicle_H, 'down')], prds=[(H, 'down')] * random.randint(10, 20), rate=1000)
    release_H_bottom = s.addReaction("release_H_bottom", subs=[(fused_vesicle_H, 'up')], prds=[(H, 'up')] * random.randint(10, 20), rate=1000)
    #H_decay = s.addReaction("H_decay", subs=[H], prds=[], rate=0.8)
    ##########

    ########## F + H <-> FH 
    F_to_H_top = s.addBidirectionalReaction("F_to_H_top", subs=[(F, 'bsoln'), (H, 'down')], prds=[(FH, 'down')], kf=2, kb=0, surface=top_PM)
    F_to_H_bottom = s.addBidirectionalReaction("F_to_H_bottom", subs=[(F, 'fsoln'), (H, 'up')], prds=[(FH, 'up')], kf=2, kb=0, surface=bottom_PM)
    ##########

    ########## R + R <-> RR
    R_to_R_top = s.addBidirectionalReaction("R_to_R_top", subs=[(R, 'down'), (R, 'down')], prds=[(RR, 'down')], kf=10, kb=0, surface=top_PM)
    R_to_R_bottom = s.addBidirectionalReaction("R_to_R_bottom", subs=[(R, 'up'), (R, 'up')], prds=[(RR, 'up')], kf=10, kb=0, surface=bottom_PM)
    ##########

    ########## FH + RR -> FHRR
    FH_to_RR_top = s.addReaction("FH_to_RR_top", subs=[(FH, 'down'), (RR, 'down')], prds=[(FHRR, 'down')], rate = 190000, surface=top_PM)
    FH_to_RR_bottom = s.addReaction("FH_to_RR_bottom", subs=[(FH, 'up'), (RR, 'up')], prds=[(FHRR, 'up')], rate = 190000, surface=bottom_PM)
    ##########

    ########## endocytosis of all molecules
    H_endo_top = s.addReaction("H_endo_top", subs=[(H, 'down')], prds=[(H, 'fsoln')], rate=endo_rate)
    H_decay_top = s.addReaction("H_decay_top", subs=[(H, 'fsoln')], prds=[], rate=100) 
    H_endo_bottom = s.addReaction("H_endo_bottom", subs=[(H, 'up')], prds=[(H, 'bsoln')], rate=endo_rate)
    H_decay_bottom = s.addReaction("H_decay_bottom", subs=[(H, 'bsoln')], prds=[], rate=100) 

    R_endo_top = s.addReaction("R_endo_top", subs=[(R, 'down')], prds=[(R, 'fsoln')], rate=endo_rate)
    R_decay_top = s.addReaction("R_decay_top", subs=[(R, 'fsoln')], prds=[], rate=100) 
    R_endo_bottom = s.addReaction("R_endo_bottom", subs=[(R, 'up')], prds=[(R, 'bsoln')], rate=endo_rate)
    R_decay_bottom = s.addReaction("R_decay_bottom", subs=[(R, 'bsoln')], prds=[], rate=100) 

    RR_endo_top = s.addReaction("RR_endo_top", subs=[(RR, 'down')], prds=[(RR, 'fsoln')], rate=endo_rate)
    RR_decay_top = s.addReaction("RR_decay_top", subs=[(RR, 'fsoln')], prds=[], rate=100) 
    RR_endo_bottom = s.addReaction("RR_endo_bottom", subs=[(RR, 'up')], prds=[(RR, 'bsoln')], rate=endo_rate)
    RR_decay_bottom = s.addReaction("RR_decay_bottom", subs=[(RR, 'bsoln')], prds=[], rate=100) 

    FH_endo_top = s.addReaction("FH_endo_top", subs=[(FH, 'down')], prds=[(FH, 'fsoln')], rate=endo_rate)
    FH_decay_top = s.addReaction("FH_decay_top", subs=[(FH, 'fsoln')], prds=[], rate=100) 
    FH_endo_bottom = s.addReaction("FH_endo_bottom", subs=[(FH, 'up')], prds=[(FH, 'bsoln')], rate=endo_rate)
    FH_decay_bottom = s.addReaction("FH_decay_bottom", subs=[(FH, 'bsoln')], prds=[], rate=100) 
    

    ########## FHRR endocytosis + decay 
    FHRR_endo_top = s.addReaction("FHRR_endo_top", subs=[(FHRR, 'down')], prds=[(FHRR, 'fsoln')], rate=endo_rate)
    FHRR_decay_top = s.addReaction("FHRR_decay_top", subs=[(FHRR, 'fsoln')], prds=[], rate=100) 
    FHRR_endo_bottom = s.addReaction("FHRR_endo_bottom", subs=[(FHRR, 'up')], prds=[(FHRR, 'bsoln')], rate=endo_rate)
    FHRR_decay_bottom = s.addReaction("FHRR_decay_bottom", subs=[(FHRR, 'bsoln')], prds=[], rate=100) 


    datafile1 = f"20240403_receptorKO.csv"
    # datafile2 = f"F_coordinates-endo_rate={endo_rate}.csv"
    # datafile3 = f"H_coordinates-endo_rate={endo_rate}.csv"
    # datafile4 = f"R_coordinates-endo_rate={endo_rate}.csv"
    # datafile5 = f"FH_coordinates-endo_rate={endo_rate}.csv"
    # datafile6 = f"RR_coordinates-endo_rate={endo_rate}.csv"
    # datafile7 = f"FHRR_coordinates-endo_rate={endo_rate}.csv"
    s.addOutputFile(datafile1, 0, 0)
    # s.addOutputFile(datafile2, 0, 0)
    # s.addOutputFile(datafile3, 0, 0)
    # s.addOutputFile(datafile4, 0, 0)
    # s.addOutputFile(datafile5, 0, 0)
    # s.addOutputFile(datafile6, 0, 0)
    # s.addOutputFile(datafile7, 0, 0)
    s.addCommand(f"listmols2 {datafile1}", "i", on=0, off=run_time, step=0.1)
    # s.addCommand(f"listmols3 F {datafile2}", "i", on=0, off=run_time, step=0.1)
    # s.addCommand(f"listmols3 H {datafile3}", "i", on=0, off=run_time, step=0.1)
    # s.addCommand(f"listmols3 R {datafile4}", "i", on=0, off=run_time, step=0.1)
    # s.addCommand(f"listmols3 FH {datafile5}", "i", on=0, off=run_time, step=0.1)
    # s.addCommand(f"listmols3 RR {datafile6}", "i", on=0, off=run_time, step=0.1)
    # s.addCommand(f"listmols3 FHRR {datafile7}", "i", on=0, off=run_time, step=0.1)



    #s.addGraphics("opengl_good", iter=10, text_display="time")	# Graphics need to be turned off
    print('[INFO] Starting simulation ...')
    s.run(stop=run_time, dt=time_step, overwrite=True)  
    end_time = datetime.now()
    delta = end_time - start_time
    print(f"Done, ran for {round(delta.seconds / 60, 2)} minutes ({endo_rate}, {exo_rate}, {time_step}, {run_time})")
    print(smoldyn.getError())


def main(args):
    (endo_rate, time_step, run_time) = args
    build_model_smoldyn(endo_rate=endo_rate, time_step=time_step, run_time=run_time)
 


if __name__ == "__main__":
    args = []
    for endo_rate in [0.067]:
       for (time_step, run_time) in [(0.01, 180.0)]:
           args.append((endo_rate, time_step, run_time))
           
    with Pool(min(len(args), MAX_CPUS)) as pool:
        pool.map(main, args)



