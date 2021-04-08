import pandas as pd
import math
from scipy.integrate import quad
from pyne import nucname
from pyne.material import Material
from pyne.transmute.chainsolve import Transmuter
from pprint import pprint
from pyne import data
import sys
# Unprocessed data pasted from nudat2 gamma list to pandas csv :
# cat gammas_te129_unprocessed.csv  | awk '{print $1 ", " $3}' > gammas_te129.csv

class gammas:
    def __init__(self, macro_filename):
        self.macro_filename = macro_filename
        self.trap_cycles = 0
        self.write_mode = 'w'
        return

#    def nuclear_activity(self, time_in_trap, half_life, particle_count):
        # Integrate over the time in the trap to get decays per trapping time
#        lamb = math.log(2) / half_life
#        my_activity = lambda time: (lamb * particle_count * math.exp(-lamb*time))
#        decay_events = quad(my_activity, 0, time_in_trap)
#        return decay_events[0]

    def write_geant4_macro(self, gamma_energy, probability_of_decay):
        with open(self.macro_filename, self.write_mode) as geantMacroFile:
            num_decays = round((probability_of_decay / 100) * (self.initial_decays))
            num_decays = num_decays * self.trap_cycles
            if num_decays > 0:
                geantMacroFile.write("# START: %s_%s\n" % (self.start_header, gamma_energy))
                geantMacroFile.write("/pga/selectGunAction 1\n")
                geantMacroFile.write("/gps/particle gamma\n")
                geantMacroFile.write("/gps/position 0 0 0\n")
                geantMacroFile.write("/gps/energy %s keV\n" % (gamma_energy))
                geantMacroFile.write("/gps/ang/type iso\n")
                geantMacroFile.write("/histo/filename %s_%s\n" % (self.start_header, gamma_energy))
                geantMacroFile.write("/run/beamOn %s\n" % (num_decays))
                geantMacroFile.write("# END\n")
        self.write_mode = 'a'

    def process_gamma_file(self, gamma_file, initial_decays, header):
        self.initial_decays = initial_decays
        self.start_header = header #"sb129m1_bg"

        df = pd.read_csv(gamma_file, sep=',')
        for i, gamma_info in df.iterrows():
            self.write_geant4_macro(gamma_info[0], gamma_info[1])
        return

    def write_gammas(self, decay_num, total_ions, trap_cycles, decay_name, decay_gamma_file):
        self.trap_cycles = trap_cycles
        decays = decay_num * total_ions
        print("To", decay_name, ":", decays)
        self.process_gamma_file(decay_gamma_file, decays, decay_name)


def parent_decay_to_gammas(mygamma, total_ions, time_in_trap_per_cycle, trap_cycles, decay_parent, decay_to_gammas):
    tm = Transmuter(tol=1e-10, log=sys.stdout)
    inp = Material({decay_parent: 1.0}, mass=1.0)
    decay_chain = tm.transmute(inp, t=time_in_trap_per_cycle, tol=1e-5)

    for decay_product in decay_chain:
        decay_name = nucname.name(decay_product)
        if decay_name in decay_to_gammas:
            mygamma.write_gammas(decay_chain[decay_product], total_ions, trap_cycles, decay_to_gammas[decay_name], "gammas_" + decay_to_gammas[decay_name] + ".csv")
    return


def sb129m1(mygamma, total_ions, time_in_trap_per_cycle, trap_cycles):
    decay_to_gammas = {'Sb129': "sb129m1_IT", 'Te129': "sb129m1_beta", 'I129': "te129"}
    parent_decay_to_gammas(mygamma, total_ions, time_in_trap_per_cycle, trap_cycles, 'Sb129m', decay_to_gammas)
    return

def sb129(mygamma, total_ions, time_in_trap_per_cycle, trap_cycles):
    decay_to_gammas = {'Te129': "sb129"}
    parent_decay_to_gammas(mygamma, total_ions, time_in_trap_per_cycle, trap_cycles, 'Sb129', decay_to_gammas)
    return

def cs129(mygamma, total_ions, time_in_trap_per_cycle, trap_cycles):
    decay_to_gammas = {'Xe129': "cs129"}
    parent_decay_to_gammas(mygamma, total_ions, time_in_trap_per_cycle, trap_cycles, 'Cs129', decay_to_gammas)
    return

def in129(mygamma, total_ions, time_in_trap_per_cycle, trap_cycles):
    decay_to_gammas = {'Sn129': "in129"}
    parent_decay_to_gammas(mygamma, total_ions, time_in_trap_per_cycle, trap_cycles, 'in129', decay_to_gammas)

def mrtof_observed_beam(total_ions, time_in_trap_per_cycle, trap_cycles):
    # trap_contents = ['I129', 'Cs129', 'Sb129', 'Sb129m', 'Sn129', 'In129', 'In129m']
    trap_contents = ['Cs129']#, 'Sb129', 'Sb129m']

    mygamma = gammas("sb129_from_rfq_60s_g4_background.mac")
    tm = Transmuter(tol=1e-10, log=sys.stdout)
    for my_isotope in trap_contents:
        print("Isotope:", my_isotope)
        inp = Material({my_isotope: 1.0}, mass=1.0)
        decay_chain = tm.transmute(inp, t=60, tol=1e-8)
        #pprint(decay_chain)
        for decay_product in decay_chain:
            print(nucname.name(decay_product), " - ", decay_chain[decay_product])

    return

mygamma = gammas("sb129_from_mrtof_2e5_1hr_g4_background.mac")
#sb129m1(mygamma, 2e5, 60, 60)
sb129(mygamma, 2e5, 60, 1)


#mrtof_observed_beam(2e8, 60, 1)
