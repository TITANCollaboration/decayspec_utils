import pandas as pd
import math
from scipy.integrate import quad
from pyne import nucname
from pyne.material import Material
from pyne.transmute.chainsolve import Transmuter
from pprint import pprint
from pyne import data
import sys
from random import randrange

# Unprocessed data pasted from nudat2 gamma list to pandas csv :
# cat gammas_te129_unprocessed.csv  | awk '{print $1 ", " $3}' > gammas_te129.csv


class gammas:
    def __init__(self, macro_filename, time_in_trap_per_cycle, trap_cycles):
        self.macro_filename = macro_filename
        self.trap_cycles = trap_cycles
        self.time_in_trap_per_cycle = time_in_trap_per_cycle

        self.write_mode = 'w'
        return

    def nuclear_activity(self, half_life, particle_count):
        # Integrate over the time in the trap to get decays per trapping time
        # This functionality has been moved over to PyNE but keeping it here
        # for sanity checking
        lamb = math.log(2) / half_life
        my_activity = lambda time: (lamb * particle_count * math.exp(-lamb*time))
        decay_events = quad(my_activity, 0, self.time_in_trap_per_cycle)
        return decay_events[0]

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
        self.start_header = header + "_" + str(randrange(0,999999)) # "sb129m1_bg"

        df = pd.read_csv(gamma_file, sep=',')
        for i, gamma_info in df.iterrows():
            self.write_geant4_macro(gamma_info[0], gamma_info[1])
        return

    def write_gammas(self, decay_num, total_ions, decay_name, decay_gamma_file):
        decays = decay_num * total_ions
        print("Decays From:", decay_name, ":", "{:.2e}".format(decays))
        self.process_gamma_file(decay_gamma_file, decays, decay_name)

    def parent_decay_to_gammas(self, total_ions, decay_parent, decay_to_gammas):
        tm = Transmuter(tol=1e-10, log=sys.stdout)
        inp = Material({decay_parent: 1.0}, mass=1.0)
        decay_chain = tm.transmute(inp, t=self.time_in_trap_per_cycle, tol=1e-5)

        for decay_product in decay_chain:
            decay_name = nucname.name(decay_product)
            if decay_name in decay_to_gammas:
                self.write_gammas(decay_chain[decay_product], total_ions, decay_to_gammas[decay_name], "gammas_" + decay_to_gammas[decay_name] + ".csv")
        return


def sb129m1(mygamma, ion_count):
    decay_to_gammas = {'Sb129': "sb129m1_IT", 'Te129': "sb129m1_beta", 'I129': "te129"}
    mygamma.parent_decay_to_gammas(ion_count, 'Sb129m', decay_to_gammas)
    return


def sb129(mygamma, ion_count):
    decay_to_gammas = {'Te129': "sb129"}
    mygamma.parent_decay_to_gammas(ion_count, 'Sb129', decay_to_gammas)
    return


def cs129(mygamma, ion_count):
    decay_to_gammas = {'Xe129': "cs129"}
    mygamma.parent_decay_to_gammas(ion_count, 'Cs129', decay_to_gammas)
    return


def in129(mygamma, ion_count):
    decay_to_gammas = {'Sn129': "in129"}
    mygamma.parent_decay_to_gammas(ion_count, 'in129', decay_to_gammas)


def from_RFQ():
    mygamma = gammas("sb129_from_rfq_60s_60cycles_g4_background.mac",
                     time_in_trap_per_cycle=3600,
                     trap_cycles=1)
    sb129m1_rate = 1e7
    trap_contents = {cs129: sb129m1_rate * (10/4), sb129: sb129m1_rate / 4, sb129m1: sb129m1_rate}
    #trap_contents = {sb129m1: 2e5}
    for my_isotope in trap_contents:
        my_isotope(mygamma, trap_contents[my_isotope])  # This is using the dict key as a function call
    return


#mygamma = gammas("sb129_from_mrtof_2e5_1hr_g4_background.mac",
#                 time_in_trap_per_cycle=60,
#                 trap_cycles=1)
#sb129m1(mygamma, 2e5)
#sb129(mygamma, 2e5)


from_RFQ()
#print(mygamma.nuclear_activity(17.7*60, 1e7))
