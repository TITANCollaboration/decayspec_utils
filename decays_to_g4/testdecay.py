import time
import sys
import numpy as np

from pyne import nucname
from pyne.material import Material
from pyne.transmute.chainsolve import Transmuter
from pprint import pprint
from pyne import data

tm = Transmuter(tol=1e-10, log=sys.stdout)
#inp = Material({'In129': 1.0}, mass=1.0)
#decay_chain = tm.transmute(inp, t=60, tol=1e-8)

trap_contents = ['I129', 'Cs129', 'Sb129', 'Sb129m', 'Sn129', 'In129', 'In129m']
for my_isotope in trap_contents:
    print("Isotope:", my_isotope)
    inp = Material({my_isotope: 1.0}, mass=1.0)
    decay_chain = tm.transmute(inp, t=60, tol=1e-8)
    #pprint(decay_chain)
    for decay_product in decay_chain:
        print(nucname.name(decay_product), " - ", decay_chain[decay_product])
    print("===================")
#        print(decay_chain[decay_product] * (1e8 / 2))

#print(data.gamma_energy(511290000))
#print("====================")
#print(data.gamma_parent(511290000))
#print(data.gamma_photon_intensity(511290000))
