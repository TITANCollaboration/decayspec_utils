#!/usr/bin/env python
# *************************************************************************************
#  Written by : Jon Ringuette
#  Purpose : Read in root file from the likes of geant4 and smear the energy out
#            so it has a more gaussian distribution and is more realistic
# *************************************************************************************
import argparse
import configparser
from array import array
import uproot
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt


def process_root_file(root_file_name):
    root_file = uproot.open(root_file_name)
    initdata = root_file["EVENT_NTUPLE"]["pulse_height"].array()
    newdata = gaussian_filter1d(initdata, sigma=.3)
    print(initdata)
    print(newdata)

    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    axs[0].hist(initdata, bins=200)
    axs[1].hist(newdata, bins=200)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description='Geant4 Macro Scheduler')

    parser.add_argument('--root_file', dest='root_file_name', required=True,
                        help="Root file to read in")
    parser.set_defaults(configFile="g4macsched.cfg")

    args, unknown = parser.parse_known_args()
    # processConfigFile(args.configFile)
    process_root_file(args.root_file_name)


if __name__ == "__main__":
    main()
