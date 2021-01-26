#!/usr/bin/env python
# *************************************************************************************
#  Written by : Jon Ringuette
#  Purpose : Read in root file from the likes of geant4 and smear the energy out
#            so it has a more gaussian distribution and is more realistic
# *************************************************************************************
import argparse
import uproot
import numpy as np
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
import pandas as pd


def process_root_file(root_file_name, output_file):
    root_file = uproot.open(root_file_name)
    initdata = root_file["EVENT_NTUPLE"]["pulse_height"].array().astype(int)
    newdata = gaussian_filter1d(initdata, sigma=.3)
    print(initdata)
    print(newdata)
    if output_file is not None:
        mydict = {'0': newdata}
        print(mydict)
        pd_particle_events = pd.DataFrame(mydict)  # convert list of dict's into pandas dataframe
        print("Writing to CSV file :", output_file)
        pd_particle_events.to_csv(output_file, sep='|', header=True, index=False, chunksize=50000, mode='w', encoding='utf-8')
        exit(0)
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    axs[0].hist(initdata, histtype='step', bins=1200)
    axs[1].hist(newdata, histtype='step', bins=1200)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description='Geant4 Macro Scheduler')

    parser.add_argument('--root_file', dest='root_file_name', required=True,
                        help="Root file to read in")
    parser.add_argument('--output', dest='output_file', default=None, required=False,
                        help="Name of output CSV file")

    parser.set_defaults(configFile="g4macsched.cfg")

    args, unknown = parser.parse_known_args()
    # processConfigFile(args.configFile)
    process_root_file(args.root_file_name, args.output_file)


if __name__ == "__main__":
    main()
