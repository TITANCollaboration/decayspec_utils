import uproot
import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse


def read_root_files(root_file_path):
    energy_vals = []
    eff_vals = []
    for file_name in glob.glob("%s/*.root" % root_file_path):
        print(file_name)
        energy = int((file_name.split('/')[-1]).split('.root')[0])
        print("Checking Energy : %i keV" % energy)
        events = uproot.open(file_name)['EVENT_NTUPLE']
        rounded_pulse_height = np.around(events.array("pulse_height"))
        #print(rounded_pulse_height)
        occurances = np.count_nonzero((rounded_pulse_height >= (energy - (energy*.01))) & (rounded_pulse_height <= (energy + (energy*.01))))
        energy_vals.append(energy)
        eff_vals.append((occurances/1000000) * 100)
        #print(eff_dict)
    #plt.xscale('log')
    plt.xlim([10, 2000])
    #plt.xticks(np.arange(0, 3200, step=100))
    plt.grid(color='k', linestyle='-', linewidth=.2)
    plt.title("EBIT 8pi array Efficiency", fontsize=20)
    plt.xlabel("Energy (keV)", fontsize=15)
    plt.ylabel("Efficiency (%)", fontsize=15)
    #plt.xticks(energy_vals)
    plt.tick_params(labelsize=15)

    plt.scatter(energy_vals, eff_vals)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description='EBIT 8Pi Efficiency grapher')

    parser.add_argument('--root_file_path', dest='root_file_path', required=True,
                        help="Directory of root files")
    args, unknown = parser.parse_known_args()

    read_root_files(args.root_file_path)
    return 0


if __name__ == "__main__":
    main()
