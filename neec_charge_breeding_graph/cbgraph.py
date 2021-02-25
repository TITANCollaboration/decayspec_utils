# Simple math to try and figure out how many particles of sb129m1 could make it to the EBIT and how many NEEC
# reactions can occur based on cross section.
import matplotlib.pyplot as plt
import numpy as np
import argparse

font_size = 'large'
cross_sec_neec_neon = 3.0e-5
#particle_count_observed_in_MRTOF = 11
attenuation = 10
iglis_supression_factor = 1.0e6
iglis_yield_reduction = 1/50
rfq_eff = .1
mrtof_eff = .3
ebit_injection_min = .05
ebit_injection_max = .3
seconds_in_hour = 3600
array_eff = 0.001
max_charge_breed_population = 0.8
max_charge_breed_population_avg = 0.75
est_PPS_in_ebit = [1.0e4, 1.0e5, 2.0e5, 3.0e5, 4.0e5, 1.0e6]  # This is at 100Hz
best_est_PPS_in_ebit = 1.0e5

def plot_stype(plt, ax):
    ax.grid()
    plt.figtext(0.15, 0.5, "Array Eff (~1MeV): " + str(array_eff*100) + "%", fontsize=font_size)
    plt.figtext(0.15, 0.45, "Neon-like Charge state in Trap (avg): " + str(max_charge_breed_population_avg*100) + "%", fontsize=font_size)
    plt.figtext(0.15, 0.4, "NEEC X-Sec: " + str("{:.1e}".format(cross_sec_neec_neon)), fontsize=font_size)

    plt.yticks(np.arange(0, 140, 5.0))
    plt.title("NEEC events over time - EGun energy tuned to NEEC : All the time, No charge breeding time")
    plt.xlabel("Time - seconds")
    plt.ylabel("# of likely NEEC events detected")
    plt.legend(loc='best')


def line_style(number):
    if number % 2:
        ls = '-'
    elif number % 3:
        ls = '--'
    elif number % 5:
        ls = ':'
    else:
        ls = '-.'
    return ls


def plot_neec_sum(args):
    neec_events_per_time = []
    cycle_time = 10
    t = np.arange(0, args.time/cycle_time, 1)

    PPS_in_ebit = [9e4, 1e5, 2e5, 3e5, 4e5]
    neec_sum_over_time = []
    for my_pps_in_ebit in PPS_in_ebit:
        particle_count_in_trap = 0
        neec_sum = 0
        for my_time in range(1, cycle_time + 1, 1):
            particle_count_in_trap = my_pps_in_ebit + (my_pps_in_ebit * my_time)
            if (my_time % 2) == 0:
                neec_sum = neec_sum + (particle_count_in_trap*array_eff*max_charge_breed_population_avg*cross_sec_neec_neon)
        neec_sum_over_time.append(neec_sum * t)

    fig, ax = plt.subplots()
    for my_pps in range(0, len(neec_sum_over_time), 1):
        ls = line_style(my_pps)
        ax.plot(t, neec_sum_over_time[my_pps], label="PPS to EBIT: " + str("{:.1e}".format(PPS_in_ebit[my_pps])), ls=ls)
    plot_stype(plt, ax)
    plt.title("NEEC events over Cycles - 10sec Cycle - EGun energy flipping between NEEC and CB each second")
    plt.xlabel("Cycles : 1 cycle = 10 seconds")
    plt.show()
    return 0


def plot_neec_rough(args):
    neec_events_per_time = []
    t = np.arange(0.0, args.time, 0.01)
    for my_pps in range(0, len(est_PPS_in_ebit), 1):
        neec_events_per_time.append(est_PPS_in_ebit[my_pps]*array_eff*max_charge_breed_population*cross_sec_neec_neon * t)
    fig, ax = plt.subplots()
    for my_pps in range(0, len(est_PPS_in_ebit), 1):
        ls = line_style(my_pps)
        ax.plot(t, neec_events_per_time[my_pps], label="PPS to EBIT: " + str("{:.1e}".format(est_PPS_in_ebit[my_pps])), ls=ls)
    plot_stype(plt, ax)
    plt.show()

    return 0


def main():
    parser = argparse.ArgumentParser(description='Geant4 Macro Scheduler')

    parser.add_argument('--time', dest='time', type=float, required=True,
                        help="Time in seconds")

    args, unknown = parser.parse_known_args()
    # processConfigFile(args.configFile)
    plot_neec_sum(args)
    # plot_neec_rough(args)

    return 0


if __name__ == "__main__":
    main()
