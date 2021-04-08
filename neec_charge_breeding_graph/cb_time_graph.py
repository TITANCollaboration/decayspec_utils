import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from cbgraph import calc_neec_xsec
import math
from scipy.integrate import quad

ions_pps = 2e5
q41_in_trap = 0
sb_129m1_halflife = 1062
sb_129_halflife = 16560

def eq_text_from_fit(poly_fit):
    my_fit_eq = ""
    poly_degree_index = len(poly_fit) - 1
    for poly_fit_term in poly_fit:
        if poly_fit_term > 0:
            my_fit_term = '+' + str(round(poly_fit_term, 3))
        else:
            my_fit_term = str(round(poly_fit_term, 3))
        my_fit_eq = my_fit_eq + my_fit_term
        if poly_degree_index > 1:
            my_fit_eq = my_fit_eq + 'x^' + str(poly_degree_index)
        elif poly_degree_index == 1:
            my_fit_eq = my_fit_eq + 'x'
        poly_degree_index = poly_degree_index - 1
    print("Terms in fit.. ", my_fit_eq)
    return my_fit_eq


def plot_cb_rate(time_lin_space, my_func_output):  #, poly_fit):
    plt.plot(time_lin_space, my_func_output, label="Polynomial Fit")
    #plt.plot(np.linspace(0, 1, len(charge)), charge, label="EBIT Sim output")
    plt.title("Fitting Sb129 +41 charge breeding", fontsize=20)
    plt.xlabel("Seconds", fontsize=20)
    plt.ylabel("% Charge state in trap", fontsize=20)
    plt.text(0.4, .4, eq_text_from_fit(poly_fit))

    plt.legend()
    plt.show()

def plot_ions_in_trap(time, ions):
    my_title = "Ions in Trap @ " + str(ions_pps) + " PPS"
    plt.plot(time, ions)
    plt.title(my_title, fontsize=20)
    plt.xlabel("Seconds", fontsize=20)
    plt.ylabel("# of q41+ ions in EBIT trap", fontsize=20)
    plt.show()
    return

def nuclear_activity(time_in_trap, half_life, particle_count):
    # Integrate over the time in the trap to get decays per trapping time
    lamb = math.log(2) / half_life
    my_activity = lambda time: (lamb * particle_count * math.exp(-lamb*time))
    decay_events = quad(my_activity, 0, time_in_trap)
    return decay_events[0]

def ions_in_trap_func(time, poly_fit):  # , q41_in_trap):
    global q41_in_trap
    q41_in_trap = q41_in_trap + ions_pps * time * cb_rate_func(time % 1, poly_fit)
    return q41_in_trap

def carefully_mod_time(time):
    if (time % 1 == 0) and (time != 0):
        return 1
    else:
        return time % 1

def calc_ions_in_trap(poly_fit, new_ions_in_trap, injections_per_time, total_breeding_time):
    # Usually used per breeding cycle
    ions_in_trap = []
    if injections_per_time == 1:
        time_lin_space = [total_breeding_time]
    else:
        time_lin_space = np.linspace(0, total_breeding_time, injections_per_time)
    for mytime in time_lin_space:
        cycle_time = carefully_mod_time(mytime)
        q41_in_trap = (new_ions_in_trap) * cb_rate_func(cycle_time, poly_fit)
        if injections_per_time != 1:
            for prev_time in time_lin_space[0:np.where(time_lin_space == mytime)[0][0]]:
                prev_cycle_time = carefully_mod_time(mytime - prev_time)
                q41_in_trap = q41_in_trap + ((new_ions_in_trap) * cb_rate_func(prev_cycle_time, poly_fit))
        ions_in_trap.append(q41_in_trap)
    return max(ions_in_trap)


def calc_over_breeding_times(poly_fit, new_ions_in_trap, breeding_time_per_cycle=1.0):
    ions_over_injection_times = []
    breeding_time_lin_space = np.linspace(0, breeding_time_per_cycle, 11)
    for breeding_time in breeding_time_lin_space:
        ions_over_injection_times.append(calc_ions_in_trap(poly_fit, new_ions_in_trap, 1, breeding_time))
    return ions_over_injection_times

def calc_over_trapping_time(poly_fit, trapping_time_per_cycle, breeding_time_per_cycle):
    # This function handles going back and forth between Breeding and NEECing
    new_ions_in_trap = 0
    cb_ions_in_trap = 0
    decay_events = 0
    trapping_time_lin_space = np.linspace(1, trapping_time_per_cycle, 10) #, dtype=int)  # This gives us [1,2, ...trapping_time_per_cycle]
    ions_over_breeding_time = []
    neec_event_count = []
    total_ions_in_trap = 0
    for trapping_time in trapping_time_lin_space:
        new_ions_in_trap = new_ions_in_trap + ions_pps
        total_ions_in_trap = total_ions_in_trap + new_ions_in_trap
        decay_events = decay_events + nuclear_activity(1, sb_129m1_halflife, total_ions_in_trap)

        if (trapping_time % 2) != 0:  # We're in Breeding time
            #print("Breeding!")
            ions_over_breeding_time.append(calc_over_breeding_times(poly_fit, new_ions_in_trap, breeding_time_per_cycle=breeding_time_per_cycle))
            cb_ions_in_trap = cb_ions_in_trap + ions_over_breeding_time[-1][-1]
            new_ions_in_trap = 0
        else:  # We're in NEEC time
            print("Neec'ing!")
            neec_time_this_cycle = 2 - breeding_time_per_cycle
            neec_event_count.append(calc_neec_xsec(cb_ions_in_trap) * neec_time_this_cycle)
    return sum(neec_event_count), total_ions_in_trap, decay_events


def calc_over_long_neec_time_short_breed_time(trapping_time_per_cycle, breeding_time, avg_charge_pop):
    # This function handles long NEEC'ing times with one short charge breeding time
    # Note this simplifies a LOT Of stuff, only use this one when your breeding time is a small fraction of your NEEC'ing time
    # on the order of 1/60 where we do not need to care about the breeding time specifically and can average the charge population
    global ions_pps
    decay_events = 0
    neec_event_count = 0
    total_ions_in_trap = 0
    total_ions_in_trap = ions_pps * trapping_time_per_cycle
    total_cb_ions_in_trap = total_ions_in_trap * avg_charge_pop
    print("Total ions in TRAP", total_ions_in_trap)
    decay_events = decay_events + nuclear_activity(trapping_time_per_cycle + breeding_time, sb_129m1_halflife, total_ions_in_trap)
    neec_event_count = calc_neec_xsec(total_cb_ions_in_trap) * trapping_time_per_cycle

    return neec_event_count, total_cb_ions_in_trap, decay_events

def plot_long_cb_rate_over_hour(trapping_time_per_cycle, breeding_time, avg_charge_pop):
    one_hr_lin_space = np.linspace(0, 60, 61, dtype=int)
    decays_over_hour = []
    neecs_over_hour = []
    total_neec_events = 0
    total_decay_events = 0
    for my_min in one_hr_lin_space:
        neec_events, total_cb_ions_in_trap, decay_events = calc_over_long_neec_time_short_breed_time(trapping_time_per_cycle, breeding_time, avg_charge_pop)
        total_neec_events = total_neec_events + neec_events
        total_decay_events = total_decay_events + decay_events
        neecs_over_hour.append(total_neec_events)
        decays_over_hour.append(total_decay_events)
    print(decays_over_hour)
    #plt.plot(one_hr_lin_space, neecs_over_hour)
    plt.plot(one_hr_lin_space, decays_over_hour)

    #plt.title("NEEC events per trapping cycle (10s)", fontsize=20)
    plt.xlabel("Time - minutes (60sec trapping time)", fontsize=20)
    plt.ylabel("Decay events", fontsize=20)
    #plt.ylabel("# of NEEC Events Possible", fontsize=20)
    plt.show()

def graph_over_trapping_times(poly_fit):
    global sb_129m1_halflife
    max_trapping_time = 200
    decay_events_per_trapping_times = []
    neec_events_per_trapping_times = []
    decay_events_per_trapping_times_for_hour = []
    decays_per_trapping_times = []
    neec_events_per_trapping_times_for_hour = []
    trapping_time_range = np.linspace(0, max_trapping_time,max_trapping_time + 1, dtype=int)
    for my_trapping_time in trapping_time_range:
        #print(my_trapping_time)
        neec_events, ions_in_trap, decay_events = calc_over_trapping_time(poly_fit, trapping_time_per_cycle=my_trapping_time)

        neec_events_per_trapping_times.append(neec_events)
        neec_events_per_trapping_times_for_hour.append((neec_events * (3600/my_trapping_time)))

        decay_events_per_trapping_times.append(decay_events)
        decay_events_per_trapping_times_for_hour.append((decay_events * (3600/my_trapping_time)))

        print("Trapping time:", my_trapping_time, "Neec Events:", neec_events)

#    for breeding_time_per_cycle in breeding_time_range:
    decay_events_per_trapping_times_for_hour_np = np.array(decay_events_per_trapping_times_for_hour)
    neec_events_per_trapping_times_for_hour_np = np.array(neec_events_per_trapping_times_for_hour)
    signal_to_noise = neec_events_per_trapping_times_for_hour_np / decay_events_per_trapping_times_for_hour_np
    print("Decay", decay_events_per_trapping_times_for_hour_np)
    print("Neec", neec_events_per_trapping_times_for_hour_np)
    print("SNR", signal_to_noise)
    #plt.plot(trapping_time_range, neec_events_per_trapping_times_for_hour)

    #plt.plot(trapping_time_range, neec_events_per_trapping_times)
    #plt.plot(trapping_time_range, decay_events_per_trapping_times_for_hour)
    plt.plot(trapping_time_range, signal_to_noise)
    #plt.title("NEEC events per trapping cycle (10s)", fontsize=20)
    plt.xlabel("Trapping time - seconds", fontsize=20)
    #plt.ylabel("Decay events per hour", fontsize=20)
    plt.ylabel("Signal to noise over hour", fontsize=20)

    #plt.ylabel("Total # of NEEC events per hour", fontsize=20)
    plt.show()


def graph_over_breeding_times(poly_fit):
    neec_rate = []
    breeding_time_range = np.linspace(0, 1, 110)
    for breeding_time_per_cycle in breeding_time_range:
        neec_rate.append(calc_over_trapping_time(poly_fit, trapping_time_per_cycle=10, breeding_time_per_cycle=breeding_time_per_cycle)[0])
    print(neec_rate)
    print("NEEC Index : ", neec_rate.index(max(neec_rate)))

    print("Time index of MAX NEEC:", breeding_time_range[neec_rate.index(max(neec_rate))])
    plt.plot(breeding_time_range, neec_rate)
#    #plt.xlim([1, 100])
    plt.title("NEEC events per trapping cycle (10s)", fontsize=20)
    plt.xlabel("Breeding time - seconds", fontsize=20)
    plt.ylabel("# NEEC events", fontsize=20)
    plt.show()

def calc_over_injection_times(poly_fit):
    injection_range = 100
    ions_over_injection_times = []

    for injection_division in range(2, injection_range + 1, 1):
        ions_over_injection_times.append(calc_ions_in_trap(poly_fit, injection_division))

    print("Max :", max(ions_over_injection_times))
    plt.plot(np.linspace(1, 100, len(ions_over_injection_times)), ions_over_injection_times)
    #plt.xlim([1, 100])
    plt.title("Ion Injection timing", fontsize=20)
    plt.xlabel("# of ion injections per breeding cycle", fontsize=20)
    plt.ylabel("# 41+ ions in trap", fontsize=20)
    plt.show()

def cb_rate_func(time, poly_fit):
    #p = np.poly1d([0.0000007, 0.00183744, 1.1])  # Griffin smearing
    p = np.poly1d(poly_fit)
    return p(time)

def calc_cb_rate(poly_fit, charge):
    charge_breed_eq = []
    time_lin_space = np.linspace(0, 1, 100)
    for mytime in time_lin_space:
        my_val = 0
        my_val = cb_rate_func(mytime, poly_fit)
        charge_breed_eq.append(my_val)
    plot_cb_rate(time_lin_space, charge_breed_eq, charge, poly_fit)
    print(eq_text_from_fit(poly_fit))
    return charge_breed_eq

def fit_me(time, charge, order):
    #time = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0]
    #charge = [0, 0, 0.024, 0.178, 0.408, 0.598, 0.717, 0.78, 0.81, 0.83, 0.84]

    #time_decay = [0, .1, .2, .3, .4, .5, .6, .7]
    #charge_decay = [.78, .62, .38, .2, .16, .15]
    poly_fit = np.polyfit(time, charge, order)

    poly_fit2 = np.polyfit(time, charge, order)
    print(poly_fit)
    return poly_fit, charge

def decays_over_hour(trapping_time):
    trappings_per_hour = 3600 / trapping_time
    #trappings_per_day = trappings_per_hour * 24
    return trappings_per_hour
#  For long NEEC'ing times we reach a charge breed equilibrium of ~67%

def simple_neec(time, charge):
    #time_decay =   [0, .1, .2, .3, .4, .5, .6, .7]
    #charge_decay = [.78, .62, .41, .35, .2, .18, .15, .13]
    my_decay_fit = fit_me(time, charge, 4)[0]
    #print(my_decay_fit)
    print(cb_rate_func(60, my_decay_fit))
    neec_total = 0
    neec_time = 60
    breed_time = .7
    time_slices = 60
    trap_particle_count = 2e5
    time_slice_length = neec_time / time_slices
    for mytime in np.linspace(0, neec_time, time_slices):
        print("Time",mytime, "CB:", cb_rate_func(mytime, my_decay_fit))
        print(trap_particle_count * cb_rate_func(mytime, my_decay_fit))
        neec_total = neec_total + (calc_neec_xsec(trap_particle_count * cb_rate_func(mytime, my_decay_fit)) * (neec_time/time_slices))
    print("Neec total:", neec_total)
    neec_total = neec_total * (3600/(neec_time + breed_time)) * neec_time # 1285.71 #* .001
    print("1hr neec", neec_total)

#time = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0]
#charge = [0, 0, 0.024, 0.178, 0.408, 0.598, 0.717, 0.78, 0.81, 0.83, 0.84]

#time = np.linspace(0, 60, 600)
time = [0, .1, 1, 1.2, 1.4, 1.6, 2.0, 4.0, 6.0, 8.0, 10, 20, 40, 60]
#charge_decay = [.78, .62, .38, .2, .16, .15]  #  This is for 8.8keV
charge = [.78, .75, .72, .67, .64, .60, .54, .41, .41, .41, .41, .41, .41, .41]  # this is for 7.1keV

#poly_fit = [ 2.06312766e-11, -8.49120473e-09,  1.10899128e-06, -6.57615944e-05,
 # 1.91855075e-03, -2.60331334e-02,  5.32627725e-01]
poly_fit = fit_me(time, charge, 4)[0]
cb_rate_func(0, poly_fit)

print(poly_fit)
my_func_output = []
for mytime in time:
    my_val = 0
    my_val = cb_rate_func(mytime, poly_fit)
    my_func_output.append(my_val)
    print(my_val)
plot_cb_rate(time, my_func_output)

simple_neec(time, charge)
#poly_fit = np.array([-7.31702769e-10,  1.41824010e-07, -1.05850904e-05,  3.79854445e-04, -6.63719139e-03,  4.95546883e-02,  5.71005542e-01])
#print(calc_over_trapping_time(poly_fit, trapping_time_per_cycle=1.4, breeding_time_per_cycle=0.69))

#calc_cb_rate(poly_fit, charge=1)
#print(calc_over_long_neec_time_short_breed_time())
#plot_long_cb_rate_over_hour()
#graph_over_trapping_times(poly_fit)  # Find optimal trapping time
#graph_over_breeding_times(poly_fit)  # Find optimal breeding time
#trapping_time = 60
#decays_per_cycle = nuclear_activity(trapping_time, sb_129m1_halflife, ions_pps)
#print(decays_over_hour(trapping_time) * decays_per_cycle)

#decays_per_cycle = nuclear_activity(trapping_time, sb_129_halflife, decays_per_cycle)
#print(decays_over_hour(trapping_time) * decays_per_cycle)
