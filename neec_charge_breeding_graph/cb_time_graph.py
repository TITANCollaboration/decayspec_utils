import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from cbgraph import calc_neec_xsec

ions_pps = 1e5
q41_in_trap = 0

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

def optimize(poly_fit):
    #bnds = ((0, None), (0, None))
    est_time = 2
    results = opt.minimize(ions_in_trap_func, est_time, args=poly_fit)
    print(results.x)


def plot_cp_rate(time_lin_space, my_func_output, charge, poly_fit):
    plt.plot(time_lin_space, my_func_output, label="Polynomial Fit")
    plt.plot(np.linspace(0, 1, len(charge)), charge, label="EBIT Sim output")
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

def ions_in_trap_func(time, poly_fit):#, q41_in_trap):
    global q41_in_trap
    q41_in_trap = q41_in_trap + ions_pps * time * cb_rate_func(time % 1, poly_fit)
    return q41_in_trap

def carefully_mod_time(time):
    if (time % 1 == 0) and (time != 0):
        return 1
    else:
        return time % 1

def calc_ions_in_trap(poly_fit, new_ions_in_trap, injections_per_time=1, total_trapping_time=1):
    ions_in_trap = []
    if injections_per_time == 1:
        time_lin_space = [total_trapping_time]
    else:
        time_lin_space = np.linspace(0, total_trapping_time, injections_per_time)
    for mytime in time_lin_space:
        #(ions_pps/injections_per_time)
        cycle_time = carefully_mod_time(mytime)
        q41_in_trap = (new_ions_in_trap) * cb_rate_func(cycle_time, poly_fit)
        if injections_per_time != 1:
            for prev_time in time_lin_space[0:np.where(time_lin_space == mytime)[0][0]]:
                prev_cycle_time = carefully_mod_time(mytime - prev_time)
                q41_in_trap = q41_in_trap + (new_ions_in_trap) * cb_rate_func(prev_cycle_time, poly_fit)
        ions_in_trap.append(q41_in_trap)
    return max(ions_in_trap)
        #print(str("{:.1e}".format(q41_in_trap)))
#    plot_ions_in_trap(time_lin_space, ions_in_trap)
#    return ions_in_trap


def calc_over_breeding_times(poly_fit, breeding_time_per_cycle=1.0):
    ions_over_injection_times = []
    ions_over_cycle_time = []
    neec_event_count = []
    neec_count = 0
    #breeding_time_per_cycle = 1.0
    trapping_time_per_cycle = 10
    breeding_time_lin_space = np.linspace(0, breeding_time_per_cycle, 110)
    trapping_time_lin_space = np.linspace(1, trapping_time_per_cycle, trapping_time_per_cycle)
#calc_neec_xsec(#ions)
    #  calc_neec_xsec(ions_in_trap)
    ion_count = 0
    new_ions_in_trap = 0
    for trapping_time in trapping_time_lin_space:
        new_ions_in_trap = new_ions_in_trap + ions_pps
        if (trapping_time % 2) == 0:  # We're in Breeding time
            for breeding_time in breeding_time_lin_space:
                ions_over_injection_times.append(calc_ions_in_trap(poly_fit, new_ions_in_trap, 1, breeding_time))
            ion_count = ion_count + max(ions_over_injection_times)
            ions_over_cycle_time.append(ion_count)
            new_ions_in_trap = 0
        else:  # We're in NEEC time
            neec_time_this_cycle = 2 - breeding_time_per_cycle
            neec_event_count.append(calc_neec_xsec(ion_count) * neec_time_this_cycle)
    #print("Possible NEEC Events: ", sum(neec_event_count))
    return sum(neec_event_count)
    #print("NEEC Time!", neec_time)
    #print("NEEC Events per second:", calc_neec_xsec())
    print(ions_over_cycle_time)
    plt.plot(np.linspace(0, trapping_time_per_cycle, len(ions_over_cycle_time)), ions_over_cycle_time)
#    #plt.xlim([1, 100])
    plt.title("Ion Injection timing", fontsize=20)
    plt.xlabel("Breeding time - seconds", fontsize=20)
    plt.ylabel("# 41+ ions in trap", fontsize=20)
    plt.show()

def graph_over_breeding_times(poly_fit):
    neec_rate = []
    breeding_time_range = np.linspace(0, 1, 101)
    for breeding_time_per_cycle in breeding_time_range:
        neec_rate.append(calc_over_breeding_times(poly_fit, breeding_time_per_cycle))
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
    p = np.poly1d(poly_fit)
    return p(time)

def calc_cb_rate(poly_fit, charge):
    charge_breed_eq = []
    time_lin_space = np.linspace(0, 1, 100)
    for mytime in time_lin_space:
        my_val = 0
        my_val = cb_rate_func(mytime, poly_fit)
        charge_breed_eq.append(my_val)
    plot_cp_rate(time_lin_space, charge_breed_eq, charge, poly_fit)
    print(eq_text_from_fit(poly_fit))
    return charge_breed_eq

def fit_me(order):
    time = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0]
    charge = [0, 0, 0.024, 0.178, 0.408, 0.598, 0.717, 0.78, 0.81, 0.83, 0.84]
    poly_fit = np.polyfit(time, charge, order)
    print(poly_fit)
    return poly_fit, charge



poly_fit, charge = fit_me(4)
#calc_over_injection_times(poly_fit)
#calc_over_breeding_times(poly_fit)
graph_over_breeding_times(poly_fit)
#calc_ions_in_trap (poly_fit)
#calc_cb_rate(poly_fit, charge)
