import math

def calc_neec_xsec(ionDensity=2e5):
    # This comes out to per second
    #cross_sec_neec_neon = 3.0e-5  # Theoretical

    cross_sec_neec_neon = 6.7e-6  # Updated
    __C__ = 3.0e10
    __EMASS__ = 5.11e5
    __ECHG__ = 1.6e-19
    beamEnergy = 7000
    beamCurrent = 0.2
    beamRadius = 1.2e-2
    barn = 1e-28 # m^2
    barn_cm = barn * 1e4 # cm^2
    volume = math.pi * beamRadius**2 * 6
    resonance_energy = .01
    currentDensity = beamCurrent / (math.pi*beamRadius**2)
    #print("Current Density: ", currentDensity)
    electronVelocity = __C__*math.sqrt(1-(beamEnergy/__EMASS__+1)**-2)
    #print("Electron Velocity:  ", electronVelocity)
    electronDensity = currentDensity / __ECHG__ / electronVelocity
    #print("Electron DENSITY:  ", electronDensity)
    calc_xsec = volume * ionDensity * electronDensity * electronVelocity * cross_sec_neec_neon * barn_cm * resonance_energy
    #print("Cross sec:  ", calc_xsec)
    #print("Electron Velocity (cm/s)", f"{electronVelocity:.2e}")
    #print("Electron Density cm^-3:", f"{electronDensity:.2e}")
    return calc_xsec  # per second


def simple_neec(total_trap_particle_count, neecing_time, breeding_time):
    array_eff = 0.001
    order_increase_list = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9]
    days = 7  # days to run
    neec_total = 0
    cb_avg = 0.41
    trap_neon_like_count = total_trap_particle_count * cb_avg
    neec_xsec_w_cb_avg = calc_neec_xsec(trap_neon_like_count)
    neec_total_count_per_neec_cycle = neec_xsec_w_cb_avg * neecing_time
    cycles_per_hour = (3600/(neecing_time + breeding_time))
    neec_total = neec_total_count_per_neec_cycle * cycles_per_hour

    print("Neon-like ions in Trap:", trap_neon_like_count)
    print("XSec for Neon-like ions in trap eV/s:", f"{neec_xsec_w_cb_avg:.2e}")
    print("Percent Avg Charge bread to Neon-like:", cb_avg)
    print("Breeding Time:", breeding_time, "NEECing Time:", neecing_time)
    print("Total NEEC counts per cycle:", neec_total_count_per_neec_cycle)
    print("Cycles per hour (3600/(neecing_time + breeding_time)):", cycles_per_hour)
    print("NEEC Counts per hour:", neec_total)
    print("NEEC count over", days, "days:", neec_total*24*days)
    print("NEEC count over", days, "days w/ Array Eff:", neec_total*24*days*array_eff)
    for my_order in order_increase_list:
        print("Order Multiplier:", f"{my_order:.0e}", "Observed NEEC evets(w/array eff):", f"{neec_total*24*days*array_eff*my_order:.2e}")

print("----------MRTOF TO EBIT----------")
simple_neec(total_trap_particle_count=2e5, neecing_time=60, breeding_time=.7)
print("----------RFQ TO EBIT----------")
simple_neec(total_trap_particle_count=1e7, neecing_time=60, breeding_time=10)
