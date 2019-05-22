import numpy as np
import time
import json
from numpy.linalg import norm as norm
from numpy.linalg import multi_dot
import matplotlib.pyplot as plt
from globalRegistry import GlobalRegistry
import orbital_mechanics
from mpl_toolkits.mplot3d import Axes3D 

ureg = GlobalRegistry.getRegistry().ureg
constants = json.load(open("constants.json"))
earths_radius = constants["earth radius"] * ureg.meters
mu = constants["mu"] * ureg.meter ** 3 / ureg.seconds ** 2 
earths_mass = constants["earth mass"] *  ureg.kilograms
j2 = constants["j2"]
earth_period = constants["earth period"] * ureg.seconds


class Orbit():
    def __init__(self, initial_location, initial_time, coordinate_system="cartesian"):
        self.initialize(initial_location, initial_time, coordinate_system)

    def initialize(self, initial_location, initial_time, coordinate_system="cartesian"):
        self.times = np.array([initial_time])
        if coordinate_system == "kepler":
            self.cartesian = orbital_mechanics.convert_kepler_to_cartesian(initial_location)
            self.kepler_params = initial_location
        else:
            self.cartesian = initial_location
            self.kepler_params = orbital_mechanics.convert_cartesian_to_kepler(initial_location)

        self.initial_orbital_corrections = self.get_orbital_corrections()
        self.T = self.get_period()
        self.n = 2 * np.pi / self.T
        if not hasattr(self, 'RAAN_dot'):
            (self.RAAN_dot, self.M_dot, self.w_dot) = self.get_orbital_corrections()
        elif(self.kepler_params.e != 0):
            (self.RAAN_dot, self.M_dot, self.w_dot) = self.get_orbital_corrections()
    
    def get_orbital_corrections(self):
        n = 2 * np.pi / self.get_period()
        kepler_params = self.kepler_params
        RAAN_dot = -1.5 * n * j2 * np.cos(kepler_params['i']) /((1 - kepler_params['e'] ** 2) ** 2) * (earths_radius / kepler_params['a'].to('meter')) ** 2 * ureg.radians
        M_dot = n + 3 * n * j2 * (3 * np.cos(kepler_params['i']) ** 2 - 1) / 4 / (( 1 - kepler_params['e'] ** 2) ** 3 / 2) * (earths_radius / kepler_params['a']) ** 2 *ureg.radians
        w_dot = -0.75 * n * j2 * (1 - 5 * np.cos(kepler_params['i']) ** 2) / (1 - kepler_params['e'] ** 2) ** -2 * (earths_radius / kepler_params['a']) ** 2 * ureg.radians
        return (RAAN_dot.to(ureg.radians / ureg.seconds), M_dot.to(ureg.radians / ureg.seconds), w_dot.to(ureg.radians / ureg.seconds))

    def get_location_after_time(self, time_delta, coordinate_system="cartesian", includ_earth_rotation="True"):
        # self.time = self.time + time_delta
        kepler_location = self.kepler_params
        if includ_earth_rotation:
            RAAN_dot = self.RAAN_dot + 2 * np.pi * ureg.radians / earth_period
        else:
            RAAN_dot = self.RAAN_dot

        kepler_location['RAAN'] = kepler_location['RAAN'] + RAAN_dot * time_delta.seconds * ureg.seconds
        kepler_location['w'] = kepler_location['w'] + self.w_dot * time_delta.seconds * ureg.seconds
        M = orbital_mechanics.get_mean_anomalie_from_true_anomalie(kepler_location['nu'], kepler_location['e'])
        kepler_location['nu'] = orbital_mechanics.approximate_eccentric_anomalie(M + self.M_dot * time_delta.seconds * ureg.seconds, kepler_location['e'])
        return orbital_mechanics.convert_kepler_to_cartesian(kepler_location)

    def get_orbit(self, total_time, steps, coordinate_system="cartesian", includ_earth_rotation="True"):
        self.times = np.arange(0, np.ceil(total_time / steps), 1) * steps.seconds
        kepler_location = self.kepler_params
        if includ_earth_rotation:
            RAAN_dot = self.RAAN_dot + 2 * np.pi * ureg.radians / earth_period
        else:
            RAAN_dot = self.RAAN_dot

        kepler_location['a'] = kepler_location['a'] * np.ones(len(self.times))
        kepler_location['e'] = kepler_location['e'] * np.ones(len(self.times))
        kepler_location['i'] = kepler_location['i'] * np.ones(len(self.times))
        kepler_location['RAAN'] = kepler_location['RAAN'] + RAAN_dot * self.times * ureg.seconds
        kepler_location['w'] = kepler_location['w'] + self.w_dot * self.times * ureg.seconds
        M = orbital_mechanics.get_mean_anomalie_from_true_anomalie(kepler_location['nu'], kepler_location['e'])
        kepler_location['nu'] = orbital_mechanics.approximate_eccentric_anomalie(M + self.M_dot * self.times * ureg.seconds, kepler_location['e'])
        self.locations = orbital_mechanics.convert_kepler_to_cartesian(kepler_location)
        return self.locations

    def plot_orbit(self):
        ax = plt.axes(projection='3d')
        ax.plot3D(self.locations['r'].magnitude[0,:], self.locations['r'].magnitude[1,:], self.locations['r'].magnitude[2,:])
        plt.show()

    def get_period(self):
        period = 2 * np.pi * np.sqrt(self.kepler_params['a'].to('meter') ** 3 / mu)
        return period

    def get_orbital_energy(self):
        r = self.get_radius()
        v = self.get_speed()
        specific_potential_energy = -mu / r
        specifuc_kinetic_energy = 1/2 * v ** 2
        specific_total_energy = specific_potential_energy + specifuc_kinetic_energy
        return specific_total_energy

    def get_orbital_angular_momentum(self):
        location = [self.cartesian['x'].to(ureg.meter).magnitude, self.cartesian['y'].to(ureg.meter).magnitude, self.cartesian['z'].to(ureg.meter).magnitude]
        velocity = [self.cartesian['vx'].to(ureg.meter / ureg.seconds).magnitude, self.cartesian['vy'].to(ureg.meter / ureg.seconds).magnitude, self.cartesian['vz'].to(ureg.meter / ureg.seconds).magnitude]
        return np.cross(location, velocity) * ureg.meter ** 2 / ureg.second

    def get_radius(self):
        location = [self.cartesian['x'].to(ureg.meter).magnitude, self.cartesian['y'].to(ureg.meter).magnitude, self.cartesian['z'].to(ureg.meter).magnitude]
        return norm(location) * ureg.meter

    def get_height(self):
        return norm(self.cartesian['r'], axis=0) * self.cartesian['r'].units - earths_radius

    def get_speed(self):
        return norm(self.cartesian['v'], axis=0) * ureg.meter / ureg.seconds

    @staticmethod
    def _get_speed(cartesian):
        return norm(cartesian['v'], axis=0) * ureg.meter / ureg.seconds
        
    @staticmethod
    def _get_orbital_energy(cartesian):
        r = Orbit._get_radius(cartesian)
        v = Orbit._get_speed(cartesian)
        specific_potential_energy = -mu / r
        specifuc_kinetic_energy = 1/2 * v ** 2
        specific_total_energy = specific_potential_energy + specifuc_kinetic_energy
        return specific_total_energy

    @staticmethod
    def _get_orbital_angular_momentum(cartesian):
        return np.cross(cartesian['r'].T, cartesian['v'].T).T * ureg.meter **2 / ureg.second

    @staticmethod
    def _get_radius(cartesian):
        return norm(cartesian['r']) * ureg.meter
