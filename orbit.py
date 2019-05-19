import numpy as np
import time
import json
from numpy.linalg import norm as norm
from numpy.linalg import multi_dot
import matplotlib.pyplot as plt
from globalRegistry import GlobalRegistry
import orbital_mechanics

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
        self.time = initial_time
        if coordinate_system == "kepler":
            self.cartesian = orbital_mechanics.convert_kepler_to_cartesian(initial_location)
            self.kepler_params = initial_location
        else:
            self.cartesian = initial_location
            self.kepler_params = orbital_mechanics.convert_cartesian_to_kepler(initial_location)

        self.cartesian_location = [self.cartesian['x'].to(ureg.meter).magnitude, self.cartesian['y'].to(ureg.meter).magnitude, self.cartesian['z'].to(ureg.meter).magnitude] * ureg.meter
        self.cartesian_velocity = [self.cartesian['vx'].to(ureg.meter /ureg.seconds).magnitude , self.cartesian['vy'].to(ureg.meter /ureg.seconds).magnitude, self.cartesian['vz'].to(ureg.meter /ureg.seconds).magnitude] * ureg.meter / ureg.seconds
        self.initial_orbital_corrections = self.get_orbital_corrections()
        self.T = self.get_period()
        self.n = 2 * np.pi / self.T
        if not hasattr(self, 'RAAN_dot'):
            (self.RAAN_dot, self.M_dot, self.w_dot) = self.get_orbital_corrections()
        elif(self.kepler_location.e != 0):
            (self.RAAN_dot, self.M_dot, self.w_dot) = self.get_orbital_corrections()
    
    def get_orbital_corrections(self):
        n = 2 * np.pi / self.get_period()
        kepler_params = orbital_mechanics.convert_cartesian_to_kepler(self.cartesian)
        RAAN_dot = -1.5 * n * j2 * np.cos(kepler_params['i']) /(1 - kepler_params['e'] ** 2) ** 2 * (earths_radius / kepler_params['a']) ** 2 * ureg.radians
        M_dot = n + 3 * n * j2 * (3 * np.cos(kepler_params['i']) ** 2 - 1) / 4 / ( 1 - kepler_params['e'] ** 2) ** 3 / 2 * (earths_radius / kepler_params['a']) ** 2 *ureg.radians
        w_dot = -0.75 * n * j2 *(1 - 5 * np.cos(kepler_params['i']) ** 2) / (1 - kepler_params['e'] ** 2) ** -2 * (earths_radius / kepler_params['a']) ** 2 * ureg.radians
        return (RAAN_dot, M_dot, w_dot)

    def get_location_after_time(self, time_delta, coordinate_system="cartesian", includ_earth_rotation="True"):
        next_time = self.time + time_delta
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
    
    def get_orbit(self, total_time, steps):
        locations = []
        for i in range(0, int(np.ceil(total_time / steps))):
            locations.append(self.get_location_after_time(i * steps)) 

        return locations


    def get_period(self):
        kepler = orbital_mechanics.convert_cartesian_to_kepler(self.cartesian)
        period = 2 * np.pi * np.sqrt(kepler['a'] ** 3 / mu)
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
       return norm(self.cartesian_location) * self.cartesian_location.units - earths_radius

    @staticmethod
    def _get_speed(cartesian):
        velocity = [cartesian['vx'].to(ureg.meter / ureg.seconds).magnitude, cartesian['vy'].to(ureg.meter / ureg.seconds).magnitude, cartesian['vz'].to(ureg.meter / ureg.seconds).magnitude]
        return norm(velocity) * ureg.meter / ureg.seconds
        
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
        location = [cartesian['x'].to(ureg.meter).magnitude, cartesian['y'].to(ureg.meter).magnitude, cartesian['z'].to(ureg.meter).magnitude]
        velocity = [cartesian['vx'].to(ureg.meter / ureg.second).magnitude, cartesian['vy'].to(ureg.meter / ureg.second).magnitude, cartesian['vz'].to(ureg.meter / ureg.second).magnitude]
        return np.cross(location, velocity) * ureg.meter **2 / ureg.second

    @staticmethod
    def _get_radius(cartesian):
        location = [cartesian['x'].to(ureg.meter).magnitude, cartesian['y'].to(ureg.meter).magnitude, cartesian['z'].to(ureg.meter).magnitude]
        return norm(location) * ureg.meter





#     def propogate_to_next_time(new_time):
#         pass

#     def propogate_by_time_delta(time_delta):
#         get_location_after_time()

# def get_radius():
#     location = [cartesian.location.x, cartesian.location.y, cartesian.location.z]
#     return norm(location)



# def get_orbital_energy(cartesian):
#     r = get_radius(cartesian)
#     v = get_speed(cartesian)
#     specific_potential_energy = -mu / r
#     specifuc_kinetic_energy = 1/2 * v ** 2
#     specific_total_energy = specific_potential_energy + specifuc_kinetic_energy

# def get_orbital_angular_momentum(cartesian):
#     location = [cartesian.location.x, cartesian.location.y, cartesian.location.z]
#     velocity = [cartesian.location.vx, cartesian.location.vy, cartesian.location.vz]
#     return np.cross(location, velocity)