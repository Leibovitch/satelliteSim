import json
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt 
from globalRegistry import GlobalRegistry
import orbital_mechanics
from mpl_toolkits.mplot3d import Axes3D
from orbit import Orbit

ureg = GlobalRegistry.getRegistry().ureg
constants = json.load(open("constants.json"))
R = constants['earth radius'] * ureg.meter
class Satellite():
    def __init__(self, res_at_nadir_at_500, look_angle, band, initial_orbit_params, initial_time, index=0, coordinate_system="kepler"):
        self.res_at_nadir_at_500 = res_at_nadir_at_500
        self.orbit = Orbit(initial_orbit_params, initial_time, coordinate_system)
        self.look_angle = look_angle
        self.band = band
    
    def calculate_range(self, look_angle):
        range_to_target = R * np.cos(look_angle) - np.sqrt(R^2*np.cos(look_angle)^2 -2)

    def get_central_angle(self, look_angle):
        h = self.orbit.get_height()
        r = h / R
        central_angle = - look_angle + np.arcsin((1 + r) * np.sin(look_angle))
        return central_angle

    def get_look_angle(self, res):
        if res <= self.res_at_nadir_at_500:
            pass
            # throw some staff
        else:
            nad_res = self.res_at_nadir_at_500
            h = self.orbit.get_height()
            range_to_target = res / nad_res * h
            look_angle = np.arccos((range_to_target ** 2 + 2 * h * R + h ** 2) / 2 / R / range_to_target)
        
        return look_angle

    def get_access_to_target(self, target, total_sim_time, sim_step):
        (lat, lon) = target.locations
        res = target.resolution 
        target_location = np.array([np.cos(lat.to('radians').magnitude) * np.cos(lon.to('radians').magnitude),
         np.cos(lat.to('radians').magnitude) * np.sin(lon.to('radians').magnitude),
          np.sin(lat.to('radians').magnitude)]) 
        central_angle = self.get_central_angle(self.get_look_angle(res))
        locations = self.orbit.get_orbit(total_sim_time, sim_step)
        angle = np.arccos(np.einsum('ij, i -> j',locations['r'].to('meters').magnitude / np.linalg.norm(locations['r'].to('meters').magnitude, axis=0), target_location))
        return  angle < central_angle

    def get_access_to_location(self, lon, lat, res, total_sim_time, sim_step):
        target_location = np.array([np.cos(lat.to('radians').magnitude) * np.cos(lon.to('radians').magnitude),
         np.cos(lat.to('radians').magnitude) * np.sin(lon.to('radians').magnitude),
          np.sin(lat.to('radians').magnitude)]) 
        central_angle = self.get_central_angle(self.get_look_angle(res))
        locations = self.orbit.get_orbit(total_sim_time, sim_step)
        angle = np.arccos(np.einsum('ij, i -> j',locations['r'].to('meters').magnitude / np.linalg.norm(locations['r'].to('meters').magnitude, axis=0), target_location))
        return  angle < central_angle

