import numpy as np
import time
import json
from numpy.linalg import norm as norm
from numpy.linalg import multi_dot
import matplotlib.pyplot as plt
from globalRegistry import GlobalRegistry

ureg = GlobalRegistry.getRegistry().ureg
constants = json.load(open("constants.json"))
earths_radius = constants["earth radius"] * ureg.meters
mu = constants["mu"] * ureg.meter ** 3 / ureg.seconds ** 2
earths_mass = constants["earth mass"] *  ureg.kilograms
j2 = constants["j2"]
earth_period = constants["earth period"]
        
#rotations
def get_x_rotation_matrix(amount):
    x_rotation_matrix = np.array(
        (
            (1, 0, 0),
            (0 ,np.cos(amount), -np.sin(amount)),
            (0, np.sin(amount), np.cos(amount))
        )
    )
    return x_rotation_matrix

def get_y_rotation_matrix(amount):
    y_rotation_matrix = np.array(
        (
            (np.cos(amount), 0, np.sin(amount)),
            (0, 1, 0),
            (-np.sin(amount),0 ,np.cos(amount))
        )
    )
    return y_rotation_matrix

def get_z_rotation_matrix(amount):
    # print(amount)
    z_rotation_matrix = np.array(
        (
            (np.cos(amount), -np.sin(amount), 0),
            (np.sin(amount), np.cos(amount, 0), 0),
            (0, 0, 1)
        )
    )
    return z_rotation_matrix

def rotate_about_x(vector, amount):
    x_rotation_matrix = get_x_rotation_matrix(amount)
    rotated_vector = np.matmul(x_rotation_matrix, vector)
    return rotated_vector

def rotate_about_y(vector, amount):
    y_rotation_matrix = get_y_rotation_matrix(amount)
    rotated_vector = np.matmul(y_rotation_matrix, vector)
    return rotated_vector

def rotate_about_z(vector, amount):
    z_rotation_matrix = get_z_rotation_matrix(amount)
    rotated_vector = np.matmul(z_rotation_matrix, vector)
    return rotated_vector

def convert_cartesian_to_kepler(cartesian):
    velocity = [cartesian['vx'].to(ureg.meter / ureg.second).magnitude, cartesian['vy'].to(ureg.meter / ureg.second).magnitude, cartesian['vz'].to(ureg.meter / ureg.second).magnitude] * ureg.meter / ureg.second
    location = [cartesian['x'].to(ureg.meter).magnitude, cartesian['y'].to(ureg.meter).magnitude, cartesian['z'].to(ureg.meter).magnitude] * ureg.meter
    from orbit import Orbit
    specific_total_energy = Orbit._get_orbital_energy(cartesian)
    h = Orbit._get_orbital_angular_momentum(cartesian) #reduced_angular_momentum
    k = [0, 0, 1]
    n = np.cross(k, h)
    ny = n[1]
    nx = n[0]
    if ny >=0 :
        RAAN = np.arccos(nx / norm(n)) * ureg.radians
    else:
        RAAN = (2 * np.pi - np.arccos(nx / norm(n))) * ureg.radians
    semimaj = -mu / (2 * specific_total_energy)
    eccentricity = np.sqrt(1 + (2 * specific_total_energy * np.inner(h, h) * (ureg.meter ** 4 / ureg.second ** 2) / mu ** 2))
    inclination = -np.arccos(h[2] / (norm(h.magnitude) * h.units))
    eccentricity_vector = np.cross(velocity, h)* h.units * velocity.units / mu - location / norm(location) / location.units
    argument_of_periapsis = np.arccos(eccentricity_vector / norm(eccentricity_vector))
    true_anomalie = np.arccos(eccentricity_vector, location)
    kepler_paremters = {
        'e': eccentricity,
        'a': semimaj,
        'i': inclination,
        'RAAN': RAAN,
        'w': argument_of_periapsis,
        'nu': true_anomalie
    }
    return kepler_paremters
    # vernal equanox is in the x direction

def convert_kepler_to_cartesian(kepler):
    planar_r = (kepler['a']) * (1 - kepler['e'] ** 2) / (1 + kepler['e'] * np.cos(kepler['nu'])) * rotate_about_z(np.array([1, 0, 0]), kepler['nu'])
    planar_v = np.sqrt(mu / kepler['a'].to(ureg.meter) / (1 - kepler['e'] ** 2)) * (-np.sin(kepler['nu']), kepler['e'] + np.cos(kepler['nu']), 0)
    rotation = multi_dot([get_z_rotation_matrix(-kepler['RAAN']), get_x_rotation_matrix(-kepler['i']), get_z_rotation_matrix(-kepler['w'])])
    r = np.matmul(rotation, planar_r.to(ureg.meters)) * ureg.meter
    v = np.matmul(rotation, planar_v.to(ureg.meter / ureg.second)) * ureg.meter / ureg.seconds
    cartesian = {
        'x': r[0],
        'y': r[1],
        'z': r[2],
        'vx': v[0],
        'vy': v[1],
        'vz': v[2]
    }
    return cartesian

def get_eccentric_from_true(nu, e):
    eccentric_anomaly = 2 * np.arctan(np.tan(nu / 2) * np.sqrt((1 - e) / (1 + e)))
    return eccentric_anomaly

def get_true_from_eccentric(E ,e):
    nu = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * tan(E / 2))
    return E

def get_mean_anomalie_from_true_anomalie(nu, e):
    E = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(nu / 2))
    M = E - e * np.sin(E)
    return M

def approximate_eccentric_anomalie(M, e, approximation_level=5):
    start = time.time()
    E = [M]
    for i in range(0, approximation_level):
        E.append(M + e * np.sin(E[-1]))
    
    end = time.time()
    # print(end - start)
    return E[-1]

def approximate_eccentric_angular_velocity(period, E, e, approximation_level=5):
    E_dot = M_dot / (1 - e * cos(E))








