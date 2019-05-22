import numpy as np
import time
import json
from numpy.linalg import norm as norm
from numpy.linalg import multi_dot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from globalRegistry import GlobalRegistry
import itertools
from math import isclose

ureg = GlobalRegistry.getRegistry().ureg
constants = json.load(open("constants.json"))
earths_radius = constants["earth radius"] * ureg.meters
mu = constants["mu"] * ureg.meter ** 3 / ureg.seconds ** 2
earths_mass = constants["earth mass"] * ureg.kilograms
j2 = constants["j2"]
earth_period = constants["earth period"]
        
    
def get_x_rotation_matrix(amount):
    c = np.cos(amount)
    s = np.sin(amount)
    z = np.zeros(np.max(amount.shape)) 
    o = np.ones(np.max(amount.shape))

    x_rotation_matrix = np.array(
        [
            np.vstack((o, z, z)),
            np.vstack((z ,c, -s)),
            np.vstack((z, s, c))
        ]
    )
    return x_rotation_matrix

def get_y_rotation_matrix(amount):
    c = np.cos(amount)
    s = np.sin(amount)
    z = np.zeros(np.max(amount.shape))
    o = np.ones(np.max(amount.shape))

    y_rotation_matrix = np.array(
        [
            np.vstack((c, z, s)),
            np.vstack((z ,o, z)),
            np.vstack((-s, z, c))
        ]
    )
    return y_rotation_matrix

def get_z_rotation_matrix(amount):
    c = np.cos(amount)
    s = np.sin(amount)
    z = np.zeros(np.max(amount.shape))
    o = np.ones(np.max(amount.shape))
    z_rotation_matrix = np.array(    
        [
            np.vstack((c, -s, z)),
            np.vstack((s ,c, z)),
            np.vstack((z, z, o))
        ]
    )

    return z_rotation_matrix

def rotate_about_x(vector, amount):
    x_rotation_matrix = get_x_rotation_matrix(amount)
    rotated_vector = np.matmul(vector, x_rotation_matrix)
    return rotated_vector

def rotate_about_y(vector, amount):
    y_rotation_matrix = get_y_rotation_matrix(amount)
    rotated_vector = np.matmul(vector, y_rotation_matrix)
    return rotated_vector

def rotate_about_z(vector, amount):
    z_rotation_matrix = get_z_rotation_matrix(amount.to(ureg.radians).magnitude)
    try:
        rotated_vector = np.einsum('ijk, i -> jk', z_rotation_matrix, vector)
    except:
        print('d')

    return rotated_vector


def convert_cartesian_to_kepler(cartesian):
    velocity = cartesian['v']
    location = cartesian['r']
    from orbit import Orbit
    specific_total_energy = Orbit._get_orbital_energy(cartesian)
    h = Orbit._get_orbital_angular_momentum(cartesian) #reduced_angular_momentum
    k = [0, 0, 1]
    n = np.cross(k, h.T).T
    if isclose(norm(n, axis=0)[0], 0, abs_tol=1e-7):
        RAAN = 0 * ureg.radians        
    else:
        ny = n[1, :]
        nx = n[0, :]
        RAAN = np.arccos(nx / norm(n, axis=0))
        RAAN = np.where(ny > 0, RAAN, 2 * np.pi - RAAN) * ureg.radians

    semimaj = -mu / (2 * specific_total_energy)
    eccentricity = np.sqrt(1 + np.round(20000 * specific_total_energy * np.einsum('ij, ij -> j', h, h) * mu ** -2) / 10000 * h.units ** 2 )
    inclination = -np.arccos(h[2, :].magnitude / norm(h.magnitude, axis=0)) * ureg.radians
    eccentricity_vector = np.cross(velocity.T, h.T) * h.units * velocity.units / mu - location.T / norm(location, axis=1) / location.units
    if not isclose(norm(n, axis=0)[0], 0, abs_tol=1e-7):
        argument_of_periapsis = np.arccos(np.einsum('ij, ij -> j', eccentricity_vector.T, n) / norm(eccentricity_vector, axis=1)/ norm(n, axis=0))
    else:
        argument_of_periapsis = np.arctan(eccentricity_vector[1, :] / eccentricity_vector[0, :])

    true_anomalie = np.arccos(np.einsum('ij, ij -> j',eccentricity_vector.T, location) / norm(eccentricity_vector, axis=1)/ norm(location, axis=0)) * ureg.radians
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
    #converted to numpy
    time_steps = np.max(kepler['e'].shape)
    planar_r = (kepler['a']) * (1 - kepler['e'] ** 2) / (1 + kepler['e'] * np.cos(kepler['nu'])) * rotate_about_z(np.array([1, 0, 0]), -kepler['nu'])
    velocity_Vec = np.array([np.array(-np.sin(kepler['nu']).magnitude), np.array(kepler['e'] + np.cos(kepler['nu']).magnitude), np.zeros([time_steps])])
    planar_v = np.sqrt(mu / kepler['a'].to(ureg.meter) / (1 - kepler['e'] ** 2)) * velocity_Vec
    X = get_x_rotation_matrix(-kepler['i'])
    Z_RAAN = get_z_rotation_matrix(-kepler['RAAN'])
    Z_w = get_z_rotation_matrix(-kepler['w'])
    rotation =  np.einsum('ijn, jkn, kln -> iln', Z_RAAN, X, Z_w)
    r = np.einsum('jin, in -> jn', rotation, planar_r.to(ureg.meters).magnitude) * ureg.meter
    v = np.einsum('jin, in ->jn', rotation, planar_v.to(ureg.meter / ureg.second)) * ureg.meter / ureg.seconds
    cartesian = {
        'r': r,
        'v': v
        # 'x': r.magnitude[0, :] * r.units,
        # 'y': r.magnitude[1, :] * r.units,
        # 'z': r.magnitude[2, :] * r.units,
        # 'vx': v.magnitude[0, :] * v.units,
        # 'vy': v.magnitude[1, :] * v.units,
        # 'vz': v.magnitude[2, :] * v.units
    }
    return cartesian

def get_eccentric_from_true(nu, e):
    eccentric_anomaly = 2 * np.arctan(np.tan(nu / 2) * np.sqrt((1 - e) / (1 + e)))
    return eccentric_anomaly

def get_true_from_eccentric(E ,e):
    nu = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(E / 2))
    return E

def get_mean_anomalie_from_true_anomalie(nu, e):
    E = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(nu / 2))
    M = E - e * np.sin(E)
    return M

def approximate_eccentric_anomalie(M, e, approximation_level=5):
    start = time.time()
    E = [M]
    for _ in itertools.repeat(None, approximation_level):
        E.append(M + e * np.sin(E[-1]))
    
    end = time.time()
    return E[-1]

def approximate_eccentric_angular_velocity(period, E, e, approximation_level=5):
    E_dot = period / (1 - e * np.cos(E))








