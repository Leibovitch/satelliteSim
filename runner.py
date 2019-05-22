from orbit import Orbit
import json
from satellite import Satellite
from datetime import timedelta, datetime
from globalRegistry import GlobalRegistry
import matplotlib.pyplot as plt
import numpy as np

ureg = GlobalRegistry.getRegistry().ureg
constants = json.load(open("constants.json"))
Re = constants["earth radius"] * ureg.meters

def convert_satellite_from_db_to_numpy_pint(database_entry):
    # test that apogee is larger then perigee
    satellite_summary = {}
    for key in database_entry:
        if (key == 'kepler_parameters'):
            database_kepler_params = database_entry[key]
            kepler_parameters = {}
            a_e_pair_existance = database_kepler_params['a']['value'] and database_kepler_params['e']
            per_apo_pair_exists = database_kepler_params['perigee']['value'] and database_kepler_params['apogee']['value']
            if (a_e_pair_existance and not per_apo_pair_exists):
                kepler_parameters['a'] = np.array([database_kepler_params['a']['value']]) * eval("ureg." + database_kepler_params['a']['units'])
                kepler_parameters[inner_key] = np.array([database_kepler_params[inner_key]])
            elif (per_apo_pair_exists and not a_e_pair_existance):
                per = database_kepler_params['perigee']['value'] * eval("ureg." + database_kepler_params['perigee']['units'])
                apo = database_kepler_params['apogee']['value'] * eval("ureg." + database_kepler_params['apogee']['units'])
                kepler_parameters['a'] = np.array(((per + apo + 2 * Re) / 2).magnitude) * eval("ureg." + database_kepler_params['a']['units'])
                kepler_parameters['e'] = np.array([((apo - per) / 2 / kepler_parameters['a']).magnitude])
            else:
                raise ValueError('y ou need to specify only a-e pair or perigee-apogee pair, not both ot nither')

            for inner_key in database_kepler_params:
                if inner_key not in ['e', 'a', 'perigee', 'apogee']:
                    kepler_parameters[inner_key] = np.array([database_kepler_params[inner_key]['value']]) * eval("ureg." + database_kepler_params[inner_key]['units'])

            satellite_summary[key] = kepler_parameters
        else:
            satellite_summary[key] = database_entry[key]['value'] * eval("ureg." + database_entry[key]['units'])

    return satellite_summary


def create_full_access_table(constellation, target_list, sim_length, sim_step):
    sat_num = len(constellation)
    target_num = len(target_list)
    access_table = np.empty((sat_num, target_num))
    for i, sat in zip(range(0, sat_num), constellation):
        for j, target in zip(range(0, target_num), target_list):
            access_table[i, j] = sat.get_access_to_target(target, sim_length, sim_step)

    return access_table

def main_sim(target_list, sim_length, sim_step):
    initial_condition = json.load(open("initial_conditions.json"))
    satellite_database = json.load(open("satellite_database.json"))
    constellation = []
    for satellite_type in initial_condition['consellation']:
        current_params = convert_satellite_from_db_to_numpy_pint(satellite_database[satellite_type])
        number_of_type = initial_condition[satellite_type]
        for i in range(0, number_of_type):
            current_kepler_params = current_params['kepler_parameters']
            nu_jumps = 2 * np.pi / number_of_type
            initial_nu = i * nu_jumps
            current_kepler_params['nu'] = np.array([initial_nu]) * ureg.radians
            current_sat = Satellite(current_params['res_at_nadir_at_500'],
            current_params['look_angle'],
            current_params['band'], current_kepler_params,
            datetime.now(), satellite_type + '_' + i)
            constellation.append(current_sat)

    access_table = create_full_access_table(constellation, target_list, sim_length, sim_step)
    return access_table
