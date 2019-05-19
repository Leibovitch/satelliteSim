from orbit import Orbit
import json
from satellite import Satellite
from datetime import timedelta, datetime
from globalRegistry import GlobalRegistry

ureg = GlobalRegistry.getRegistry().ureg

def convert_satellite_from_db_to_pint(database_entry):
    satellite_summary = {}
    for key in database_entry:
        print(key)
        if (key == 'kepler_parameters'):
            database_kepler_params = database_entry[key]
            kepler_parameters = {}
            for inner_key in database_kepler_params:
                print(inner_key)
                if inner_key == 'e':
                    kepler_parameters[inner_key] = database_kepler_params[inner_key]
                else:
                    kepler_parameters[inner_key] = database_kepler_params[inner_key]['value'] * eval("ureg." + database_kepler_params[inner_key]['units'])

            satellite_summary[key] = kepler_parameters
        else:
            satellite_summary[key] = database_entry[key]['value'] * eval("ureg." + database_entry[key]['units'])

    return satellite_summary

initial_condition = json.load(open("initial_conditions.json"))
satellite_database = json.load(open("satellite_database.json"))

constellation = []
for key in initial_condition:
    current_params = convert_satellite_from_db_to_pint(satellite_database[key])
    current_sat = Satellite(current_params['res_at_nadir_at_500'], current_params['look_angle'], current_params['band'], current_params['kepler_parameters'], datetime.now())
    start_time = datetime.now()
    current_sat.get_access_to_location(0 * ureg.degree, 0 * ureg.degree, 1 * ureg.meter, timedelta(weeks = 1), timedelta(minutes = 1))
    constellation.append(current_sat)
    end_time =datetime.now()

print(end_time - start_time)