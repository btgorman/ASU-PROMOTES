# Copyright 2017 Brandon Tom Gorman

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import CoolProp.CoolProp as CP
import components
from openpyxl import load_workbook
import math

initial_conditions = 'input 100MW.xlsx'

def main():
    dictionary = {}

    balance_design_point(dictionary)

    return 0

#Design point calculations are currently not working properly, instead we are using time series at 900 for design point values. TODO: figure out how to make it work for less effort
def balance_design_point(dictionary):
    wb = load_workbook(initial_conditions, data_only=True)
    ws = wb['Sheet1']  # ws is now an IterableWorksheet
    for x in ws.rows:
        dictionary[x[0].value] = x[4].value

    SR3_temps = [[],[],[],[],[]] # 0->inner wall, 1->interm inner wall, 2->interm outer wall, 3->outer wall, 4->ambient
    ROx_temps = [[],[],[],[]] # 0->inner, 1->inner wall, 2->outer wall, 3-> ambient
    HS_temps = [[],[],[],[]] # 0->inner, 1->inner wall, 2->outer wall, 3-> ambient
    CS_temps = [[],[],[],[]] # 0->inner, 1->inner wall, 2->outer wall, 3-> ambient

    mol_abo3_max = dictionary['mol_abo3_max']
    solar_multiple = dictionary['solar_multiple']
    dictionary['A_sf'] = dictionary['A_sf'] * solar_multiple
    components.solve_delta(dictionary)

    dictionary['temp_6'] = 1044.0 + 273.15
    dictionary['temp_12'] = 386.9 + 273.15

    dictionary['mol_abo3_12'] = mol_abo3_max * 0.5
    dictionary['mol_air_12'] = ((mol_abo3_max - dictionary['mol_abo3_12']) * dictionary['rho_abo3'] / dictionary['dens_abo3']) * (1.0 + dictionary['ullage']) * dictionary['p0air'] / (dictionary['temp_12']*dictionary['R'])
    cp_air_12 = CP.PropsSI('CPMOLAR', 'T', dictionary['temp_12'],'P', dictionary['p0air'], 'Air')
    dictionary['energy_12'] = dictionary['mol_abo3_12'] * dictionary['cp_abo3'] * (dictionary['temp_12']-dictionary['temp_0']) + dictionary['mol_air_12'] * cp_air_12 * (dictionary['temp_12']-dictionary['temp_0'])
    dictionary['energy_12'] = dictionary['energy_12'] / 3600.0

    dictionary['mol_abo3_6'] = mol_abo3_max * 0.5
    dictionary['mol_N_6'] = ((mol_abo3_max - dictionary['mol_abo3_6']) * dictionary['rho_abo3'] / dictionary['dens_abo3']) * (1.0 + dictionary['ullage']) * dictionary['p0air'] / (dictionary['temp_6']*dictionary['R'])
    cp_N_6 = CP.PropsSI('CPMOLAR', 'T', dictionary['temp_6'],'P', dictionary['p0air'], 'Nitrogen')
    dictionary['energy_6'] = dictionary['mol_abo3_6'] * dictionary['cp_abo3'] * (dictionary['temp_6']-dictionary['temp_0']) + dictionary['mol_N_6'] * cp_N_6 * (dictionary['temp_6']-dictionary['temp_0']) + dictionary['mol_abo3_6'] * 0.5 * dictionary['delta'] * dictionary['H_rxn']
    dictionary['energy_6'] = dictionary['energy_6'] / 3600.0

    dictionary['i'] = 0
    dictionary['temp_7'] = dictionary['temp_6']
    dictionary['temp_13'] = dictionary['temp_12']
    
    dni = 900.0
    
    # Solve power block
    components.solve_power_block(dictionary)
    
    # Step 3. Solve ROx
    components.calibrate_ROx(dictionary)

    # Step 2. Balance SR3 and HX
    mols_delta = 1.0
    while mols_delta > 0.001:
        delta0 = dictionary['mol_abo3_1413']
        components.solve_SR3(dictionary, dni, SR3_temps)
        UA_HX, HX_eff = components.calibrate_HX(dictionary)
        delta1 = dictionary['mol_abo3_1413']
        mols_delta = math.fabs(delta0 - delta1)
    
    # Step 4. Solve HS + nitrogen blanket
    r_sb, A_hs, volume_hs = components.solve_HS(dictionary, mol_abo3_max, HS_temps)
   
    # Step 5. Solve CS
    components.solve_cold_storage(dictionary, CS_temps)
    
    # Step 6. Solve lift
    components.solve_work_lift(dictionary)
    
    # Step 7. Solve oxygen pump
    components.vacuum_pump(dictionary)

    print('SR3: Solar field area is', dictionary['A_sf'], 'm^2 at solar multiple of', solar_multiple)
    print('SR3: The number of receivers is', dictionary['number_of_receivers'], '(rounded up) at the input concentration factor of', dictionary['concentration_factor'], 'MW/m^2')
    print('HS/CS: If the max number of particle mols is', mol_abo3_max, ', the volume of the storage units are', volume_hs, 'm^3')
    print('Storage: Expected number of moles for 1 hour capacity', dictionary['mol_abo3_710']*3600)
    print('ROx: At particle T_high of', dictionary['temp_7']-273.15, 'C and particle T_low of', dictionary['temp_10']-273.15,'C one pipe has length and diameter of', dictionary['L_ROx'], dictionary['D_ROx'], 'm')
    print('HX: At DNI', dni, 'and particle inlet temp', dictionary['temp_13']-273.15, 'C and effectiveness of', HX_eff, 'the UA value is', UA_HX)
    print('')
    print('System: SR3-ROx productivity ratio', 100.0*dictionary['mol_abo3_1413']/dictionary['mol_abo3_710'], '%')
    print('')

    return 0

if __name__ == '__main__':
    main()