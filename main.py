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

import csv
import CoolProp.CoolProp as CP
import components
from openpyxl import load_workbook
import math

import numpy as np
import matplotlib.pyplot as plt

import scipy as sp
import scipy.stats

initial_conditions = 'input 100MW.xlsx'
solar_field_area = 'not_perfect'

def main():
    dictionary = {}
    time_balance_output_list = {}
    design_balance_output_list = {}

    time_balance_output_list, dni_array, volume_hs, energy_16max, r_sb, A_hs, mol_abo3_max, A_sf, SR3_temps, ROx_temps, HS_temps, CS_temps = balance_timeseries(dictionary)
    SA_rox = dictionary['L_ROx'] * 2.0 * dictionary['D_ROx']

    cost_list = cost_calculations(SA_rox, dictionary, volume_hs, mol_abo3_max, A_sf, energy_16max, r_sb, A_hs)

    tot_approx_op_hours = predict_dispatch_schedule(time_balance_output_list, dni_array)

    # Error checking time balance output list
    for listelem in time_balance_output_list:
        for elem in listelem:
            if elem < -0.01:
                print('time balance elem = {!s}'.format(elem))
                print('')

    #Step 5.0 Write Excel output
    fl = open('OUTPUTPROMOTES.csv', 'w')
    writer = csv.writer(fl)
    for values in time_balance_output_list:
        writer.writerow(values)
        # writer.writerow('')
    for values in enumerate(cost_list):
        writer.writerow(values)
        # writer.writerow('')
    fl.close()
    
    print('')
    print('Current solar field area of', dictionary['A_sf'], 'm^2 at solar multiple of', dictionary['solar_multiple'])
    print('')
    SR3_eta_num = sum(time_balance_output_list[3][:]) - sum(time_balance_output_list[12][:])
    SR3_eta_den = sum(time_balance_output_list[1][:])
    ROx_eta_num = sum(time_balance_output_list[10][:]) - sum(time_balance_output_list[8][:])
    ROx_eta_den = sum(time_balance_output_list[6][:]) - sum(time_balance_output_list[9][:])
    ROx_eta = ROx_eta_num / ROx_eta_den
    print('SR3 efficiency is', 100.0 * SR3_eta_num / SR3_eta_den, '%')
    print('ROx efficiency is', 100.0 * ROx_eta, '%')
    
    sys_th_eta_num = sum(time_balance_output_list[10][:]) - sum(time_balance_output_list[8][:]) + ROx_eta * dictionary['e_hs_end'] + dictionary['e_cs_end']
    sys_th_eta_den = sum(time_balance_output_list[1][:])  + ROx_eta * dictionary['e_hs_begin'] + dictionary['e_cs_begin']
    sys_th_eta = sys_th_eta_num / sys_th_eta_den
    sys_ap_eta_num = sum(time_balance_output_list[14][:]) + sum(time_balance_output_list[15][:])
    sys_ap_eta_den = tot_approx_op_hours * 111.7 * 1000.0 * 1000.0 # 111.7 MWe rated power
    sys_ap_eta = 1.0 - (sys_ap_eta_num / sys_ap_eta_den)
    print('System-boundary efficiency is', 100.0 * sys_th_eta * sys_ap_eta, '%')

    # plot_temperatures(SR3_temps, ROx_temps, HS_temps, CS_temps)

    return 0

#Step 3: Calculate the Time Series balance - this is the meet and bones of the code
def balance_timeseries(dictionary):
    wb = load_workbook(initial_conditions, data_only=True)
    ws = wb['Sheet1']  # ws is now an IterableWorksheet
    for x in ws.rows:
        dictionary[x[0].value] = x[4].value

    dispatch = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for elem in range(0, len(dispatch)):
        dispatch[elem] = 1
    # dispatch[12] = 0
    print('Running with dispatch schedule:', dispatch)
    print('')
  
    time_balance_outputlist = [[0 for i in range(24)] for y in range(42)]
    SR3_temps = [[],[],[],[],[]] # 0->inner wall, 1->interm inner wall, 2->interm outer wall, 3->outer wall, 4->ambient
    ROx_temps = [[],[],[],[]] # 0->inner, 1->inner wall, 2->outer wall, 3-> ambient
    HS_temps = [[],[],[],[]] # 0->inner, 1->inner wall, 2->outer wall, 3-> ambient
    CS_temps = [[],[],[],[]] # 0->inner, 1->inner wall, 2->outer wall, 3-> ambient

    mol_abo3_max = dictionary['mol_abo3_max']
    solar_multiple = dictionary['solar_multiple']
    dictionary['A_sf'] = dictionary['A_sf'] * solar_multiple
    A_sf = dictionary['A_sf']
    components.solve_delta(dictionary)

    dictionary['mol_abo3_12'] = mol_abo3_max
    dictionary['mol_air_12'] = ((mol_abo3_max - dictionary['mol_abo3_12']) * dictionary['rho_abo3'] / dictionary['dens_abo3']) * (1.0 + dictionary['ullage']) * dictionary['p0air'] / (dictionary['temp_12']*dictionary['R'])
    cp_air_12 = CP.PropsSI('CPMOLAR', 'T', dictionary['temp_12'],'P', dictionary['p0air'], 'Air')
    dictionary['energy_12'] = dictionary['mol_abo3_12'] * dictionary['cp_abo3'] * (dictionary['temp_12']-dictionary['temp_0']) + dictionary['mol_air_12'] * cp_air_12 * (dictionary['temp_12']-dictionary['temp_0'])
    dictionary['energy_12'] = dictionary['energy_12'] / 3600.0

    dictionary['mol_abo3_6'] = 0.0
    dictionary['mol_N_6'] = ((mol_abo3_max - dictionary['mol_abo3_6']) * dictionary['rho_abo3'] / dictionary['dens_abo3']) * (1.0 + dictionary['ullage']) * dictionary['p0air'] / (dictionary['temp_6']*dictionary['R'])
    cp_N_6 = CP.PropsSI('CPMOLAR', 'T', dictionary['temp_6'],'P', dictionary['p0air'], 'Nitrogen')
    dictionary['energy_6'] = dictionary['mol_abo3_6'] * dictionary['cp_abo3'] * (dictionary['temp_6']-dictionary['temp_0']) + dictionary['mol_N_6'] * cp_N_6 * (dictionary['temp_6']-dictionary['temp_0']) + dictionary['mol_abo3_6'] * 0.5 * dictionary['delta'] * dictionary['H_rxn']
    dictionary['energy_6'] = dictionary['energy_6'] / 3600.0

    dictionary['e_hs_begin'] = dictionary['energy_6']
    dictionary['e_cs_begin'] = dictionary['energy_12']
    
    # STEP 1. CALL FOR 24 HOUR LOOP to solve for the "real time" analysis for a full day
    for i in range(24):
        dictionary['i'] = i
        dictionary['temp_7'] = dictionary['temp_6']
        dictionary['temp_13'] = dictionary['temp_12']
        dictionary['mol_abo3_1413'] = 0.0
        dictionary['mol_abo3_710'] = 0.0

        dni, dni_array = components.set_dni(i)

        # Solve power block
        components.solve_power_block(dictionary)
        
        # Step 3. Solve ROx
        components.solve_ROx(dictionary, ROx_temps, dispatch)

        # Step 2. Balance SR3 and HX
        mols_delta = 1.0
        while mols_delta > 0.001:
            delta0 = dictionary['mol_abo3_1413']
            components.solve_SR3(dictionary, dni, SR3_temps)
            components.solve_HX(dictionary)
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
        energy_16max = dictionary['energy_16max']

        # Step 8. Record everything
        time_balance_outputlist[0][i]=dictionary['energy_1']
        time_balance_outputlist[1][i]=dictionary['energy_2']
        time_balance_outputlist[2][i]=dictionary['energy_3']
        time_balance_outputlist[3][i]=dictionary['energy_4']
        time_balance_outputlist[4][i]=dictionary['energy_5']
        time_balance_outputlist[5][i]=dictionary['energy_6']
        time_balance_outputlist[6][i]=dictionary['energy_7']
        time_balance_outputlist[7][i]=dictionary['energy_8']
        time_balance_outputlist[8][i]=dictionary['energy_9']
        time_balance_outputlist[9][i]=dictionary['energy_10']
        time_balance_outputlist[10][i]=dictionary['energy_11']
        time_balance_outputlist[11][i]=dictionary['energy_12']
        time_balance_outputlist[12][i]=dictionary['energy_13']
        time_balance_outputlist[13][i]=dictionary['energy_14']
        time_balance_outputlist[14][i]=dictionary['energy_15']
        time_balance_outputlist[15][i]=dictionary['energy_16']
        time_balance_outputlist[16][i]=dictionary['energy_17']
        time_balance_outputlist[17][i]=dictionary['energy_18']
        time_balance_outputlist[18][i]=dictionary['energy_19']
        time_balance_outputlist[19][i]=dictionary['energy_20']
        time_balance_outputlist[20][i]=dictionary['energy_21']
        time_balance_outputlist[21][i]=dictionary['energy_22']
        time_balance_outputlist[22][i]=dictionary['energy_23']
        time_balance_outputlist[23][i]=dictionary['energy_24']
        time_balance_outputlist[24][i]=dictionary['energy_25']
        time_balance_outputlist[25][i]=dictionary['temp_1']
        time_balance_outputlist[26][i]=dictionary['temp_3']
        time_balance_outputlist[27][i]=dictionary['temp_4']
        time_balance_outputlist[28][i]=dictionary['temp_6']
        time_balance_outputlist[29][i]=dictionary['temp_7']
        time_balance_outputlist[30][i]=dictionary['temp_10']
        time_balance_outputlist[31][i]=dictionary['temp_12']
        time_balance_outputlist[32][i]=dictionary['temp_13']
        time_balance_outputlist[33][i]=dictionary['mol_abo3_1413']
        time_balance_outputlist[34][i]=dictionary['mol_abo3_710']
        time_balance_outputlist[35][i]=dictionary['mol_abo3_6']
        time_balance_outputlist[36][i]=dictionary['mol_abo3_12']
        time_balance_outputlist[37][i]=dictionary['mol_N_22']
        time_balance_outputlist[38][i]=dictionary['mol_N_23']
        time_balance_outputlist[39][i]=dictionary['mol_air_24']
        time_balance_outputlist[40][i]=dictionary['mol_air_25']
        time_balance_outputlist[41][i]=dni

    dictionary['e_hs_end'] = time_balance_outputlist[5][-1]
    dictionary['e_cs_end'] = time_balance_outputlist[11][-1]
    mol_abo3_6max_a=dictionary['mol_abo3_6max']
    return time_balance_outputlist, dni_array, volume_hs, energy_16max, r_sb, A_hs, mol_abo3_max, A_sf, SR3_temps, ROx_temps, HS_temps, CS_temps

#Step 4. CALCULATE the cost of everything.
def cost_calculations(SA_rox,dictionary,volume_hs, mol_abo3_6max,A_sf,energy_16max,r_sb, A_hs):
    cost_SR3,cost_vp,cost_particles,cost_tower,cost_elevator,cost_HX,cost_sH,cost_sLH,cost_sUH,cost_Rox,cost_control,cost_sf,cost_powerblock,cost_bal,cost_contingency,cost_own,cost_mctot,cost_total_direct,cost_tower_receiver,cost_thermalstorage,cost_solar_control,cost_gen_rox,cost_balance_of_Plant_total,cost_ppe,capacity_factor,estimated_production,cost_om,cost_particle_replacement,capital_cost_year,energy_density, energy_stored, cost_energy, cost_perkWth_TOTAL  =components.solve_cost(SA_rox, r_sb, A_hs,dictionary['P_st'],mol_abo3_6max,dictionary['R_rec'],dictionary['D_ap'],dictionary['pi'],dictionary['F_con'], dictionary['F_rec_par'],dictionary['mol_mass_abo3'],dictionary['hd'],volume_hs, dictionary['w0'],dictionary['w1'],dictionary['w2'],dictionary['w3'],dictionary['w4'],dictionary['cost_lay'],dictionary['cost_insfi'],dictionary['cost_pertcon'],dictionary['cost_expbo'],dictionary['cost_reincon'],dictionary['F_ms'], dictionary['C_n2g'] ,dictionary['cost_SR3_a'], dictionary['M_SR3'],dictionary['vp_0'], dictionary['vp_1'],dictionary['Per_s'],dictionary['M_ec'],dictionary['M_p'],dictionary['M_SP'], dictionary['C_ppart'], dictionary['F_tpre'],dictionary['F_sca'],dictionary['C_ele_i'],dictionary['hx_0'],dictionary['hx_1'],dictionary['area_hx'], dictionary['F_sv'], dictionary['F_SR3C'],dictionary['CF_uh'],dictionary['F_v'],dictionary['r_rox'],dictionary['L_rox'],dictionary['M_roxcomp'],dictionary['P_r'],dictionary['cost_omi'],dictionary['F_cown'],dictionary['contingency'],dictionary['B_psf'],dictionary['B_ps'],dictionary['F_ti'],dictionary['F_tc'],dictionary['F_tp'],dictionary['F_ts'],dictionary['F_comp'],dictionary['cost_field'],dictionary['cost_fieldprep'],A_sf,dictionary['area_receiver'],energy_16max,dictionary['cp_abo3'],dictionary['temp_7'])
    #set a list to store costs
    cost_list = [0 for i in range(35)]
    #listing all of the costs
    cost_list[0]= cost_SR3
    cost_list[1]=cost_vp
    cost_list[2]=cost_particles
    cost_list[3]=cost_tower
    cost_list[4]=cost_elevator
    cost_list[5]=cost_HX
    cost_list[6]=cost_sH
    cost_list[7]=cost_sLH
    cost_list[8]=cost_sUH
    cost_list[9]=cost_Rox
    cost_list[10]=cost_control
    cost_list[11]=cost_sf
    cost_list[12]=cost_powerblock
    cost_list[13]=cost_bal
    cost_list[14]=cost_contingency
    cost_list[15]=cost_own
    cost_list[16]=cost_mctot
    cost_list[17]=0.0
    cost_list[18]=cost_total_direct
    cost_list[19]=cost_tower_receiver
    cost_list[20]=cost_thermalstorage
    cost_list[21]=cost_solar_control
    cost_list[22]=cost_gen_rox
    cost_list[23]=cost_balance_of_Plant_total
    cost_list[24]=cost_ppe
    cost_list[25]=capacity_factor
    cost_list[26]=estimated_production
    cost_list[27]= cost_om
    cost_list[28]=cost_particle_replacement
    cost_list[29]= capital_cost_year
    cost_list[30]= volume_hs
    cost_list[31] = energy_density
    cost_list[32] = energy_stored
    cost_list[33]= cost_energy
    cost_list [34] = cost_perkWth_TOTAL

    return cost_list

def predict_dispatch_schedule(system_states, dni_array):

    off_sun_dispatch = [0] * len(dni_array)
    on_sun_dispatch = [0] * len(dni_array)

    storage_target = 6
    approx_op_hours = sum(system_states[33][:]) / max(system_states[34][:])
    tot_approx_op_hours = approx_op_hours
    print('Approx operational hours:', approx_op_hours, '\n')

    index_off_sun = 12
    for iter in range(1, len(dni_array)+1):
        if dni_array[-iter] > 350.0:
            index_off_sun = len(dni_array) - iter
            break

    stop_off_gen = min(len(dni_array) - index_off_sun, storage_target)
    next_day_gen = storage_target - stop_off_gen
    approx_op_hours += next_day_gen

    for iter in range(index_off_sun, index_off_sun + stop_off_gen):
        if approx_op_hours >= 1:
            off_sun_dispatch[iter+1] = 1.0 * 100
            approx_op_hours -= 1
        else:
            break

    print('OFF SUN DISPATCH', off_sun_dispatch)

    for iter in range(0, index_off_sun+1):
        dec = index_off_sun - iter
        if approx_op_hours >= 1:
            on_sun_dispatch[dec] = 1.0 * 100
            approx_op_hours -= 1
        else:
            on_sun_dispatch[dec] = approx_op_hours * 100
            break

    print('ON SUN DISPATCH', on_sun_dispatch)

    return tot_approx_op_hours

def plot_temperatures(SR3_temps, ROx_temps, HS_temps, CS_temps):

    def mean_confidence_interval(data, confidence=0.95):
        a = 1.0 * np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * sp.stats.t.ppf((1.0+confidence)/2., n-1)
        return m, h

    def median(data):
        a = 1.0 * np.array(data)
        m = np.median(a)
        return m

    SR3_X = np.linspace(0, 1, 5, endpoint=True)
    ROx_X = np.linspace(0, 1, 4, endpoint=True)
    HS_X = np.linspace(0, 1, 4, endpoint=True)
    CS_X = np.linspace(0, 1, 4, endpoint=True)

    SR3_X_tick = [0.0, 1.0/4.0, 2.0/4.0, 3.0/4.0, 1.0]
    SR3_X_tick_label = ['Inner Surface', 'Evacuated Inner Surface', 'Evacuated Outer Surface', 'Exterior Surface', 'Ambient']
    ROx_X_tick = [0.0, 1.0/3.0, 2.0/3.0, 1.0]
    ROx_X_tick_label = ['Interior', 'Inner Surface', 'Exterior Surface', 'Ambient']
    HS_X_tick = [0.0, 1.0/3.0, 2.0/3.0, 1.0]
    HS_X_tick_label = ['Interior', 'Inner Surface', 'Exterior Surface', 'Ambient']
    CS_X_tick = [0.0, 1.0/3.0, 2.0/3.0, 1.0]
    CS_X_tick_label = ['Interior', 'Inner Surface', 'Exterior Surface', 'Ambient']

    SR3_Y = []
    SR3_Y_med = []
    SR3_Y_err = []
    ROx_Y = []
    ROx_Y_med = []
    ROx_Y_err = []
    HS_Y = []
    HS_Y_med = []
    HS_Y_err = []
    CS_Y = []
    CS_Y_med = []
    CS_Y_err = []

    for label in SR3_temps:
        Y, error = mean_confidence_interval(label)
        Y_med = median(label)
        SR3_Y.append(Y)
        SR3_Y_med.append(Y_med)
        SR3_Y_err.append(error)

    for label in ROx_temps:
        Y, error = mean_confidence_interval(label)
        Y_med = median(label)
        ROx_Y.append(Y)
        ROx_Y_med.append(Y_med)
        ROx_Y_err.append(error)

    for label in HS_temps:
        Y, error = mean_confidence_interval(label)
        Y_med = median(label)
        HS_Y.append(Y)
        HS_Y_med.append(Y_med)
        HS_Y_err.append(error)

    for label in CS_temps:
        Y, error = mean_confidence_interval(label)
        Y_med = median(label)
        CS_Y.append(Y)
        CS_Y_med.append(Y_med)
        CS_Y_err.append(error)

    plt.xticks(SR3_X_tick, SR3_X_tick_label)
    plt.xlim(-0.1, 1.1)
    plt.ylim(0.0, 1200.0)
    plt.ylabel('Celsius')
    plt.title('SR3 Temperature Map')
    plt.errorbar(SR3_X, SR3_Y, yerr=SR3_Y_err, linewidth=1.0, label='mean', color='#00acd3')
    plt.plot(SR3_X, SR3_Y_med, linewidth=1.0, label='median')
    plt.style.use('seaborn-pastel')
    plt.legend(loc='upper right', frameon=False, fontsize='small')
    plt.show()

    plt.xticks(ROx_X_tick, ROx_X_tick_label)
    plt.xlim(-0.1, 1.1)
    plt.ylim(0.0, 1200.0)
    plt.ylabel('Celsius')
    plt.title('ROx Temperature Map')
    plt.errorbar(ROx_X, ROx_Y, yerr=ROx_Y_err, linewidth=1.0, label='mean',color='#00acd3')
    plt.plot(ROx_X, ROx_Y_med, linewidth=1.0, label='median')
    plt.style.use('seaborn-pastel')
    plt.legend(loc='upper right', frameon=False, fontsize='small')
    plt.show()

    plt.xticks(HS_X_tick, HS_X_tick_label)
    plt.xlim(-0.1, 1.1)
    plt.ylim(0.0, 1200.0)
    plt.ylabel('Celsius')
    plt.title('HS Temperature Map')
    plt.errorbar(HS_X, HS_Y, yerr=HS_Y_err, linewidth=1.0, label='mean',color='#00acd3')
    plt.plot(HS_X, HS_Y_med, linewidth=1.0, label='median')
    plt.style.use('seaborn-pastel')
    plt.legend(loc='upper right', frameon=False, fontsize='small')
    plt.savefig('image.png')
    plt.show()

    plt.xticks(CS_X_tick, CS_X_tick_label)
    plt.xlim(-0.1, 1.1)
    plt.ylim(0.0, 1200.0)
    plt.ylabel('Celsius')
    plt.title('CS Temperature Map')
    plt.errorbar(CS_X, CS_Y, yerr=CS_Y_err, linewidth=1.0, label='mean',color='#00acd3')
    plt.plot(CS_X, CS_Y_med, linewidth=1.0, label='median')
    plt.style.use('seaborn-pastel')
    plt.legend(loc='upper right', frameon=False, fontsize='small')
    plt.show()

def solveTH(dictionary):
    a = math.log((dictionary['p0o2']+0.0)/dictionary['p502']);
    b = (dictionary['H_rxn']+0.0) / (dictionary['R'] * dictionary['temp_4']);
    z = (7.0 - 2.0*a + 2.0*b) / (7.0 + 2.0*b);

    TH = dictionary['temp_4'] / z;
    return TH

if __name__ == '__main__':
    main()