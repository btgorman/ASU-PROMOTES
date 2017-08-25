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

import math
import CoolProp.CoolProp as CP

# def set_dni(i):
#     #dni_array =[0.0, 0.0, 0.0, 0.0, 00.0, 007.0, 225.0, 573.0, 440.0, 608.0, 746.0, 638.0, 772.0, 947.0, 818.0, 867.0, 637.0, 142.0, 000.0, 0.0, 0.0, 0.0, 0.0, 0.0] # The higher value is 947 unless we are trying to get the Design Point DNI, then put it to 900
#     #dni_array = [0,0,0,0,0,0,24,225,463,663,802,833,782,666,632,329,229,82,2,0,0,0,0,0] # south africa - need to figure out if this is right 833
#     #dni_array =[0.0, 0.0, 0.0, 0.0, 00.0, 00.0, 00.0, 900.0, 900, 900, 900, 900, 900, 900, 900, 900, 00.0, 00.0, 00.0, 0.0, 0.0, 0.0, 0.0, 0.0]

#     dni_array = [0.0, 0.0, 0.0, 0.0, 0.0, 201.0, 536.0, 284.0, 549.0, 675.0, 584.0, 714.0, 881.0, 942.0, 622.0, 696.0, 534.0, 361.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Barstow Spring and Autumn equinox -- April 17
#     # dni_array = [0.0, 0.0, 0.0, 0.0, 4.0, 302.0, 502.0, 708.0, 732.0, 807.0, 826.0, 835.0, 888.0, 695.0, 786.0, 773.0, 677.0, 499.0, 172.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Barstow Summer solstice -- July 19
#     # dni_array = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 420.0, 793.0, 834.0, 658.0, 627.0, 664.0, 457.0, 74.0, 322.0, 199.0, 86.0, 340.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Barstow Winter solstice -- March 12
#     dni = dni_array[i]

#     return dni, dni_array

def solve_power_block(dictionary):
    mass_flow_air = dictionary['mass_flow_air']
    molar_mass_air = dictionary['molar_mass_air']
    temp_0 = dictionary['temp_0']
    temp_9 = dictionary['temp_9']
    temp_11 = dictionary['temp_11']
    p0air = dictionary['p0air']
    pr_turbine = dictionary['pr_turbine']

    mol_air = mass_flow_air / molar_mass_air

    cp_air = CP.PropsSI('CPMOLAR', 'T', 0.5*(temp_9+temp_11), 'P', p0air*pr_turbine, 'Air')

    energy_9 = mol_air * cp_air * (temp_9-temp_0)
    energy_11 = mol_air * cp_air * (temp_11-temp_0)
    
    dictionary['energy_9'] = energy_9
    dictionary['energy_11'] = energy_11
    dictionary['mol_air'] = mol_air

    return 0

def solve_delta(dictionary):
    temp_4 = dictionary['temp_4']
    p0o2 = dictionary['p0o2']
    pSR3 = dictionary['pSR3']
    f0 = dictionary['f0']
    f1 = dictionary['f1']
    f2 = dictionary['f2']
    f3 = dictionary['f3']
    f4 = dictionary['f4']
    f6 = dictionary['f6']
    f7 = dictionary['f7']
    f8 = dictionary['f8']
    t_reference = 400.0 + 273.15

    U = 0.5 * math.log(p0o2/pSR3)
    B = (t_reference) / temp_4
    X = -math.log(B)
    Z = (f0+f1*X+f2*(math.exp(X)-1.0-X)+f3*U+f7*U*X)/(1.0+f4*X+f6*U+f8*U*X)

    delta = math.exp(-Z) - 0.000574

    dictionary['delta'] = delta

    return 0

def solve_SR3_qloss(dictionary, aperture_area, SR3_temps):
    dia_aper = dictionary['dia_aper']
    SR3_area_ratio = dictionary['SR3_area_ratio']
    SR3_insul_thick = dictionary['SR3_insul_thick']
    SR3_insul_mult = dictionary['SR3_insul_mult'] # UNUSED
    SR3_insul_cond = dictionary['SR3_insul_cond']
    SR3_insul_emiss = dictionary['SR3_insul_emiss']
    SR3_mainb_thick = dictionary['SR3_mainb_thick']
    SR3_mainb_mult = dictionary['SR3_mainb_mult'] # UNUSED
    SR3_mainb_cond = dictionary['SR3_mainb_cond']
    SR3_mainb_emiss = dictionary['SR3_mainb_emiss']
    temp_4 = dictionary['temp_4']
    temp_0 = dictionary['temp_0']
    p0air = dictionary['p0air']
    sigma = dictionary['sigma']
    i = dictionary['i']

    num_aper = math.ceil(aperture_area / (math.pi * (dia_aper*0.5)**2))
    dictionary['number_of_receivers'] = num_aper
    aperture_area = aperture_area / num_aper
    convection_amb_air = 11.1
    convection_amb_air_adjust = 12.0

    while math.fabs(convection_amb_air - convection_amb_air_adjust) > 0.005:
        convection_amb_air = convection_amb_air_adjust

        # INSULATION CONDUCTION
        sr3_radi_in = math.sqrt((SR3_area_ratio+1.0)*aperture_area / (6.0 * math.pi))
        sr3_radi_insul = sr3_radi_in + SR3_insul_thick
        SR3_length = 2.0 * sr3_radi_in
        C1 = temp_4
        C2 = math.log(sr3_radi_insul/sr3_radi_in) / (SR3_insul_cond * 2.0 * math.pi * SR3_length)
        
        # EVACUATED SPACE RADIATION
        sr3_radi_mainb = math.pi * sr3_radi_insul / (2.0 * math.sqrt(2.0))
        C3 = sigma * (2.0 * math.pi * SR3_length * sr3_radi_insul)
        C4 = 1.0 / SR3_insul_emiss
        C5 = (1.0-SR3_mainb_emiss) * (sr3_radi_insul/sr3_radi_mainb) / SR3_mainb_emiss
        
        # MAIN BODY CONDUCTION AND AMBIENT AIR CONVECTION
        sr3_radi_out = sr3_radi_mainb + SR3_mainb_thick
        A4 = 2.0 * math.pi * sr3_radi_out * SR3_length
        C6 = math.log(sr3_radi_out/sr3_radi_mainb) / (SR3_mainb_cond * 2.0 * math.pi * SR3_length)
        C7 = 1.0/(convection_amb_air * A4)
        C8 = temp_0
        
        C9 = C4 + C5
        C10 = C6 + C7
        
        inc = 1.0
        T2 = temp_0
        val1 = math.fabs( (C1-T2)/C2 - (C3*T2**4 - C3*( (C10*C1+C2*C8-C10*T2) /C2)**4)/C9)
        T2 = T2 + inc
        val2 = math.fabs( (C1-T2)/C2 - (C3*T2**4 - C3*( (C10*C1+C2*C8-C10*T2) /C2)**4)/C9)
        while val2 < val1:
            val1 = val2
            T2 = T2 + inc
            val2 = math.fabs( (C1-T2)/C2 - (C3*T2**4 - C3*( (C10*C1+C2*C8-C10*T2) /C2)**4)/C9)
        
        inc = -0.001
        val1 = val2
        T2 = T2 + inc
        val2 = math.fabs( (C1-T2)/C2 - (C3*T2**4 - C3*( (C10*C1+C2*C8-C10*T2) /C2)**4)/C9)
        while val2 < val1:
            val1 = val2
            T2 = T2 + inc
            val2 = math.fabs( (C1-T2)/C2 - (C3*T2**4 - C3*( (C10*C1+C2*C8-C10*T2) /C2)**4)/C9)
        
        T2 = T2 - inc
        T3 = (C10*C1 + C2*C8 - C10*T2) / C2
        T4 = T3 + C6 * (T2 - C1) / C2
        Tf_out = 0.5 * (T4 + temp_0)
        
        accel = 9.8
        diam = sr3_radi_out * 2.0
        visc = CP.PropsSI('VISCOSITY', 'T', Tf_out, 'P', p0air, 'Air')
        conduct = CP.PropsSI('CONDUCTIVITY', 'T', Tf_out, 'P', p0air, 'Air') 
        diff = conduct / (CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air') * CP.PropsSI('CPMASS', 'T', Tf_out, 'P', p0air, 'Air'))
        prand = CP.PropsSI('PRANDTL', 'T', Tf_out, 'P', p0air, 'Air')
        kvisc = visc / CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air')

        rayl = accel * (1.0/Tf_out) * (T4 - temp_0) * diam * diam * diam / (diff * kvisc)
        nuss = (0.6 + (0.387 * rayl**(1.0/6.0)) / ( ( 1.0 + (0.559/prand)**(9.0/16.0) )**(8.0/27.0) )) ** 2
        convection_amb_air_adjust = nuss * conduct / diam

    print(convection_amb_air)

    SR3_temps[0].append(temp_4-273.15)
    SR3_temps[1].append(T2-273.15)
    SR3_temps[2].append(T3-273.15)
    SR3_temps[3].append(T4-273.15)
    SR3_temps[4].append(temp_0-273.15)

    qloss = (C1-T2) / C2
    qloss = qloss * num_aper

    return qloss

def solve_SR3(dictionary, dni, SR3_temps):
    sigma = dictionary['sigma']
    emiss = dictionary['emiss']
    # UNUSED area_receiver = dictionary['area_receiver']
    dni = dni
    cp_abo3 = dictionary['cp_abo3']
    temp_0 = dictionary['temp_0']
    temp_1 = dictionary['temp_1']
    temp_4 = dictionary['temp_4']
    temp_13 = dictionary['temp_13']
    H_rxn = dictionary['H_rxn']
    cp_o2 = dictionary['cp_o2']
    mol_abo3_6 = dictionary['mol_abo3_6']
    mol_abo3_710 = dictionary['mol_abo3_710']
    A_sf = dictionary['A_sf']
    pSR3 = dictionary['pSR3']
    eta_sf = dictionary['eta_sf']
    mol_abo3_max = dictionary['mol_abo3_max']
    rho_abo3 = dictionary['rho_abo3']
    dens_abo3 = dictionary['dens_abo3']
    ullage = dictionary['ullage']
    delta = dictionary['delta']
    concentration_factor = float(dictionary['concentration_factor']) # MW/m^2
    
    cp_o2_5 = CP.PropsSI('CPMOLAR', 'T', temp_4, 'P', pSR3, 'Oxygen')
    volume_hs = (mol_abo3_max * rho_abo3 / dens_abo3) * (1.0 + ullage)
    
    concentration_factor = concentration_factor * (10.0**6.0) # W/m^2
    design_point = 900.0 # W/m^2
    aperture_area = A_sf * design_point * eta_sf / concentration_factor

    if dni > 350.0: # enough solar radiation
        energy_2 = dni * A_sf * eta_sf
        energy_18 = sigma * emiss * aperture_area * (temp_4**4)
        energy_18 += solve_SR3_qloss(dictionary, aperture_area, SR3_temps)
        mol_abo3_1413 = (energy_2-energy_18) / (cp_abo3 * (temp_4-temp_1) + delta* 0.5 * H_rxn + delta * 0.5 * cp_o2_5 * (temp_4-temp_0))

        if ((mol_abo3_1413 - mol_abo3_710) * 3600.0 + mol_abo3_6) * rho_abo3/dens_abo3 * (1.0+ullage) > volume_hs:
            mol_abo3_1413 = mol_abo3_710 + (((volume_hs/(1.0+ullage)) * dens_abo3/rho_abo3 - mol_abo3_6) / 3600.0) - 0.001
            energy_2 = ( mol_abo3_1413 * (cp_abo3 * (temp_4-temp_1) + delta * 0.5 * H_rxn + delta * 0.5 * cp_o2_5 * (temp_4-temp_0)) ) + energy_18
    else: # not enough solar radiation
        energy_2 = 0.0
        energy_18 = 0.0
        mol_abo3_1413 = 0.0

    mol_o2 = mol_abo3_1413 * delta * 0.5

    energy_1 = mol_abo3_1413 * cp_abo3 * (temp_1-temp_0)
    energy_4 = mol_abo3_1413 * cp_abo3 * (temp_4-temp_0) + mol_abo3_1413 * delta * 0.5 * H_rxn
    energy_5 = mol_o2 * cp_o2_5 * (temp_4-temp_0)
    energy_13 = mol_abo3_1413*cp_abo3*(temp_13-temp_0)

    if math.fabs(energy_1 + energy_2 - energy_4 - energy_5 - energy_18) > 1.0:
        print('SR3 energy imbalance warning {!s}'.format(energy_1 + energy_2 - energy_4 - energy_5 - energy_18))

    dictionary['mol_abo3_1413'] = mol_abo3_1413
    dictionary['mol_o2'] = mol_o2
    dictionary['energy_1'] = energy_1
    dictionary['energy_2'] = energy_2
    dictionary['energy_4'] = energy_4
    dictionary['energy_5'] = energy_5
    dictionary['energy_13'] = energy_13
    dictionary['energy_18'] = energy_18

    return 0

def solve_ROx_qloss(dictionary, ROx_temps):
    temp_0 = dictionary['temp_0']
    temp_7 = dictionary['temp_7']
    temp_9 = dictionary['temp_9']
    temp_11 = dictionary['temp_11']
    temp_10 = dictionary['temp_10']
    p0air = dictionary['p0air']
    pr_turbine = dictionary['pr_turbine']
    ROx_insul_thick = dictionary['ROx_insul_thick']
    ROx_insul_mult = dictionary['ROx_insul_mult'] # UNUSED
    ROx_insul_cond = dictionary['ROx_insul_cond']
    L_ROx = dictionary['L_ROx']
    D_ROx = dictionary['D_ROx']
    m_dot_air = dictionary['mass_flow_air'] * 0.001 / dictionary['ROx_pipes']
    i = dictionary['i']
    
    r_ROx_in = D_ROx * 0.5
    r_ROx_out = r_ROx_in + ROx_insul_thick
    A_ROx_in = 2.0 * math.pi * r_ROx_in * L_ROx
    A_ROx_out = 2.0 * math.pi * r_ROx_out * L_ROx

    temp_inf_ROx = 0.25 * (temp_7 + temp_9 + temp_10 + temp_11)
    convection_air_in = 4.0
    convection_air_in_adjust = 5.0
    convection_air_out = 5.0
    convection_air_out_adjust = 6.0

    while math.fabs(convection_air_in - convection_air_in_adjust) > 0.005 or math.fabs(convection_air_out - convection_air_out_adjust) > 0.005:
        convection_air_in = convection_air_in_adjust
        convection_air_out = convection_air_out_adjust

        R_tot = (convection_air_in * A_ROx_in)**-1
        R_tot += math.log(r_ROx_out/r_ROx_in)/(ROx_insul_cond * 2.0 * math.pi * L_ROx)
        R_tot += (convection_air_out * A_ROx_out)**-1

        q_loss = (temp_inf_ROx - temp_0) / R_tot

        Ts_in = temp_inf_ROx - q_loss * ((convection_air_in * A_ROx_in)**-1)
        Ts_out = temp_inf_ROx - q_loss * ((convection_air_in * A_ROx_in)**-1 + math.log(r_ROx_out/r_ROx_in)/(ROx_insul_cond * 2.0 * math.pi * L_ROx))
        Tf_in = 0.5 * (temp_inf_ROx + Ts_in)
        Tf_out = 0.5 * (Ts_out + temp_0)
        
        accel = 9.8 # m/s^2
        visc_air_in = CP.PropsSI('VISCOSITY', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air')
        cond_air_in = CP.PropsSI('CONDUCTIVITY', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air')
        diff_air_in = cond_air_in / (CP.PropsSI('DMASS', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air') * CP.PropsSI('CPMASS', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air'))
        prand_air_in = CP.PropsSI('PRANDTL', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air')
        kvisc_air_in = visc_air_in / CP.PropsSI('DMASS', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air')

        visc_air_out = CP.PropsSI('VISCOSITY', 'T', Tf_out, 'P', p0air, 'Air')
        cond_air_out = CP.PropsSI('CONDUCTIVITY', 'T', Tf_out, 'P', p0air, 'Air')
        diff_air_out = cond_air_out / (CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air') * CP.PropsSI('CPMASS', 'T', Tf_out, 'P', p0air, 'Air'))
        prand_air_out = CP.PropsSI('PRANDTL', 'T', Tf_out, 'P', p0air, 'Air')
        kvisc_air_out = visc_air_out / CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air')

        Ray_air_in = 4.0 * m_dot_air / (math.pi * D_ROx * visc_air_in)
        nuss_in = 0.023 * Ray_air_in**(4.0/5.0) * prand_air_in**(0.3)
        convection_air_in_adjust = (cond_air_in / D_ROx) * nuss_in

        Ray_air_out = (accel * (1.0/Tf_out) * (Ts_out - temp_0) * L_ROx**3) / (kvisc_air_out * diff_air_out)
        nuss_out = (0.825 + ( (0.387 * Ray_air_out**(1.0/6.0)) / ((1.0 + (0.492/prand_air_out)**(9.0/16.0))**(8.0/27.0)) ))**2
        convection_air_out_adjust = (cond_air_out / L_ROx) * nuss_out

    ROx_temps[0].append(temp_inf_ROx-273.15)
    ROx_temps[1].append(Ts_in-273.15)
    ROx_temps[2].append(Ts_out-273.15)
    ROx_temps[3].append(temp_0-273.15)

    return q_loss

def solve_ROx(dictionary, ROx_temps, dispatch):
    mol_abo3_6 = dictionary['mol_abo3_6']
    mol_air = dictionary['mol_air']
    cp_air = dictionary['cp_air']
    cp_abo3 = dictionary['cp_abo3']
    temp_0 = dictionary['temp_0']
    temp_7 = dictionary['temp_7']
    temp_9 = dictionary['temp_9']
    temp_11 = dictionary['temp_11']
    temp_10 = dictionary['temp_10']
    p0air = dictionary['p0air']
    pr_turbine = dictionary['pr_turbine']
    delta = dictionary['delta']
    H_rxn = dictionary['H_rxn']
    time_step = dictionary['time_step']
    L_ROx = dictionary['L_ROx']
    D_ROx = dictionary['D_ROx']
    air_mol_mass = dictionary['molar_mass_air'] * 0.001 #kg/mol
    part_mol_mass = dictionary['mol_mass_abo3'] * 0.001 #kg/mol
    part_mol_vol = dictionary['rho_abo3'] #m^3/mol
    part_pack_dens = dictionary['part_pack_dens'] #%
    D_ROx = dictionary['D_ROx']
    L_ROx = dictionary['L_ROx']
    rox_pipes = dictionary['ROx_pipes']
    time_res = dictionary['time_res']
    i = dictionary['i']

    if temp_10 < temp_9:
        print("ERROR: temp_10 is less than temp_9\n")
        exit(0)
    temp_10_2 = 200.0 + 273.15 # CS temp inlet. Used to keep molar flow of system low
    temp_inf_ROx = 0.25 * (temp_7 + temp_9 + temp_10 + temp_11)

    cp_air = CP.PropsSI('CPMOLAR', 'T', 0.5*(temp_9+temp_11), 'P', p0air*pr_turbine, 'Air')
    visc = CP.PropsSI('VISCOSITY', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')
    prand = CP.PropsSI('PRANDTL', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')
    cond = CP.PropsSI('CONDUCTIVITY', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')
    air_dens = CP.PropsSI('DMASS', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air') #kg/m^3
    kvisc = visc / air_dens

    D_p = 0.000250 #m
    A_p = math.pi * ((D_p*0.5)**2) #m^2
    V_p = 4.0*math.pi*((D_p*0.5)**3)/3.0 #m^3

    accel = 9.8
    part_mass = V_p * part_mol_mass / part_mol_vol #kg
    drag_c = 0.5

    energy_9 = mol_air * cp_air * (temp_9-temp_0)
    energy_11 = mol_air * cp_air * (temp_11-temp_0)

    energy_17 = solve_ROx_qloss(dictionary, ROx_temps)
    energy_17 = rox_pipes * energy_17
    mol_abo3_710 = (energy_11 + energy_17 - energy_9) / (cp_abo3*(temp_7-temp_10) + delta * 0.5 * H_rxn)

    if dispatch[i] == 1 and mol_abo3_710 <= mol_abo3_6 / time_step and mol_abo3_710 > 0:
        mol_abo3_710 = mol_abo3_710 / rox_pipes
        mol_air = mol_air / rox_pipes

        m_dot_pipe_part = mol_abo3_710 * part_mol_mass #kg/s
        m_dot_pipe_air = mol_air * air_mol_mass #kg/s

        max_vel = math.sqrt(2.0*part_mass*accel / (air_dens*A_p*drag_c))
        part_dens = part_mol_mass * part_pack_dens / part_mol_vol
        dens = air_dens * (1.0 - part_pack_dens)
        pipe_vol_dot_part = m_dot_pipe_part/part_dens
        pipe_vol_dot_air = m_dot_pipe_air/dens
        part_pack_vol = pipe_vol_dot_part / (pipe_vol_dot_air + pipe_vol_dot_part)
        pipe_cyl_cross_area = (pipe_vol_dot_part + pipe_vol_dot_air) / max_vel
        pipe_cyl_diam = 2.0 * math.sqrt(pipe_cyl_cross_area/math.pi)

        if pipe_cyl_diam - D_ROx > 0.0005:
            print('pipe_cyl_diam greater by ', pipe_cyl_diam - D_ROx, ' which is ', 100.0*(pipe_cyl_diam - D_ROx)/D_ROx, '%')
            print('Temp 7 is', temp_7-273.15, 'C')
        elif pipe_cyl_diam - D_ROx < -0.1:
            print('D_ROx may be sized too long!', pipe_cyl_diam, D_ROx)

        Reyn = max_vel * D_p / kvisc
        Nuss = Nuss = 2.0 + 0.6 * (Reyn**0.5) * (prand**(1.0/3.0))
        U = Nuss * cond / D_p

        pipe_cyl_length = time_res * max_vel * 0.5

        if pipe_cyl_length - L_ROx > 0.005:
            print('pipe_cyl_length greater by ', pipe_cyl_length - L_ROx, ' which is ', 100.0*(pipe_cyl_length - L_ROx)/L_ROx, '%')
            print('Temp 7 is', temp_7-273.15, 'C')
        elif pipe_cyl_length - L_ROx < -0.1:
            print('L_ROx may be sized too long!')

    if dispatch[i] != 1 or mol_abo3_710 > mol_abo3_6 / time_step or mol_abo3_710 < 0:
        mol_abo3_710 = 0.0
        mol_air = 0.0
        energy_9 = 0.0
        energy_11 = 0.0
        energy_17 = 0.0
        L_ROx = 0.0
        D_ROx = 0.0
        U = 0.0
        ROx_HX_eff = 0.0
        NTU = 0.0

    mol_abo3_710 = mol_abo3_710 * rox_pipes
    mol_air = mol_air * rox_pipes

    energy_7 = mol_abo3_710 * cp_abo3 * (temp_7-temp_0) + mol_abo3_710 * delta * 0.5 * H_rxn
    energy_10 = mol_abo3_710 * cp_abo3 * (temp_10-temp_0)

    if math.fabs(energy_7 + energy_9 - energy_10 - energy_11 - energy_17) > 1.0:
        print('ROx energy imbalance')

    dictionary['mol_air'] = mol_air
    dictionary['mol_abo3_710'] = mol_abo3_710
    dictionary['temp_10'] = temp_10
    dictionary['temp_10_2'] = temp_10_2
    dictionary['energy_7'] = energy_7
    dictionary['energy_9'] = energy_9
    dictionary['energy_10'] = energy_10
    dictionary['energy_11'] = energy_11
    dictionary['energy_17'] = energy_17

    return 0

def calibrate_ROx(dictionary):
    mol_abo3_6 = dictionary['mol_abo3_6']
    mol_air = dictionary['mol_air']
    cp_air = dictionary['cp_air']
    cp_abo3 = dictionary['cp_abo3']
    temp_0 = dictionary['temp_0']
    temp_7 = dictionary['temp_7']
    temp_9 = dictionary['temp_9']
    temp_11 = dictionary['temp_11']
    temp_10 = dictionary['temp_10']
    delta = dictionary['delta']
    H_rxn = dictionary['H_rxn']
    time_step = dictionary['time_step']
    p0air = dictionary['p0air']
    pr_turbine = dictionary['pr_turbine']
    ROx_insul_thick = dictionary['ROx_insul_thick']
    ROx_insul_mult = dictionary['ROx_insul_mult'] # UNUSED
    ROx_insul_cond = dictionary['ROx_insul_cond']
    L_ROx = dictionary['L_ROx']
    D_ROx = dictionary['D_ROx']
    rox_pipes = dictionary['ROx_pipes']
    air_mol_mass = dictionary['molar_mass_air'] * 0.001 #kg/mol
    part_mol_mass = dictionary['mol_mass_abo3'] * 0.001 #kg/mol
    part_mol_vol = dictionary['rho_abo3'] #m^3/mol
    part_pack_dens = dictionary['part_pack_dens'] #%
    time_res = dictionary['time_res']

    cp_air = CP.PropsSI('CPMOLAR', 'T', 0.5*(temp_9+temp_11), 'P', p0air*pr_turbine, 'Air')
    if temp_10 < temp_9:
        print("ERROR: temp_10 is less than temp_9\n")
        exit(0)
    temp_10_2 = 200.0 + 273.15 # CS temp inlet. Used to keep molar flow of system low
    temp_inf_ROx = 0.25 * (temp_7 + temp_9 + temp_10 + temp_11)

    D_p = 0.000250 #m
    A_p = math.pi * ((D_p*0.5)**2) #m^2
    V_p = 4.0*math.pi*((D_p*0.5)**3)/3.0 #m^3
    accel = 9.8 # m/s^2
    drag_c = 0.5
    part_mass = V_p * part_mol_mass / part_mol_vol#kg

    visc = CP.PropsSI('VISCOSITY', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')
    prand = CP.PropsSI('PRANDTL', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')
    cond = CP.PropsSI('CONDUCTIVITY', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')
    air_dens = CP.PropsSI('DMASS', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air') #kg/m^3
    kvisc = visc / air_dens

    visc_air_in = CP.PropsSI('VISCOSITY', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')
    cond_air_in = CP.PropsSI('CONDUCTIVITY', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')
    diff_air_in = cond_air_in / (CP.PropsSI('DMASS', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air') * CP.PropsSI('CPMASS', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air'))
    prand_air_in = CP.PropsSI('PRANDTL', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')
    kvisc_air_in = visc_air_in / CP.PropsSI('DMASS', 'T', temp_inf_ROx, 'P', p0air*pr_turbine, 'Air')

    visc_air_out = CP.PropsSI('VISCOSITY', 'T', temp_0, 'P', p0air, 'Air')
    cond_air_out = CP.PropsSI('CONDUCTIVITY', 'T', temp_0, 'P', p0air, 'Air')
    diff_air_out = cond_air_out / (CP.PropsSI('DMASS', 'T', temp_0, 'P', p0air, 'Air') * CP.PropsSI('CPMASS', 'T', temp_0, 'P', p0air, 'Air'))
    prand_air_out = CP.PropsSI('PRANDTL', 'T', temp_0, 'P', p0air, 'Air')
    kvisc_air_out = visc_air_out / CP.PropsSI('DMASS', 'T', temp_0, 'P', p0air, 'Air')

    mol_abo3_710 = 0.0
    mol_abo3_710_old = 1000.0

    while math.fabs(mol_abo3_710 - mol_abo3_710_old) > 0.001:
        mol_abo3_710 = mol_abo3_710_old

        energy_7 = mol_abo3_710 * cp_abo3 * (temp_7-temp_0) + mol_abo3_710 * 0.5 * delta * H_rxn
        energy_9 = mol_air * cp_air * (temp_9-temp_0)
        energy_10 = mol_abo3_710 * cp_abo3 * (temp_10-temp_0)
        energy_11 = mol_air * cp_air * (temp_11-temp_0)

        mol_abo3_710 = mol_abo3_710 / rox_pipes
        mol_air = mol_air / rox_pipes
        m_dot_pipe_part = mol_abo3_710 * part_mol_mass #kg/s
        m_dot_pipe_air = mol_air * air_mol_mass #kg/s

        max_vel = math.sqrt(2.0*part_mass*accel / (air_dens*A_p*drag_c))
        print('max velocity', max_vel)
        part_dens = part_mol_mass * part_pack_dens / part_mol_vol
        dens = air_dens * (1.0 - part_pack_dens)
        pipe_vol_dot_part = m_dot_pipe_part/part_dens
        pipe_vol_dot_air = m_dot_pipe_air/dens
        part_pack_vol = pipe_vol_dot_part / (pipe_vol_dot_air + pipe_vol_dot_part)
        pipe_cyl_cross_area = (pipe_vol_dot_part + pipe_vol_dot_air) / max_vel
        pipe_cyl_diam = 2.0 * math.sqrt(pipe_cyl_cross_area/math.pi)

        Reyn = max_vel * D_p / kvisc
        Nuss = 2.0 + 0.6 * (Reyn**0.5) * (prand**(1.0/3.0))
        U = Nuss * cond / D_p

        pipe_cyl_length = time_res * max_vel * 0.5

        r_ROx = pipe_cyl_diam * 0.5
        r_ROx_out = r_ROx + ROx_insul_thick
        A_ROx_in = 2.0 * math.pi * r_ROx * pipe_cyl_length
        A_ROx_out = 2.0 * math.pi * r_ROx_out * pipe_cyl_length

        Ray_air_in = 4.0 * m_dot_pipe_air / (math.pi * pipe_cyl_diam * visc_air_in)
        nuss_in = 0.023 * Ray_air_in**(4.0/5.0) * prand_air_in**(0.3)
        h_ROx_in = (cond_air_in / pipe_cyl_diam) * nuss_in

        Ray_air_out = (accel * 1.0/(temp_0+12.5) * ((temp_0+25.0)-temp_0) * pipe_cyl_length**3) / (kvisc_air_out * diff_air_out)
        nuss_out = (0.825 + ( (0.387 * Ray_air_out**(1.0/6.0)) / ((1.0 + (0.492/prand_air_out)**(9.0/16.0))**(8.0/27.0)) ))**2
        h_ROx_out = (cond_air_out / pipe_cyl_length) * nuss_out

        R_tot = (h_ROx_in * A_ROx_in)**-1 + math.log(r_ROx_out/r_ROx)/(ROx_insul_cond * 2.0 * math.pi * pipe_cyl_length) + (h_ROx_out * A_ROx_out)**-1
        energy_17 = (temp_inf_ROx - temp_0) / R_tot
        energy_17 = rox_pipes * energy_17

        mol_abo3_710_old = (energy_11 + energy_17 - energy_9) / (cp_abo3*(temp_7-temp_10) + delta * 0.5 * H_rxn)
        mol_air = mol_air * rox_pipes
        mol_abo3_710 = mol_abo3_710 * rox_pipes

    qloss = energy_17 / rox_pipes
    Ts_in = temp_inf_ROx - qloss * ((h_ROx_in * A_ROx_in)**-1)
    Ts_out = temp_inf_ROx - qloss * ((h_ROx_in * A_ROx_in)**-1 + math.log(r_ROx_out/r_ROx)/(ROx_insul_cond * 2.0 * math.pi * pipe_cyl_length))
    Tf_in = 0.5 * (temp_inf_ROx + Ts_in)
    Tf_out = 0.5 * (Ts_out + temp_0)

    visc_air_in = CP.PropsSI('VISCOSITY', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air')
    cond_air_in = CP.PropsSI('CONDUCTIVITY', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air')
    diff_air_in = cond_air_in / (CP.PropsSI('DMASS', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air') * CP.PropsSI('CPMASS', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air'))
    prand_air_in = CP.PropsSI('PRANDTL', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air')
    kvisc_air_in = visc_air_in / CP.PropsSI('DMASS', 'T', Tf_in, 'P', p0air*pr_turbine, 'Air')

    visc_air_out = CP.PropsSI('VISCOSITY', 'T', Tf_out, 'P', p0air, 'Air')
    cond_air_out = CP.PropsSI('CONDUCTIVITY', 'T', Tf_out, 'P', p0air, 'Air')
    diff_air_out = cond_air_out / (CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air') * CP.PropsSI('CPMASS', 'T', Tf_out, 'P', p0air, 'Air'))
    prand_air_out = CP.PropsSI('PRANDTL', 'T', Tf_out, 'P', p0air, 'Air')
    kvisc_air_out = visc_air_out / CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air')

    mol_abo3_710 = 0.0
    while math.fabs(mol_abo3_710 - mol_abo3_710_old) > 0.001:
        mol_abo3_710 = mol_abo3_710_old

        energy_7 = mol_abo3_710 * cp_abo3 * (temp_7-temp_0) + mol_abo3_710 * 0.5 * delta * H_rxn
        energy_9 = mol_air * cp_air * (temp_9-temp_0)
        energy_10 = mol_abo3_710 * cp_abo3 * (temp_10-temp_0)
        energy_11 = mol_air * cp_air * (temp_11-temp_0)

        mol_abo3_710 = mol_abo3_710 / rox_pipes
        mol_air = mol_air / rox_pipes
        m_dot_pipe_part = mol_abo3_710 * part_mol_mass #kg/s
        m_dot_pipe_air = mol_air * air_mol_mass #kg/s

        max_vel = math.sqrt(2.0*part_mass*accel / (air_dens*A_p*drag_c))
        part_dens = part_mol_mass * part_pack_dens / part_mol_vol
        dens = air_dens * (1.0 - part_pack_dens)
        pipe_vol_dot_part = m_dot_pipe_part/part_dens
        pipe_vol_dot_air = m_dot_pipe_air/dens
        part_pack_vol = pipe_vol_dot_part / (pipe_vol_dot_air + pipe_vol_dot_part)
        pipe_cyl_cross_area = (pipe_vol_dot_part + pipe_vol_dot_air) / max_vel
        pipe_cyl_diam = 2.0 * math.sqrt(pipe_cyl_cross_area/math.pi)

        Reyn = max_vel * D_p / kvisc
        Nuss = 2.0 + 0.6 * (Reyn**0.5) * (prand**(1.0/3.0))
        U = Nuss * cond / D_p

        pipe_cyl_length = time_res * max_vel * 0.5

        r_ROx = pipe_cyl_diam * 0.5
        r_ROx_out = r_ROx + ROx_insul_thick
        A_ROx_in = 2.0 * math.pi * r_ROx * pipe_cyl_length
        A_ROx_out = 2.0 * math.pi * r_ROx_out * pipe_cyl_length

        Ray_air_in = 4.0 * m_dot_pipe_air / (math.pi * pipe_cyl_diam * visc_air_in)
        nuss_in = 0.023 * Ray_air_in**(4.0/5.0) * prand_air_in**(0.3)
        h_ROx_in = (cond_air_in / pipe_cyl_diam) * nuss_in

        Ray_air_out = (accel * 1.0/(Tf_out) * (Ts_out - temp_0) * pipe_cyl_length**3) / (kvisc_air_out * diff_air_out)
        nuss_out = (0.825 + ( (0.387 * Ray_air_out**(1.0/6.0)) / ((1.0 + (0.492/prand_air_out)**(9.0/16.0))**(8.0/27.0)) ))**2
        h_ROx_out = (cond_air_out / pipe_cyl_length) * nuss_out

        R_tot = (h_ROx_in * A_ROx_in)**-1 + math.log(r_ROx_out/r_ROx)/(ROx_insul_cond * 2.0 * math.pi * pipe_cyl_length) + (h_ROx_out * A_ROx_out)**-1
        energy_17 = (temp_inf_ROx - temp_0) / R_tot
        energy_17 = rox_pipes * energy_17

        mol_abo3_710_old = (energy_11 + energy_17 - energy_9) / (cp_abo3*(temp_7-temp_10) + delta * 0.5 * H_rxn)
        mol_air = mol_air * rox_pipes
        mol_abo3_710 = mol_abo3_710 * rox_pipes
    mol_abo3_710 = mol_abo3_710_old

    if mol_abo3_710 > mol_abo3_6 / time_step or mol_abo3_710 < 0:
        mol_abo3_710 = 0.0
        mol_air = 0.0
        energy_9 = 0.0
        energy_11 = 0.0
        energy_17 = 0.0
        L_ROx = 0.0
        D_ROx = 0.0
        U = 0.0
        ROx_HX_eff = 0.0
        NTU = 0.0

    energy_7 = mol_abo3_710 * cp_abo3 * (temp_7-temp_0) + mol_abo3_710 * delta * 0.5 * H_rxn
    energy_9 = mol_air * cp_air * (temp_9-temp_0)
    energy_10 = mol_abo3_710 * cp_abo3 * (temp_10-temp_0)
    energy_11 = mol_air * cp_air * (temp_11-temp_0)

    dictionary['mol_air'] = mol_air
    dictionary['mol_abo3_710'] = mol_abo3_710
    dictionary['temp_10'] = temp_10
    dictionary['temp_10_2'] = temp_10_2
    dictionary['energy_7'] = energy_7
    dictionary['energy_9'] = energy_9
    dictionary['energy_10'] = energy_10
    dictionary['energy_11'] = energy_11
    dictionary['energy_17'] = energy_17
    dictionary['L_ROx'] = pipe_cyl_length
    dictionary['D_ROx'] = pipe_cyl_diam

    return 0

def solve_HS_qloss(dictionary, r_sb, r_sb_out, ht_hs_in, D_sb_out, ht_hs_out, A_circle_h, A_circle_c, A_side_c, HS_temps):
    HS_insul_thick = dictionary['HS_insul_thick']
    HS_insul_mult = dictionary['HS_insul_mult']
    HS_insul_cond = dictionary['HS_insul_cond']
    p0air = dictionary['p0air']
    temp_6 = dictionary['temp_6']
    temp_0 = dictionary['temp_0']
    i = dictionary['i']

    convection_up = 2.0
    convection_side = 2.0
    convection_low = 2.0
    convection_up_adjust = 3.0
    convection_side_adjust = 3.0
    convection_low_adjust = 3.0

    while math.fabs(convection_up - convection_up_adjust) > 0.005 or math.fabs(convection_side - convection_side_adjust) > 0.005 or math.fabs(convection_low - convection_low_adjust) > 0.005:
        convection_up = convection_up_adjust
        convection_side = convection_side_adjust
        convection_low = convection_low_adjust

        R_up = (HS_insul_thick / (HS_insul_cond*A_circle_h)) + (1.0/(convection_up*A_circle_c))
        R_side = (math.log(r_sb_out/r_sb) / (HS_insul_cond * 2.0 * math.pi * ht_hs_in)) + (1.0/(convection_side*A_side_c))
        R_low = (HS_insul_thick / (HS_insul_cond*A_circle_h)) + (1.0/(convection_low*A_circle_c))

        R_tot = ((R_up**-1) + (R_side**-1) + (R_low**-1))**-1

        qloss = (temp_6 - temp_0) / R_tot

        Ts_in = temp_6
        Ts_out = temp_0 + qloss * 1.0 / (convection_side*A_side_c)
        Tf_in = 0.5 * (temp_6 + Ts_in)
        Tf_out = 0.5 * (temp_0 + Ts_out)

        accel = 9.8 # m/s^2

        visc_air = CP.PropsSI('VISCOSITY', 'T', Tf_out, 'P', p0air, 'Air')
        cond_air = CP.PropsSI('CONDUCTIVITY', 'T', Tf_out, 'P', p0air, 'Air')
        diff_air = cond_air / (CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air') * CP.PropsSI('CPMASS', 'T', Tf_out, 'P', p0air, 'Air'))
        prand_air = CP.PropsSI('PRANDTL', 'T', Tf_out, 'P', p0air, 'Air')
        kvisc_air = visc_air / CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air')

        Ray_air_up = (accel * (1.0/Tf_out) * (Ts_out - temp_0) * D_sb_out**3) / (kvisc_air * diff_air)
        nuss_up = 0.54 * (Ray_air_up**0.25)
        convection_up_adjust = (cond_air / D_sb_out) * nuss_up

        Ray_air_side = (accel * (1.0/Tf_out) * (Ts_out - temp_0) * ht_hs_out**3) / (kvisc_air * diff_air)
        nuss_side = (0.825 + ( (0.387 * Ray_air_side**(1.0/6.0)) / ((1.0 + ((0.492/prand_air)**(9.0/16.0)) )**(8.0/27.0)) ))**2
        convection_side_adjust = (cond_air / ht_hs_out) * nuss_side

        Ray_air_low = (accel * (1.0/Tf_out) * (Ts_out - temp_0) * D_sb_out**3) / (kvisc_air * diff_air)
        nuss_low = 0.52 * (Ray_air_low**0.2)
        convection_low_adjust = (cond_air / D_sb_out) * nuss_low

    HS_temps[0].append(temp_6-273.15)
    HS_temps[1].append(Ts_in-273.15)
    HS_temps[2].append(Ts_out-273.15)
    HS_temps[3].append(temp_0-273.15)

    return qloss

def solve_HS(dictionary, mol_abo3_max, HS_temps):
    hd = dictionary['hd']
    mol_N_6 = dictionary['mol_N_6']
    mol_abo3_710 = dictionary['mol_abo3_710']
    mol_abo3_1413 = dictionary['mol_abo3_1413']
    mol_abo3_6 = dictionary['mol_abo3_6']
    cp_abo3 = dictionary['cp_abo3']
    temp_0 = dictionary['temp_0']
    temp_4 = dictionary['temp_4']
    temp_6 = dictionary['temp_6']
    temp_7 = dictionary['temp_7']
    time_step = float(dictionary['time_step'])
    delta = dictionary['delta']
    H_rxn = dictionary['H_rxn']
    mol_abo3_max = dictionary['mol_abo3_max']
    rho_abo3 = dictionary['rho_abo3']
    dens_abo3 = dictionary['dens_abo3']
    ullage = dictionary['ullage']
    p0air = dictionary['p0air']
    R = dictionary['R']
    HS_insul_thick = dictionary['HS_insul_thick']
    energy_6 = dictionary['energy_6']
    energy_4 = dictionary['energy_4']
    energy_7 = dictionary['energy_7']
    cp_N = CP.PropsSI('CPMOLAR', 'T', temp_6, 'P', p0air, 'Nitrogen')

    mol_abo3_6last = mol_abo3_6
    mol_abo3_6 = (mol_abo3_1413 - mol_abo3_710) * time_step + mol_abo3_6last
    volume_abo3_6_avg = (mol_abo3_6last + mol_abo3_6) * 0.5 * rho_abo3 / dens_abo3

    volume_hs = mol_abo3_max * rho_abo3 / dens_abo3 * (1.0 + ullage)
    r_sb = (volume_hs / (math.pi * 2.0 * hd))**(1.0/3.0)
    ht_hs_in = volume_hs / (math.pi * r_sb**2)
    A_hs = ht_hs_in * 2.0 * math.pi * r_sb + 2.0 * math.pi * r_sb**2

    volume_N_6 = volume_hs - mol_abo3_6 * rho_abo3 / dens_abo3
    mol_N_6last = mol_N_6
    mol_N_6 = (volume_N_6 * p0air) / (R * temp_6)

    if mol_N_6 < mol_N_6last:
        mol_N_22 = 0.0
        mol_N_23 = (mol_N_6last - mol_N_6) / time_step
    elif mol_N_6 == mol_N_6last:
        mol_N_22 = 0.0
        mol_N_23 = 0.0
    elif mol_N_6 > mol_N_6last:
        mol_N_23 = 0.0
        mol_N_22 = (mol_N_6-mol_N_6last) / time_step

    r_sb_out = r_sb + HS_insul_thick
    D_sb_out = r_sb_out * 2.0
    ht_hs_out = ht_hs_in + 2.0 * HS_insul_thick

    A_circle_h = math.pi * (r_sb**2)
    A_circle_c = math.pi * (r_sb_out**2)
    A_side_c = math.pi * D_sb_out * ht_hs_out
    
    energy_8 = solve_HS_qloss(dictionary, r_sb, r_sb_out, D_sb_out, ht_hs_in, ht_hs_out, A_circle_h, A_circle_c, A_side_c, HS_temps)
    energy_23 = mol_N_23*cp_N*(temp_6-temp_0)
    energy_22 = 0.0

    numerator = energy_6 + energy_4 - energy_7 - energy_23 - energy_8 - (mol_abo3_6*delta*0.5*H_rxn)/time_step
    denominator = (mol_abo3_6*cp_abo3 + mol_N_6*cp_N)/time_step
    temp_6_new = temp_0 + numerator/denominator

    if temp_6_new < temp_0:
        energy_8 = energy_6 + energy_4 - energy_7 - energy_23 - (mol_abo3_6*delta*0.5*H_rxn)/time_step - 0.001
        numerator = energy_6 + energy_4 - energy_7 - energy_23 - energy_8 - (mol_abo3_6*delta*0.5*H_rxn)/time_step
        denominator = (mol_abo3_6*cp_abo3 + mol_N_6*cp_N)/time_step
        temp_6_new = temp_0 + numerator/denominator

    energy_6_old = energy_6
    energy_6 = (mol_abo3_6*cp_abo3*(temp_6_new-temp_0) + mol_abo3_6*0.5*delta*H_rxn + mol_N_6*cp_N*(temp_6_new-temp_0))/time_step

    hserror = energy_6_old + energy_4 - energy_7 - energy_23 - energy_8 - energy_6
    if hserror > 2.0:
        print('HS error is ', hserror, 'W')
    elif hserror < -2.0:
        print('HS error is ', hserror, 'W')

    dictionary['mol_N_6'] = mol_N_6
    dictionary['energy_8'] = energy_8
    dictionary['mol_abo3_6'] = mol_abo3_6 
    dictionary['energy_6'] = energy_6 
    dictionary['mol_abo3_max'] = mol_abo3_max 
    dictionary['mol_N_22'] = mol_N_22 
    dictionary['mol_N_23'] = mol_N_23
    dictionary['energy_22'] = energy_22
    dictionary['energy_23'] = energy_23
    dictionary['temp_6'] = temp_6_new
    
    return r_sb, A_hs, volume_hs

def solve_CS_qloss(dictionary, r_cs, r_cs_out, D_cs_out, ht_cs_in, ht_cs_out, A_circle_h, A_circle_c, A_side_c, CS_temps):
    CS_insul_thick = dictionary['CS_insul_thick']
    CS_insul_mult = dictionary['CS_insul_mult']
    CS_insul_cond = dictionary['CS_insul_cond']
    p0air = dictionary['p0air']
    temp_12 = dictionary['temp_12']
    temp_0 = dictionary['temp_0']
    i = dictionary['i']

    convection_up = 2.0
    convection_side = 2.0
    convection_low = 2.0
    convection_up_adjust = 3.0
    convection_side_adjust = 3.0
    convection_low_adjust = 3.0

    while math.fabs(convection_up - convection_up_adjust) > 0.005 or math.fabs(convection_side - convection_side_adjust) > 0.005 or math.fabs(convection_low - convection_low_adjust) > 0.005:
        convection_up = convection_up_adjust
        convection_side = convection_side_adjust
        convection_low = convection_low_adjust

        R_up = (CS_insul_thick / (CS_insul_cond*A_circle_h)) + (1.0/(convection_up*A_circle_c))
        R_side = (math.log(r_cs_out/r_cs) / (CS_insul_cond * 2.0 * math.pi * ht_cs_in)) + (1.0/(convection_side*A_side_c))
        R_low = (CS_insul_thick / (CS_insul_cond*A_circle_h)) + (1.0/(convection_low*A_circle_c))

        R_tot = ((R_up**-1) + (R_side**-1) + (R_low**-1))**-1

        qloss = (temp_12 - temp_0) / R_tot

        Ts_in = temp_12
        Ts_out = temp_0 + qloss * 1.0 / (convection_side*A_side_c) 
        Tf_in = 0.5 * (temp_12 + Ts_in)
        Tf_out = 0.5 * (temp_0 + Ts_out)

        accel = 9.8 # m/s^2
    
        visc_air = CP.PropsSI('VISCOSITY', 'T', Tf_out, 'P', p0air, 'Air')
        cond_air = CP.PropsSI('CONDUCTIVITY', 'T', Tf_out, 'P', p0air, 'Air')
        diff_air = cond_air / (CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air') * CP.PropsSI('CPMASS', 'T', Tf_out, 'P', p0air, 'Air'))
        prand_air = CP.PropsSI('PRANDTL', 'T', Tf_out, 'P', p0air, 'Air')
        kvisc_air = visc_air / CP.PropsSI('DMASS', 'T', Tf_out, 'P', p0air, 'Air')

        Ray_air_up = (accel * (1.0/Tf_out) * (Ts_out - temp_0) * D_cs_out**3) / (kvisc_air * diff_air)
        nuss_up = 0.54 * (Ray_air_up**0.25)
        convection_up_adjust = (cond_air / D_cs_out) * nuss_up

        Ray_air_side = (accel * (1.0/Tf_out) * (Ts_out - temp_0) * ht_cs_out**3) / (kvisc_air * diff_air)
        nuss_side = (0.825 + ( (0.387 * Ray_air_side**(1.0/6.0)) / ((1.0 + ((0.492/prand_air)**(9.0/16.0)) )**(8.0/27.0)) ))**2
        convection_side_adjust = (cond_air / ht_cs_out) * nuss_side

        Ray_air_low = (accel * (1.0/Tf_out) * (Ts_out - temp_0) * D_cs_out**3) / (kvisc_air * diff_air)
        nuss_low = 0.52 * (Ray_air_low**0.2)
        convection_low_adjust = (cond_air / D_cs_out) * nuss_low
    
    CS_temps[0].append(temp_12-273.15)
    CS_temps[1].append(Ts_in-273.15)
    CS_temps[2].append(Ts_out-273.15)
    CS_temps[3].append(temp_0-273.15)
    
    return qloss

def solve_cold_storage(dictionary, CS_temps):
    hd = dictionary['hd']
    mol_air_12 = dictionary['mol_air_12']
    mol_abo3_710 = dictionary['mol_abo3_710']
    mol_abo3_1413 = dictionary['mol_abo3_1413']
    mol_abo3_12 = dictionary['mol_abo3_12']
    cp_abo3 = dictionary['cp_abo3']
    temp_0 = dictionary['temp_0']
    temp_10 = dictionary['temp_10']
    temp_12 = dictionary['temp_12']
    temp_13 = dictionary['temp_13']
    time_step = float(dictionary['time_step'])
    mol_abo3_max = dictionary['mol_abo3_max']
    rho_abo3 = dictionary['rho_abo3']
    dens_abo3 = dictionary['dens_abo3']
    ullage = dictionary['ullage']
    p0air = dictionary['p0air']
    R = dictionary['R']
    CS_insul_thick = dictionary['CS_insul_thick']
    energy_10 = dictionary['energy_10']
    energy_12 = dictionary['energy_12']
    energy_13 = dictionary['energy_13']
    cp_air_in = CP.PropsSI('CPMOLAR', 'T', temp_12, 'P', p0air, 'Air')

    mol_abo3_12last = mol_abo3_12
    mol_abo3_12 = (-mol_abo3_1413 + mol_abo3_710) * time_step + mol_abo3_12last
    volume_abo3_12_avg = (mol_abo3_12last + mol_abo3_12) * 0.5 * rho_abo3 / dens_abo3

    volume_cs = (mol_abo3_max * rho_abo3 / dens_abo3) * (1.0 + ullage)
    r_cs = (volume_cs / (math.pi * 2.0 * hd))**(1.0/3.0)
    ht_cs_in = volume_cs / (math.pi * r_cs**2)

    volume_air_12 = volume_cs - mol_abo3_12 * rho_abo3 / dens_abo3
    mol_air_12last = mol_air_12
    mol_air_12 = (volume_air_12 * p0air) / (R * temp_12)

    if mol_air_12 < mol_air_12last:
        mol_air_24 = 0.0
        mol_air_25 = (mol_air_12last - mol_air_12) / time_step
    elif mol_air_12 == mol_air_12last:
        mol_air_24 = 0.0
        mol_air_25 = 0.0
    elif mol_air_12 > mol_air_12last:
        mol_air_25 = 0.0
        mol_air_24 = (mol_air_12 - mol_air_12last) / time_step

    r_cs_out = r_cs + CS_insul_thick
    D_cs_out = r_cs_out * 2.0
    ht_cs_out = ht_cs_in + 2.0 * CS_insul_thick

    A_circle_h = math.pi * (r_cs**2)
    A_circle_c = math.pi * (r_cs_out**2)
    A_side_c = math.pi * D_cs_out * ht_cs_out

    energy_14 = solve_CS_qloss(dictionary, r_cs, r_cs_out, D_cs_out, ht_cs_in, ht_cs_out, A_circle_h, A_circle_c, A_side_c, CS_temps)
    energy_25 = mol_air_25*cp_air_in*(temp_12-temp_0)
    energy_24 = 0.0

    numerator = energy_12 + energy_10 - energy_13 - energy_25 - energy_14
    denominator = (mol_abo3_12*cp_abo3 + mol_air_12*cp_air_in)/time_step
    temp_12_new = temp_0 + numerator/denominator

    if temp_12_new < temp_0:
        energy_14 = energy_12 + energy_10 - energy_13 - energy_25 - 0.001
        numerator = energy_12 + energy_10 - energy_13 - energy_25 - energy_14
        denominator = (mol_abo3_12*cp_abo3 + mol_air_12*cp_air_in)/time_step
        temp_12_new = temp_0 + numerator/denominator

    if temp_12_new > temp_10:
        if temp_12_new - 1.0 >= temp_10:
            print('Possible error in CS heat loss calculation')
        energy_14_old_high = energy_14
        energy_14 = energy_12 + energy_10 - energy_13 - energy_25 - denominator * (temp_10 - temp_0)
        temp_12_new = temp_10

    energy_12_old = energy_12
    energy_12 = (mol_abo3_12*cp_abo3*(temp_12_new-temp_0) + mol_air_12*cp_air_in*(temp_12_new-temp_0))/time_step

    cserror = energy_12_old + energy_10 - energy_13 - energy_25 - energy_14 - energy_12
    if cserror > 2.0:
        print('CS error is ', cserror, 'W')
    elif cserror < -2.0:
        print('CS error is ', cserror, 'W')

    dictionary['mol_air_12'] = mol_air_12
    dictionary['mol_air_24'] = mol_air_24 
    dictionary['mol_air_25'] = mol_air_25
    dictionary['temp_12'] = temp_12_new
    dictionary['mol_abo3_12'] = mol_abo3_12
    dictionary['energy_12'] = energy_12
    dictionary['energy_14'] = energy_14
    dictionary['energy_24'] = energy_24
    dictionary['energy_25'] = energy_25
    potato = mol_abo3_1413 * cp_abo3 * (temp_13-temp_0)
    dictionary['potato'] = potato
    
    return 0

def solve_work_lift(dictionary):
    mol_abo3_1413 = dictionary['mol_abo3_1413']
    mol_mass_abo3 = dictionary['mol_mass_abo3']
    g_lift = dictionary['g_lift']
    ht_lift = dictionary['ht_lift']
    eta_lift = dictionary['eta_lift']
    parasitics = dictionary['lift_parasitic']
    g_to_kg = 0.001

    energy_15 = (mol_abo3_1413 * mol_mass_abo3 * g_to_kg) * g_lift * ht_lift / eta_lift + parasitics

    dictionary['energy_15'] = energy_15

    return 0

def solve_HX(dictionary):
    mol_o2 = dictionary['mol_o2']
    mol_abo3_1413 = dictionary['mol_abo3_1413']
    cp_abo3 = dictionary['cp_abo3']
    U_hx = dictionary['U_hx']
    area_hx = dictionary['area_hx']
    temp_0 = dictionary['temp_0']
    temp_1 = dictionary['temp_1']
    temp_3 = dictionary['temp_3']
    temp_4 = dictionary['temp_4']
    pSR3 = dictionary['pSR3']

    cp_o2_5 = CP.PropsSI('CPMOLAR', 'T', temp_4, 'P', pSR3, 'Oxygen')
    cp_o2_3 = CP.PropsSI('CPMOLAR', 'T', temp_3, 'P', pSR3, 'Oxygen')
    cp_o2 = 0.5*(cp_o2_5 + cp_o2_3)

    if mol_abo3_1413 > 0.0:
        mol_dot_abo3_12 = dictionary['mol_abo3_12'] / 3600.0
        temp_12 = dictionary['temp_12']
        temp_10 = dictionary['temp_10']
        temp_fraction_CS = min(mol_dot_abo3_12 / mol_abo3_1413, 1.0)
        temp_13 = temp_12 + (1.0 - temp_fraction_CS) * (temp_10 - temp_12)
    else:
    	temp_13 = dictionary['temp_13']

    if mol_abo3_1413 > 0.0:
        c_h = mol_o2 * cp_o2_5
        c_c = mol_abo3_1413 * cp_abo3
        cr = min(c_h,c_c) / max(c_h,c_c)
        NTU = (U_hx * area_hx) / min(c_h,c_c)

        # effectiveness for counterflow HX
        eta_eff = (1.0 - math.exp(-NTU * (1.0-cr))) / (1.0 - cr * math.exp(-NTU * (1.0-cr)))

        # T1 is the temperature of particles exiting the HX
        temp_1 = temp_13 + eta_eff * min(c_h,c_c) * (temp_4-temp_13) / c_c

        # T3 is the temperature of Oxygen exiting the HX
        temp_3 = temp_4 - eta_eff * min(c_h,c_c) * (temp_4-temp_13) / c_h
    else:
        temp_1 = temp_13
        temp_3 = temp_4
        eta_eff = 0.0

    energy_1 = mol_abo3_1413 * cp_abo3 * (temp_1-temp_0)
    energy_3 = mol_o2 * cp_o2_5 * (temp_3-temp_0)
    energy_5 = mol_o2 * cp_o2_5 * (temp_4-temp_0)
    energy_13 = mol_abo3_1413 * cp_abo3 * (temp_13-temp_0)
    energy_19 = energy_5 + energy_13 - energy_1 - energy_3

    if math.fabs(energy_19) > 0.5:
        energy_19 = 0.0

    if math.fabs(energy_5 + energy_13 - energy_1 - energy_3 - energy_19) > 1.0:
        print('HX energy imbalance warning {!s}'.format(energy_5 + energy_13 - energy_1 - energy_3 - energy_19))

    dictionary['temp_1'] = temp_1
    dictionary['temp_3'] = temp_3
    dictionary['energy_1'] = energy_1
    dictionary['energy_3'] = energy_3
    dictionary['energy_5'] = energy_5
    dictionary['energy_13'] = energy_13
    dictionary['energy_19'] = energy_19

    return 0

def calibrate_HX(dictionary):
    mol_o2 = dictionary['mol_o2']
    mol_abo3_1413 = dictionary['mol_abo3_1413']
    cp_abo3 = dictionary['cp_abo3']
    U_hx = dictionary['U_hx']
    area_hx = dictionary['area_hx']
    temp_13 = dictionary['temp_13']
    temp_0 = dictionary['temp_0']
    temp_1 = dictionary['temp_1']
    temp_3 = dictionary['temp_3']
    temp_4 = dictionary['temp_4']
    pSR3 = dictionary['pSR3']

    cp_o2_5 = CP.PropsSI('CPMOLAR', 'T', temp_4, 'P', pSR3, 'Oxygen')

    if mol_abo3_1413 > 0.0:
        c_h = mol_o2 * cp_o2_5
        c_c = mol_abo3_1413 * cp_abo3
        cr = min(c_h,c_c) / max(c_h,c_c)

        # effectiveness for counterflow HX
        eta_eff = 0.85
        
        NTUnum = -math.log((1.0-eta_eff) / (1.0-eta_eff*cr))
        NTUden = 1.0 - cr
        NTU = NTUnum/NTUden
        UA = NTU * min(c_h,c_c)

        # T1 is the temperature of particles exiting the HX
        temp_1 = temp_13 + eta_eff * min(c_h,c_c) * (temp_4-temp_13) / c_c

        # T3 is the temperature of Oxygen exiting the HX
        temp_3 = temp_4 - eta_eff * min(c_h,c_c) * (temp_4-temp_13) / c_h
    else:
        temp_1 = temp_13
        temp_3 = temp_4
        eta_eff = 0.0

    energy_1 = mol_abo3_1413 * cp_abo3 * (temp_1-temp_0)
    energy_3 = mol_o2 * cp_o2_5 * (temp_3-temp_0)
    energy_5 = mol_o2 * cp_o2_5 * (temp_4-temp_0)
    energy_13 = mol_abo3_1413 * cp_abo3 * (temp_13-temp_0)
    energy_19 = energy_5 + energy_13 - energy_1 - energy_3

    if math.fabs(energy_19) < 0.5:
        energy_19 = 0.0

    dictionary['temp_1'] = temp_1
    dictionary['temp_3'] = temp_3
    dictionary['energy_1'] = energy_1
    dictionary['energy_3'] = energy_3
    dictionary['energy_5'] = energy_5
    dictionary['energy_13'] = energy_13
    dictionary['energy_19'] = energy_19

    return UA, eta_eff

def vacuum_pump(dictionary):
    mol_o2 = dictionary['mol_o2']
    R = dictionary['R']
    temp_pump = dictionary['temp_pump']
    p0air = dictionary['p0air']
    p502 = dictionary['p502']
    eta_pump = dictionary['eta_pump']
    energy_16max = dictionary['energy_16max']
    temp_0 = dictionary['temp_0']
    parasitics = dictionary['pump_parasitic']

    energy_16 = mol_o2 * R * temp_0 * math.log(p0air / p502) / eta_pump + parasitics
    energy_16max = max(energy_16, energy_16max)

    dictionary['energy_16'] = energy_16
    dictionary['energy_16max'] = energy_16max

    return 0

def solve_cost(SA_rox,r_sb, A_hs, P_st,mol_abo3_6max,R_rec,D_ap,pi,F_con, F_rec_par,mol_mass_abo3,hd,volume_hs, w0,w1,w2,w3,w4,cost_lay,cost_insfi,cost_pertcon,cost_expbo,cost_reincon,F_ms, C_n2g ,cost_SR3_a, M_SR3,vp_0, vp_1,Per_s,M_ec,M_p,M_SP, C_ppart, F_tpre,F_sca,C_ele_i,hx_0,hx_1,area_hx, F_sv, F_SR3C,CF_uh,F_v,r_rox,L_rox,M_roxcomp,P_r,cost_omi,F_cown,contingency,B_psf,B_ps,F_ti,F_tc,F_tp,F_ts,F_comp,cost_field,cost_fieldprep,A_sf,area_receiver,energy_16max,cp_abo3,temp_7):
    R_pap = 5.702#(D_ap**2)*pi*F_con/4000
    ''' THIS IS 78 for 111.7 MWh system'''
    N_rec= 4 #R_rec/R_pap
    #S_en= mass_abo3_required* rho_E /1000000/h_st
    #A_sf=(R_rec*10^6)/eta_field/DP TODO: make A_sf = to this to start with
    #m_part=(1+F_rec_par)*mass_abo3_needed
    mass_abo3_required = (mol_abo3_6max*mol_mass_abo3)/1000
    #F_ps = 1/(1+F_rec_par)
    #r_sb= math.pow(4*volume_hs/(pi*hd),(.333333))/2
    #A_hs= 4*hd*pi*r_sb**2
    h_gen=hd *2*r_sb
    r0=r_sb + w0
    r1=r0+ w1
    r2=r1+ w2
    r3= r2+ w3
    r4= r3+ w4
    V_0=pi*(r0**2)*h_gen
    V_1=pi*(r1**2)*h_gen
    V_2=pi*(r2**2)*h_gen
    V_3=pi*(r3**2)*h_gen
    V_4=pi*(r4**2)*h_gen
    c0 =(V_0-volume_hs)*cost_lay
    c1=(V_1 - volume_hs)*cost_insfi
    c2=(V_2-V_1)*cost_pertcon
    c3=(V_3-V_2)*cost_expbo
    c4=(V_4-V_3)*cost_reincon
    Cmisc=(c0+c1+c2+c3+c4)*F_ms
    cost_tHB =c0+c1+c2+c3+c4+Cmisc+C_n2g*(mol_abo3_6max/35347180.19)
    #cost_pA=cost_tHB/A_hs
    N_recab= math.ceil(N_rec) #numbers of receivers
    #S_en= mass_abo3_required* rho_E /1000000/h_st
    #1.0 Solve cost of individual components

    cost_SR3= cost_SR3_a*M_SR3*area_receiver*N_recab*(1+Per_s)
    cost_vp=(vp_0+vp_1*(energy_16max/1000)/N_recab)*(1+Per_s+M_ec+M_p)*N_recab #TODO:make sure energy_16 is in kW
    cost_particles=mass_abo3_required*M_SP*C_ppart
    cost_tower=(Per_s+1+M_ec)*F_tpre*R_rec**F_sca
    cost_elevator=R_rec*C_ele_i
    #TODO: fix this and make a way to size the heat exchanger)
    cost_HX=(hx_0+hx_1*area_hx)*(mol_abo3_6max/35347180.19)
    cost_sH=cost_tHB
    cost_sLH= A_hs*F_sv*cost_SR3_a*M_SR3*(1+Per_s)*F_SR3C
    cost_sUH=CF_uh*cost_sLH*F_v

    #cost_Rox_a= ((pi*2*r_rox* L_rox) *cost_SR3_a*M_SR3*M_roxcomp*(1+ Per_s))
    #cost_Rox= cost_Rox_a*.7*mol_abo3_6max/35347180.19
    S1=pi*2*r_rox* L_rox
    S2 =SA_rox
    cost_Rox = SA_rox*cost_SR3_a*M_SR3*M_roxcomp*(1+ Per_s)

    cost_control=(cost_Rox+cost_sUH+cost_sLH+ cost_sH+cost_HX+cost_SR3+ cost_vp+ cost_elevator)*F_comp
    cost_sf=A_sf*(cost_field+cost_fieldprep)
    cost_powerblock =F_ti*F_tc*P_r*1000*F_tp*((P_r*1000)**F_ts)*(1+Per_s+M_ec+M_p)
    cost_bal=P_r*1000*B_ps*P_st**B_psf
    cost_contingency=(cost_SR3+ cost_vp+cost_particles+cost_tower+cost_elevator+cost_powerblock+cost_sf+cost_control+ cost_Rox+cost_sUH+cost_sLH+cost_sH+cost_HX)*contingency
    cost_own =(cost_bal+cost_powerblock+cost_sf+cost_control+ cost_Rox+cost_sUH+cost_sLH+ cost_sH+cost_HX+cost_SR3+ cost_vp+cost_particles+cost_tower+ cost_elevator)* F_cown
    cost_mctot = cost_own+ cost_contingency+cost_bal+cost_powerblock+cost_sf+ cost_control + cost_Rox+cost_sUH+cost_sLH+ cost_sH+cost_HX+cost_SR3+ cost_vp+cost_particles +cost_tower+cost_elevator
   #2.0 DO A BALANCE OF PLANT : Calculate exactly what the name suggests

    cost_total_direct = cost_mctot-cost_own-cost_contingency
    cost_tower_receiver=cost_SR3+cost_tower+cost_elevator+cost_sUH
    cost_thermalstorage=cost_sLH+cost_sH+cost_particles
    cost_solar_control =cost_sf+ cost_control
    cost_gen_rox=cost_Rox+cost_powerblock
    cost_balance_of_Plant_total= cost_bal+ cost_HX+ cost_vp
    #total per year
    cost_ppe=cost_mctot/P_r/1000
    capacity_factor = .531
    estimated_production = P_r*capacity_factor*365*24/1000
    cost_om = estimated_production*(cost_omi*10^6)/(24*365)
    cost_particle_replacement= cost_particles*.1
    capital_cost_year = cost_mctot* .08
    #calculating energy density

    estimated_production = .531*P_r*365*24/1000

    cost_perkWth_om= cost_om/estimated_production/1000000
    cost_perkWth_partrep=cost_particle_replacement/estimated_production/1000000
    cost_perkWth_totalcost=capital_cost_year/estimated_production/1000000
    cost_perkWth_TOTAL =  cost_perkWth_om+cost_perkWth_partrep+cost_perkWth_totalcost

    energy_density = cp_abo3*(1050-(200))/mol_mass_abo3/3600 + .205/2*320/3600/.136 #gives kWh/kg the .000277 is a conversion from kj to kwh #energy stored
    energy_actual = cp_abo3*(1047-(388))/mol_mass_abo3/3600 + .205/2*320/3600/.136
    energy_stored = energy_density*mass_abo3_required #Kwh
    cost_energy = cost_particles*.91/energy_stored

    return cost_SR3,cost_vp,cost_particles,cost_tower,cost_elevator,cost_HX,cost_sH,cost_sLH,cost_sUH,cost_Rox,cost_control,cost_sf,cost_powerblock,cost_bal,cost_contingency,cost_own,cost_mctot,cost_total_direct,cost_tower_receiver,cost_thermalstorage,cost_solar_control,cost_gen_rox,cost_balance_of_Plant_total,cost_ppe,capacity_factor,estimated_production,cost_om,cost_particle_replacement,capital_cost_year, energy_density, energy_stored, cost_energy,cost_perkWth_TOTAL
