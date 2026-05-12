#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 14:28:46 2026

@author: sondrehuse

General description

This file presents the proof calculation of the lateral resistance of a 
a Eurotec angle bracket propesed as a connection between the top beam and CLT wall
for testing at SoFSI in Bristol as part of the ERIES-HYSTERESIS project. 
Verifications have been done according to EN1995-1-1:2025. 

Using:
Full pattern 230x120 angle bracket from Eurotec:
41 Ø5 x 60 Angle-bracket screws + 6 Ø5 x 120 inclined screws (23 degree angle)
Screws are Paneltwistec Countersunk head wood screws
"""

import math

# CLT wall - Stora Enso LS5s-2* C24

wall_t = 160                                      # mm, Thickness of CLT wall, 60-40-60
wall_l = 2000                                     # mm, length of CLT wall
p_k = 350                                         # kg/m^3, C24 density

t_plate = 3                                       # mm, thickness of steel plate
d_screws = 5.0                                    # mm, diameter of the screws
d1_screws = 3.3                                   # mm, inner diameter of screws
d_head_screws = 10                                # mm, diameter of screw head
l_screws_perp = 60                                # mm, length of perpendicular screws
l_screws_inclined = 120                           # mm, length of inclined screws
lw_perp = 36                                      # mm, withdrawal length of perpendicular screws - ref. Figure 11.1 EN1995-1-1:2025
lw_inclined = 70                                  # mm, withdrawal length of inclined screws - ref. Figure 11.1 EN1995-1-1:2025
alfa = 23                                         # degrees, angle of incline 
th1 = t_plate                                     # mm, embedment depth of member 1 (steel plate)
th2_perp = l_screws_perp - t_plate - 2            # mm, embedment depth of member 2 (CLT)
th2_inclined = l_screws_inclined - 9.8            # mm, embedment depth of member 2 (CLT)

#%% Computing the resistance of a single fastener following chapter 11.2 in EN1995-1-1:2025 

# Factors
kmod = 1.1                                                  # Instantaneous load duration, Table 5.4 EN1995-1-1:2025
gamma_r = 1.25                                              # Partial factor for CLT  

k_screw = 8.2                                               # Table 11.2 EN1995-1-1:2025
k_w = 1                                                     # Table 11.2 EN1995-1-1:2025
k_mat = 1                                                   # Table 11.2 EN1995-1-1:2025    
k_p_perp = 0.7                                              # Table 11.2 EN1995-1-1:2025, value for perpendicular
k_p_inclined = 1.1                                          # Table 11.2 EN1995-1-1:2025, value for inclined

f_u_k = 1240                                                # Table 11.3 EN1995-1-1:2025 - Martensitic stainless steel
M_y_k = 5910                                                # Nmm, Value based on ETA-11/0024
k_plt = t_plate / d_screws                                  # Table 11.6 EN1995-1-1:2025 - interpolated value  
f_h1_k = k_plt * 600                                        # N/mm^2, Table 11.6 EN1995-1-1:2025 - for steel plate

f_h2_k_perp = (0.082 * p_k * d_screws**-0.3)                # N/mm^2, Table 11.6 EN1995-1-1:2025 - for perpendicular screws gives 1 in denominator

beta_perp = f_h2_k_perp / f_h1_k                            # Embedment ratio - for perpendicular screws

# Head pull-through resistance 11.2.2.2 EN1995-1-1:2025
F_pull_k = (15 * math.exp(-d_head_screws / 50)              # N
*d_head_screws**2*(p_k/350)**0.8)

# Characteristic witdrawal resistance, Value based on ETA-11/0024
f_w_k_perp = 12.1                                           # N/mm^2

# Characteristic withdrawal resistance 11.2.2.3 EN1995-1-1:2025
F_w_k_perp = math.pi * d_screws * lw_perp * f_w_k_perp      # N, for perpendicular screws

# Steel tensile resistance, 11.2.2.4 EN1995-1-1:2025
F_t_k = 0.9 * math.pi * (d1_screws/2)**2 * f_u_k            # N

# Charactersitc axial resistance
F_ax_k_perp = min(max(F_pull_k, F_w_k_perp), F_t_k)         # N, for perpendicular screws 

# Failure modes per 11.2.3.2(7) - smallest value governs
th_ratio = th2_perp/th1

mode_a_perp = f_h1_k * th1 * d_screws

mode_b_perp = f_h2_k_perp * th2_perp * d_screws

mode_c_perp = (f_h1_k*th1*d_screws/(1+beta_perp)) * (math.sqrt(beta_perp+2*beta_perp**2 * (1 + th_ratio + th_ratio**2) + beta_perp**3*th_ratio**2) - beta_perp *(1+th_ratio))

mode_d_perp = 1.05 * (f_h1_k*th1*d_screws/(2+beta_perp)) * (math.sqrt(2*beta_perp * (1+beta_perp) + (4*beta_perp*(2+beta_perp)*M_y_k) / (f_h1_k * d_screws * th1**2)) - beta_perp)

mode_e_perp = 1.05 * (f_h1_k*th2_perp*d_screws/(1+2*beta_perp)) * (math.sqrt(2*beta_perp**2 * (1+beta_perp) + ((4*beta_perp*(1+2*beta_perp)*M_y_k) / (f_h1_k * d_screws * th2_perp**2))) - beta_perp)

mode_f_perp = 1.15 * math.sqrt(2*beta_perp / (1+beta_perp)) * math.sqrt(2*M_y_k * f_h1_k * d_screws)

mode_check_perp = [mode_a_perp, mode_b_perp, mode_c_perp, mode_d_perp, mode_e_perp, mode_f_perp]
mode_names_perp = ['mode_a', 'mode_b', 'mode_c', 'mode_d', 'mode_e', 'mode_f']

min_index_perp = mode_check_perp.index(min(mode_check_perp))
gov_mode_name_perp = mode_names_perp[min_index_perp]
gov_mode_value_perp = mode_check_perp[min_index_perp]

print("Computing lateral resistance of a single perpendicular screw:")
print(f"Governing failure mode is {gov_mode_name_perp}: {gov_mode_value_perp:.0f} N")

k_rp1 = 0.25     # Table 11.8 EN1995-1-1:2025
k_rp2 = 1        # Table 11.9 EN1995-1-1:2025
    
# 11.2.3.3 (1) - Rope effect (f)
if gov_mode_name_perp == 'mode_f':  
    # Dowel-effect contribution 11.2.3.3(7)
    F_D_k_perp = mode_f_perp   # N
    # Rope-effect contribution - 11.2.3.8 EN1995-1-1:2025
    F_rp_k_perp = min(k_rp1 * F_ax_k_perp, k_rp2 * F_D_k_perp)
elif gov_mode_name_perp == 'mode_e':  
    # Dowel-effect contribution 11.2.3.3(7)
    F_D_k_perp = mode_e_perp   # N
    # Rope-effect contribution - 11.2.3.8 EN1995-1-1:2025
    F_rp_k_perp = min(k_rp1 * F_ax_k_perp, k_rp2 * F_D_k_perp)
elif gov_mode_name_perp == 'mode_d':  
    # Dowel-effect contribution 11.2.3.3(7)
    F_D_k_perp = mode_d_perp   # N
    # Rope-effect contribution - 11.2.3.8 EN1995-1-1:2025
    F_rp_k_perp = min(k_rp1 * F_ax_k_perp, k_rp2 * F_D_k_perp)
else:
    print("Connection is not ductile, failure mode (d), (e) or (f) not governing")
    F_D_k_perp = 0
    F_rp_k_perp = 0
    
# Characteristic lateral resistance per shear plane of single fastener, 11.2.3.1 EN 1995-1-1:2025
F_vk_perp = (F_D_k_perp + F_rp_k_perp) # N
print(f"F_vk_perp = {F_vk_perp:.0f} N")



f_h2_k_inclined =( (0.082 * p_k * d_screws**-0.3) /               # N/mm^2, Table 11.6 EN1995-1-1:2025 - for inclined screws
(2.5 * math.cos(math.radians(90-alfa))**2 
+ math.sin(math.radians(90-alfa))**2) )  

beta_inclined = f_h2_k_inclined / f_h1_k                          # Embedment ratio - for inclined screws

# Characteristic witdrawal strength, Table 11.2 EN1995-1-1:2025
f_w_k_inclined = ( k_screw * k_w * k_mat * d_screws**(-0.33) 
* (p_k/350)**k_p_inclined )                                       # N/mm^2

# Characteristic withdrawal resistance 11.2.2.3 EN1995-1-1:2025
F_w_k_inclined = math.pi * d_screws * lw_inclined * f_w_k_inclined         # N

# Charactersitc axial resistance
F_ax_k_inclined = min(max(F_pull_k, F_w_k_inclined), F_t_k)       # N 

# Failure modes per 11.2.3.2(7)
th_ratio_inclined = th2_inclined/th1

mode_a_inclined = f_h1_k * th1 * d_screws

mode_b_inclined = f_h2_k_inclined * th2_inclined * d_screws

mode_c_inclined = (f_h1_k*th1*d_screws/(1+beta_inclined)) * (math.sqrt(beta_inclined+2*beta_inclined**2 * (1 + th_ratio_inclined + th_ratio_inclined**2) + beta_inclined**3*th_ratio_inclined**2) - beta_inclined *(1+th_ratio_inclined))

mode_d_inclined = 1.05 * (f_h1_k*th1*d_screws/(2+beta_inclined)) * (math.sqrt(2*beta_inclined * (1+beta_inclined) + (4*beta_inclined*(2+beta_inclined)*M_y_k) / (f_h1_k * d_screws * th1**2)) - beta_inclined)

mode_e_inclined = 1.05 * (f_h1_k*th2_inclined*d_screws/(1+2*beta_inclined)) * (math.sqrt(2*beta_inclined**2 * (1+beta_inclined) + ((4*beta_inclined*(1+2*beta_inclined)*M_y_k) / (f_h1_k * d_screws * th2_inclined**2))) - beta_inclined)

mode_f_inclined = 1.15 * math.sqrt(2*beta_inclined / (1+beta_inclined)) * math.sqrt(2*M_y_k * f_h1_k * d_screws)

mode_check_inclined = [mode_a_inclined, mode_b_inclined, mode_c_inclined, mode_d_inclined, mode_e_inclined, mode_f_inclined]
mode_names_inclined = ['mode_a', 'mode_b', 'mode_c', 'mode_d', 'mode_e', 'mode_f']

min_index_inclined = mode_check_inclined.index(min(mode_check_inclined))
gov_mode_name_inclined = mode_names_inclined[min_index_inclined]
gov_mode_value_inclined = mode_check_inclined[min_index_inclined]

print("---------------------------")
print("Computing lateral resistance of a single inclined screw:")
print(f"Governing failure mode is {gov_mode_name_inclined}: {gov_mode_value_inclined:.0f} N")

# 11.2.3.3 (1) - Rope effect (f)
if gov_mode_name_inclined == 'mode_f':  
    # Dowel-effect contribution 11.2.3.3(7)
    F_D_k_inclined = mode_f_inclined   # N
    # Rope-effect contribution - 11.2.3.8 EN1995-1-1:2025
    F_rp_k_inclined = min(k_rp1 * F_ax_k_inclined, k_rp2 * F_D_k_inclined)
elif gov_mode_name_inclined == 'mode_e':  
    # Dowel-effect contribution 11.2.3.3(7)
    F_D_k_inclined = mode_e_inclined   # N
    # Rope-effect contribution - 11.2.3.8 EN1995-1-1:2025
    F_rp_k_inclined = min(k_rp1 * F_ax_k_inclined, k_rp2 * F_D_k_inclined)
elif gov_mode_name_inclined == 'mode_d':  
    # Dowel-effect contribution 11.2.3.3(7)
    F_D_k_inclined = mode_d_inclined   # N
    # Rope-effect contribution - 11.2.3.8 EN1995-1-1:2025
    F_rp_k_inclined = min(k_rp1 * F_ax_k_inclined, k_rp2 * F_D_k_inclined)
else:
    print("Connection is not ductile, failure mode (d), (e) or (f) not governing")
    F_D_k_inclined = 0
    F_rp_k_inclined = 0
    
# Characteristic lateral resistance per shear plane of single fastener, 11.2.3.1 EN 1995-1-1:2025
F_vk_inclined = (F_D_k_inclined + F_rp_k_inclined ) # N
print(f"F_vk_inclined = {F_vk_inclined:.0f} N")


F_vk_ETA = 42       #kN


#%% Total fastener yielding capacity Eurocode 5

n_perp = 41
n_inclined = 4

n_eff_perp = max(n_perp**0.9, 0.9*n_perp)
n_eff_inclined = max(n_inclined**0.9, 0.9*n_inclined)

F_vk_group = n_eff_perp*F_vk_perp + n_eff_inclined*F_vk_inclined





#%% Brittle failure perpendicular to grain - 11.6.1 (4) EN1995-1-1:2025 

kmat_1 = 0.6
kG = 0.05*p_k + 2
bef = 2000
kcon0 = 1
kcon90 = 1
he = 103
h = 2000
Fsp_Rd = kmod / gamma_r * kmat_1 * kG * bef * kcon0 * kcon90 * math.sqrt(he / (1-(he/h)))

print("---------------------------")
print("Splitting capacity:")
print (f"Fsp_Rd = {Fsp_Rd:.0f} N")

#%% Computing the Design resistance for the Angle Bracket

Fv_rk_ETA = 2430                # N, Lateral resistance per fastener, ETA-19/0020

if F_vk_perp > Fv_rk_ETA:    
    Fv_Rk = 47.6                # kN, Value based on screw configuration, ETA-19/0020
    Fv_Rd = kmod * Fv_Rk / gamma_r
    print("---------------------------")
    print(f"Design lateral resistance per Shear Angle Bracket: {Fv_Rd:.0f} kN")

#%% Computing steel and bolt resistance for the Angle Bracket

F_steel_rk = 116   #kN
F_steel_Rd = F_steel_rk/gamma_r


d_bolt = 12

e2 = 1.5*d_bolt
p2 = 3*d_bolt


k1 = min(2.8*e2/d_bolt - 1.7, 1.4*p2/d_bolt - 1.7, 1.5) 
alfa_b = 1
fu = 1000
t = 3
alfa_v = 0.5
f_ub = 1000

A_bolt = math.pi*(d_bolt**2)/4




F_v_bolt_rd = 2 * alfa_v * f_ub * A_bolt / gamma_r


F_b_rd = 2 * k1*alfa_b*f_ub*d_bolt*t/gamma_r

