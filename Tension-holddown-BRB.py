#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 12:09:16 2026

@author: sondrehuse
"""
import math

wall_t = 200                                      # mm, Thickness of CLT wall, 60-40-60

wall_l = 2000                                     # mm, length of CLT wall

p_k = 350                                        # kg/m^3, C24 density

 

t_plate = 3                                       # mm, thickness of steel plate

d_screws = 5.0                                    # mm, diameter of the screws

d1_screws = 3.2                                   # mm, inner diameter of screws

d_head_screws = 7.2                           # mm, diameter of screw head

l_screws = 60                                # mm, length of perpendicular screws


lw = 51                           # mm, withdrawal length of perpendicular screws - ref. Figure 11.1 EN1995-1-1:2025


th1 = t_plate                                     # mm, embedment depth of member 1 (steel plate)

th2 = l_screws - t_plate - 3            # mm, embedment depth of member 2 (CLT)


 

#%% Computing the resistance of a single fastener following chapter 11.2 in EN1995-1-1:2025


# Factors

kmod = 1.1                                                  # Instantaneous load duration, Table 5.4 EN1995-1-1:2025

gamma_r = 1.25                                              # Partial factor for CLT 

 

k_screw = 8.2                                               # Table 11.2 EN1995-1-1:2025

k_w = 1                                                     # Table 11.2 EN1995-1-1:2025

k_mat = 1                                                   # Table 11.2 EN1995-1-1:2025   

k_p = 1.1                                              # Table 11.2 EN1995-1-1:2025, value for perpendicular

k_p_inclined = 0.7                                         # Table 11.2 EN1995-1-1:2025, value for inclined

 
As = math.pi*(d1_screws/2)**2

F_t_k = 7900                                                # N from ETA-11/0024 of 2025/08/20

f_u_k =  F_t_k/As                                           # Based on F_t_k from ETA f_u_k = F_t_k/As

M_y_k = 4300          # ETA 11/0024
k_plt = t_plate / d_screws                                  # Table 11.6 EN1995-1-1:2025 - interpolated value 

f_h1_k = k_plt * 600                                        # N/mm^2, Table 11.6 EN1995-1-1:2025 - for steel plate

 

f_h2_k = (0.082 * p_k * d_screws**-0.3)/(2.5*math.cos(math.radians(90))**2 + math.sin(math.radians(90))**2) # N/mm^2, Table 11.6 EN1995-1-1:2025 - for perpendicular screws gives 1 in denominator

 
beta = f_h2_k / f_h1_k                            # Embedment ratio - for perpendicular screws

 

# Head pull-through resistance 11.2.2.2 EN1995-1-1:2025

F_pull_k = (15 * math.exp(-d_head_screws / 50)              # N
            *d_head_screws**2*(p_k/350)**0.8)

 

# Characteristic witdrawal strength, Table 11.2 EN1995-1-1:2025

f_w_k = ( k_screw * k_w * k_mat * d_screws**(-0.33)

* (p_k/350)**k_p )                                     # N/mm^2

 

# Characteristic withdrawal resistance 11.2.2.3 EN1995-1-1:2025

F_w_k = math.pi * d_screws * lw * f_w_k      # N, for perpendicular screws

 

# Charactersitc axial resistance

F_ax_k = min(max(F_pull_k, F_w_k), F_t_k)         # N, for perpendicular screws

 

# Failure modes per - smallest value governs

th_ratio = th2/th1

 

mode_a = f_h1_k * th1 * d_screws

 

mode_b = f_h2_k * th2 * d_screws

 

mode_c = (f_h1_k*th1*d_screws/(1+beta)) * (math.sqrt(beta+2*beta**2 * (1 + th_ratio + th_ratio**2) + beta**3*th_ratio**2) - beta *(1+th_ratio))

 

mode_d = 1.05 * (f_h1_k*th1*d_screws/(2+beta)) * (math.sqrt(2*beta * (1+beta) + (4*beta*(2+beta)*M_y_k) / (f_h1_k * d_screws * th1**2)) - beta)

 

mode_e = 1.05 * (f_h1_k*th2*d_screws/(1+2*beta)) * (math.sqrt(2*beta**2 * (1+beta) + ((4*beta*(1+2*beta)*M_y_k) / (f_h1_k * d_screws * th2**2))) - beta)

 

mode_f = 1.15 * math.sqrt(2*beta / (1+beta)) * math.sqrt(2*M_y_k * f_h1_k * d_screws)

 

mode_check = [mode_a, mode_b, mode_c, mode_d, mode_e, mode_f]

mode_names = ['mode_a', 'mode_b', 'mode_c', 'mode_d', 'mode_e', 'mode_f']

 

min_index = mode_check.index(min(mode_check))

gov_mode_name = mode_names[min_index]

gov_mode_value = mode_check[min_index]

 

print("Computing lateral resistance of a single perpendicular screw:")

print(f"governing failure mode is {gov_mode_name}: {gov_mode_value:.0f} N")

 

k_rp1 = 0.25     # Table 11.8 EN1995-1-1:2025

k_rp2 = 1        # Table 11.9 EN1995-1-1:2025

   

# 11.2.3.3 (1) - Rope effect (f)

if gov_mode_name == 'mode_f': 

    # Dowel-effect contribution 11.2.3.3(7)

    F_D_k = mode_f   # N

    # Rope-effect contribution - 11.2.3.8 EN1995-1-1:2025

    F_rp_k = min(k_rp1 * F_ax_k, k_rp2 * F_D_k)

elif gov_mode_name == 'mode_e': 

    # Dowel-effect contribution 11.2.3.3(7)

    F_D_k = mode_e   # N

    # Rope-effect contribution - 11.2.3.8 EN1995-1-1:2025

    F_rp_k = min(k_rp1 * F_ax_k, k_rp2 * F_D_k)

elif gov_mode_name == 'mode_d': 

    # Dowel-effect contribution 11.2.3.3(7)

    F_D_k = mode_d   # N

    # Rope-effect contribution - 11.2.3.8 EN1995-1-1:2025

    F_rp_k = min(k_rp1 * F_ax_k, k_rp2 * F_D_k)

else:

    print("Connection is not ductile, failure mode (d), (e) or (f) not governing")

    F_D_k = 0

    F_rp_k = 0

   

# Characteristic lateral resistance per shear plane of single fastener, 11.2.3.1 EN 1995-1-1:2025

F_vk = (F_D_k + F_rp_k) # N

print(f"F_vk = {F_vk:.0f} N")

#Resistance of a whole fastner
n = 25
n_eff = max(n**0.9, 0.9*n)                                   # Effective number of fastners for WBS screws, from ETA-23/0813
n90 = 3 

F_vd = kmod*n_eff*F_vk/gamma_r                      # EC5 11.3.3(3) Splitting resistance

# Resistance of steel plate

F_1_steel_Rk = 57100                                #N  Characteristic load-carrying capacity of the hold-down, from ETA-23/08
F_1_steel_Rd = F_1_steel_Rk / gamma_r # N 



# Row shear resistance
F_vd_t = kmod*3.9/gamma_r

k_br_v = 0.75  #EC5 Table 11.23

alfa_cl = 0.7*(t_plate/d_screws)+0.3                #EC5 11.5.10(5)

t_ef = min(th2, alfa_cl*7*th2/(3+th2/d_screws))

a1 = 17.5

n0 = 11
a_3t = 150

l_con = a1 * (n0 - 1) + a_3t                                     #mm

Fv_la_d = k_br_v*t_ef*l_con*F_vd_t                    #EC5 11.5.9.1(1) Side shear plane resistance
F_rs_Rd = 2*n90*Fv_la_d                               #EC5 11.5.5(1) Row shear resistance

# Plug shear resistance
k_br_t = 1.25

b_net = 40 - 2*d_screws                                #mm
f_t_0d = kmod*14.5/gamma_r
b_con = b_net

F_t_d = k_br_t*b_net*t_ef*f_t_0d
F_v_bot_d = k_br_v*l_con*b_con*F_vd_t
F_ps_Rd = max(2*Fv_la_d, F_t_d + F_v_bot_d)         #EC 11.5.7 Plug shear resistance

# Tension failure of the threaded steel rod
d_bolt = 16
A_s = math.pi*(d_bolt/2)**2
f_u_b = 1000 
R_t_d = 0.9*f_u_b*A_s/gamma_r


print("\n--- GLOBAL CONNECTION CAPACITY CHECK ---")

capacities = {
    "Steel rod tension R_t_d": R_t_d,
    "Plug shear resistance F_ps_Rd": F_ps_Rd,
    "Row shear resistance F_rs_Rd": F_rs_Rd,
    "Steel plate capacity F_1_steel_Rd": F_1_steel_Rd,
    "Splitting capacity F_vd": F_vd
}


gov_capacity_name = min(capacities, key=capacities.get)
gov_capacity_value = capacities[gov_capacity_name]


for name, value in capacities.items():
    print(f"{name:<40s} = {value:10.0f} N")

print("\n----------------------------------------")
print(f"GOVERNING CAPACITY: {gov_capacity_name}")
print(f"Design resistance R_d = {gov_capacity_value:.0f} N")
print("----------------------------------------")


"""
Plug shear resistance in CLT
"""
"""
t_ef: Two different values. When t_ef_el is governing - Pure brittle failure. When t_ef_pl is governing - mixed failure
"""

#t_ef = 2*math.sqrt(M_y_k/(f_h2_k*d1_screws))
t_ef = 42.69

t1 = 40 
t2 = 17

wc = 45
wm = 100
t_par = 40
t_per = 40
L = 340
l2 = 177.5

E0 = 12000
Gr = 50 #N/mm^2
G = 690 #N/mm^2

ft = 14.5 
fv = 3.9
fv_r = 0.33 * fv


Ath = wc * t1


kh = (2*E0*Ath)/l2

kr = Gr*L*(wc + wm)/(2*t_per)

ka = 4*G*t_per*L/wc

# Because the shear deformation is assumed uniform across its thickness, the stiffness becomes:
    
ka_c = 2 * ka * (t_ef - t_par)/t_per
kr_c = kr*t_per/(t_par + t_per - t_ef)

# Joint resistance governed by the head tensile plane failure

pw_h = ft*wc*t_par*(kh + kr_c + ka_c)/kh


# Joint resistance governed by the bottom shear plane failure
kd = kr_c

pw_d = fv_r*wc*L*(kh + kr_c + ka_c)/kd


# Joint resistance governed by the adjacent shear planes failure

A_va = 2 * (wm - wc)*L

pw_a = fv_r * A_va * (kh + kr_c + ka_c)/ka_c


n_p = 1

pw_c = min(pw_h, pw_d, pw_a)

pw_c_d = kmod*pw_c/gamma_r


#P_w_C governs and the joint resistance governed by the head plane and is the first layer to fail --> k_h = 0

Kd2 = kr
Ka2 = ka
SUM_K_D2 = kr + ka
SUM_K_A2 = kr + ka
A_vb = wc*L

pw_d2 = fv_r*A_vb*SUM_K_D2/Kd2

pw_a2 = fv_r*A_va*SUM_K_A2/Ka2

pw_extra_capacity = min(pw_a2, pw_d2)

if pw_extra_capacity > pw_c:
    pw_final = pw_extra_capacity
    print(pw_final)
    
else:
    pw_final = pw_c
    print(pw_final)
    
        
# Pw needs to be less than Pr with reduced penetration length + rolling shear resistance


L_p = 57 #mm

L_p_red = 57 - t_par

th1_red = L_p_red

th2_red = 40


f_h2k_red = f_h2_k
f_h1k_red = f_h2_k

beta_red = f_h2k_red/f_h1k_red
th_ratio_red = th2_red/th1_red


mode_a_red = f_h1k_red * th1_red * d_screws

 

mode_b_red = f_h2k_red * th2_red * d_screws

 

mode_c_red = (f_h1k_red*th1_red*d_screws/(1+beta_red)) * (math.sqrt(beta_red+2*beta_red**2 * (1 + th_ratio_red + th_ratio_red**2) + beta_red**3*th_ratio_red**2) - beta_red *(1+th_ratio_red))

 

mode_d_red = 1.05 * (f_h1k_red*th1_red*d_screws/(2+beta_red)) * (math.sqrt(2*beta_red * (1+beta_red) + (4*beta_red*(2+beta_red)*M_y_k) / (f_h1k_red * d_screws * th1_red**2)) - beta_red)

 

mode_e_red = 1.05 * (f_h1k_red*th2_red*d_screws/(1+2*beta_red)) * (math.sqrt(2*beta_red**2 * (1+beta_red) + ((4*beta_red*(1+2*beta_red)*M_y_k) / (f_h1k_red * d_screws * th2_red**2))) - beta_red)

 

mode_f_red = 1.15 * math.sqrt(2*beta_red / (1+beta_red)) * math.sqrt(2*M_y_k * f_h1k_red * d_screws)


gov_red = min(mode_a_red, mode_b_red, mode_c_red, mode_d_red, mode_e_red, mode_f_red)*n_eff


limit_red = gov_red + fv_r*A_vb















