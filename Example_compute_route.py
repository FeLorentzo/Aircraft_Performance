import numpy as np
import matplotlib.pyplot as plt
import Aircraft_performance_functions as apf

# Nomenclatura

# etapa 1 - subida de 0 a FL200 + descida de FL200 a FL100 + espera de 30 min em máximo endurance + descida de 
# FL100 a 0 + cruzeiro na distância restante

# etapa 2 - subida de 0 a FL410 + descida de FL410 a 0 + cruzeiro no restante

# etapa 3 - recalculando na ordem = etapa 2 -> etapa 1, com consumo contínuo de peso

# Abreviações para plot

# RC = rate os climb
# DT = distance traveled
# FS = fuel spent

# Inputs

S = 92.5 # m²
AR = 8.9
e = 0.85
CD0 = 0.025
Thrust0 = 92300 # N
TSFC = 0.85 / 3600 # Kg/N s
MTOW = 450300
Wf = 130000
W_0fuel = MTOW - Wf
BOW = 320300
K = 1 / (np.pi * AR * e)
g = 9.81
Emax, CL_Emax = apf.compute_Emax(K, CD0)

# Cálculo de subida

# Parâmetros para subida

h_disc_ph1 = 30
W = MTOW
h1 = 0
h2 = 20000
delta_isa_ph1 = 20

v_climb_ph1, gamma_climb_ph1, RC_climb_ph1, time_climb_ph1, DT_climb_ph1, FS_climb_ph1, h_climb_ph1 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph1, delta_isa_ph1, Thrust0, W, S, K, CD0, Emax, CL_Emax, TSFC 
)

# Armazena distância horizontal percorrida e combustível consumido para subida

dist_climb_phase1 = np.cumsum(DT_climb_ph1)[-1]/1000
dfuel_climb_phase1 = np.cumsum(FS_climb_ph1)[-1]
t_climb_phase1 = np.cumsum(time_climb_ph1)[-1]/60

# Visualização subida etapa 1

# apf.show_results(RC_climb_ph1, h_climb_ph1, 'Razão de subida [m/s]', 'Altitude [m]', 'Razão de subida por\nvelocidade calibrada, para alternativa')
# apf.show_results(np.cumsum(time_climb_ph1)/60, np.cumsum(FS_climb_ph1), 'Tempo [min]', 'Combustível consumido [N]', 'Combustível consumido por\ntempo, para alternativa')
# apf.show_results(np.cumsum(DT_climb_ph1)/1000, np.cumsum(FS_climb_ph1), 'Distância horizontal percorrida [km]', 'Combustível consumido [N]', 'Combustível consumido por distância \npercorrida, para alternativa')
# apf.show_results(np.cumsum(DT_climb_ph1)/1000, gamma_climb_ph1, 'Distância horizontal percorrida [km]', 'ângulo de subida [º]', 'Ângulo de subida por \ndistância percorrida, para alternativa')

# Cálculo de descida até 10000 ft

h1 = 20000
h2 = 10000
T_flight_iddle = 0.05 * Thrust0
W = MTOW - dfuel_climb_phase1

v_descent_ph1, gamma_descent_ph1, RC_descent_ph1, time_descent_ph1, DT_descent_ph1, FS_descent_ph1, h_descent_ph1 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph1, delta_isa_ph1, T_flight_iddle, W, S, K, CD0, Emax, CL_Emax, TSFC
)

# Armazena distância horizontal percorrida e combustível consumido para descida

dist_descent_phase1_ft10000 = np.cumsum(DT_descent_ph1)[-1]/1000
dfuel_descent_phase1_ft_10000 = np.cumsum(FS_descent_ph1)[-1]
t_descent_phase1_ft10000 = np.cumsum(time_descent_ph1)[-1]/60

# apf.show_results(RC_descent_ph1, h_descent_ph1, 'Razão de descida [m/s]', 'Altitude [m]', 'Razão de descida por\nvelocidade calibrada, para alternativa')
# apf.show_results(np.cumsum(time_descent_ph1)/60, np.cumsum(FS_descent_ph1), 'Tempo [min]', 'Combustível consumido [N]', 'Combustível consumido por\ntempo, para alternativa')
# apf.show_results(np.cumsum(DT_descent_ph1)/1000, np.cumsum(FS_descent_ph1), 'Distância horizontal percorrida [km]', 'Combustível consumido [N]', 'Combustível consumido por distância \npercorrida, para alternativa')
# apf.show_results(np.cumsum(DT_descent_ph1)/1000, gamma_descent_ph1, 'Distância horizontal percorrida [km]', 'ângulo de descida [º]', 'Ângulo de descida por \ndistância percorrida, para alternativa')

# Cálculo de espera na alternativa

_, _, rho, _ = apf.atmosfera_real(apf.ft2m(10000), delta_isa=20)
W = MTOW - dfuel_climb_phase1 - dfuel_descent_phase1_ft_10000

t_plot, fuel_spent_plot_standby = apf.max_endurance_cruise(rho, S, CD0, K, W, TSFC, t_total=30*60)

# Armazenando combustível consumido e tempo na espera

dfuel_standby_phase1 = np.cumsum(fuel_spent_plot_standby)[-1]
t_standby_phase1 = t_plot[-1]

# apf.show_results(t_plot, np.cumsum(fuel_spent_plot_standby), 'Tempo [min]', 'Combustível consumido [N]', 'Combustível consumido por\ncombustível consumido, para espera na alternativa')

# Restante da descida até pouso

h1 = 10000
h2 = 2000
W = MTOW - dfuel_descent_phase1_ft_10000 - dfuel_climb_phase1 - dfuel_standby_phase1

v_descent2_ph1, gamma_descent2_ph1, RC_descent2_ph1, time_descent2_ph1, DT_descent2_ph1, FS_descent2_ph1, h_descent2_ph1 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph1, delta_isa_ph1, T_flight_iddle, W, S, K, CD0, Emax, CL_Emax, TSFC
)

# Armazena distância horizontal percorrida, combustível consumido e tempo despendido para descida

dist_descent_phase1_ft0 = np.cumsum(DT_descent2_ph1)[-1]/1000
dfuel_descent_phase1_ft_0 = np.cumsum(FS_descent2_ph1)[-1]
t_descent_phase1_ft0 = np.cumsum(time_descent2_ph1)[-1]/60

# apf.show_results(RC_descent2_ph1, h_descent2_ph1, 'Razão de descida [m/s]', 'Altitude [m]', 'Razão de descida por\nvelocidade calibrada, para alternativa')
# apf.show_results(np.cumsum(time_descent2_ph1)/60, np.cumsum(FS_descent2_ph1), 'Tempo [min]', 'Combustível consumido [N]', 'Combustível consumido por\ntempo, para alternativa')
# apf.show_results(np.cumsum(DT_descent2_ph1)/1000, np.cumsum(FS_descent2_ph1), 'Distância horizontal percorrida [km]', 'Combustível consumido [N]', 'Combustível consumido por distância \npercorrida, para alternativa')
# apf.show_results(np.cumsum(DT_descent2_ph1)/1000, gamma_descent2_ph1, 'Distância horizontal percorrida [km]', 'ângulo de descida [º]', 'Ângulo de descida por \ndistância percorrida, para alternativa')

# Cálculo de cruzeiro

h = 20000
W = MTOW - dfuel_climb_phase1 - dfuel_descent_phase1_ft_10000 - dfuel_standby_phase1 - dfuel_descent_phase1_ft_0
delta_x_first_phase = apf.nm2m(200) - dist_climb_phase1 - dist_descent_phase1_ft10000

t_plot_cruise_ph1, x_plot_cruise_ph1, fuel_spent_plot_cruise_ph1 = apf.long_range_cruise_fixed_distance(
    h, delta_isa_ph1, S, K, CD0, TSFC, delta_x_first_phase, W, rho
    )    

# Armazenando combustível consumido, distância percorrida e tempo despedido

dist_cruise_phase1 = x_plot_cruise_ph1[-1]
dfuel_cruise_phase1 = np.cumsum(fuel_spent_plot_cruise_ph1)[-1]
t_cruise_phase1 = t_plot_cruise_ph1[-1]

# apf.show_results(t_plot_cruise_ph1, x_plot_cruise_ph1, 'Tempo [min]', 'Distância percorrida [km]', 'Distância percorrida por\ntempo, para cruzeiro na alternativa')
# apf.show_results(t_plot_cruise_ph1, np.cumsum(fuel_spent_plot_cruise_ph1), 'Tempo [min]', 'Combustível consumido [N]', 'Distância percorrida por\ncombustível consumido, para cruzeiro na alternativa')

######################################################################################################

# Etapa 2 - aeroporto de saída até destino

# Subida

# Parâmetros para subida

h_disc_ph2 = 40
h1 = 0
h2 = 41000
delta_isa_ph2 = 15
W = MTOW - 100*g # Desconta-se 100 Kg de combustível para decolagem

v_climb_ph2, gamma_climb_ph2, RC_climb_ph2, time_climb_ph2, DT_climb_ph2, FS_climb_ph2, h_climb_ph2 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph2, delta_isa_ph2, Thrust0, W, S, K, CD0, Emax, CL_Emax, TSFC
)
    
# Armazena distância horizontal percorrida e combustível consumido para subida

dist_climb_phase2 = np.cumsum(DT_climb_ph2)[-1]/1000
dfuel_climb_phase2 = np.cumsum(FS_climb_ph2)[-1]
t_climb_phase2 = np.cumsum(time_climb_ph2)[-1]/60

# apf.show_results(RC_climb_ph2, h_climb_ph2, 'Razão de subida [m/s]', 'Altitude [m]', 'Razão de subida por\nvelocidade calibrada, para destino')
# apf.show_results(np.cumsum(time_climb_ph2)/60, np.cumsum(FS_climb_ph2), 'Tempo [min]', 'Combustível consumido [N]', 'Combustível consumido por\ntempo, para destino')
# apf.show_results(np.cumsum(DT_climb_ph2)/1000, np.cumsum(FS_climb_ph2), 'Distância horizontal percorrida [km]', 'Combustível consumido [N]', 'Combustível por distância \npercorrida, para destino')
# apf.show_results(np.cumsum(DT_climb_ph2)/1000, gamma_climb_ph2, 'Distância horizontal percorrida [km]', 'ângulo de subida [º]', 'Ângulo de subida por \ndistância percorrida, para destino')

# Descida para etapa 2

h1 = 41000
h2 = 0
W = MTOW - 100*g - dfuel_climb_phase2
T_flight_iddle = 0.05 * Thrust0

v_descent_ph2, gamma_descent_ph2, RC_descent_ph2, time_descent_ph2, DT_descent_ph2, FS_descent_ph2, h_descent_ph2 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph2, delta_isa_ph2, T_flight_iddle, W, S, K, CD0, Emax, CL_Emax, TSFC
)

# Armazena distância horizontal percorrida e combustível consumido para descida

dist_descent_phase2 = np.cumsum(DT_descent_ph2)[-1]/1000
dfuel_descent_phase2 = np.cumsum(FS_descent_ph2)[-1]
t_descent_phase2 = np.cumsum(time_descent_ph2)[-1]/60

# apf.show_results(RC_descent_ph2, h_descent_ph2, 'Razão de descida [m/s]', 'Altitude [m]', 'Razão de descida por\nvelocidade calibrada, para destino')
# apf.show_results(np.cumsum(time_descent_ph2)/60, np.cumsum(FS_descent_ph2), 'Tempo [min]', 'Combustível consumido [N]', 'Combustível consumido por\ntempo, para destino')
# apf.show_results(np.cumsum(DT_descent_ph2)/1000, np.cumsum(FS_descent_ph2), 'Distância horizontal percorrida [km]', 'Combustível consumido [N]', 'Combustível consumido por distância \npercorrida, para destino')
# apf.show_results(np.cumsum(DT_descent_ph2)/1000, gamma_descent_ph2, 'Distância horizontal percorrida [km]', 'ângulo de descida [º]', 'Ângulo de descida por \ndistância percorrida, para destino')

# Cruzeiro para etapa 2

h = 41000
W = (MTOW - dfuel_climb_phase1 - dfuel_descent_phase1_ft_10000 - dfuel_standby_phase1 
     - dfuel_descent_phase1_ft_0 - dfuel_cruise_phase1 - dfuel_climb_phase2 - 
     dfuel_descent_phase2)

leftover_fuel = 0.1 # Porcentagem de combustível de sobra desejado

t_plot_cruise_ph2, x_plot_cruise_ph2, fuel_spent_plot_cruise_ph2 = apf.long_range_cruise_fixed_W(
    h, delta_isa_ph2, S, K, CD0, TSFC, W, W_0fuel, leftover_fuel, Wf, rho
    )    

# Armazenando combustível consumido

dist_cruise_phase2 = x_plot_cruise_ph2[-1]
dfuel_cruise_phase2 = np.cumsum(fuel_spent_plot_cruise_ph2)[-1]
t_cruise_phase2 = t_plot[-1]

# apf.show_results(t_plot_cruise_ph2, x_plot_cruise_ph2, 'Tempo [min]', 'Distância percorrida [km]', 'Distância percorrida por\ntempo, para cruzeiro na alternativa')
# apf.show_results(t_plot_cruise_ph2, np.cumsum(fuel_spent_plot_cruise_ph2), 'Tempo [min]', 'Combustível consumido [N]', 'Distância percorrida por\ncombustível consumido, para cruzeiro na alternativa')

# Resultados totais após estimetiva

dfuel_estimated = (dfuel_climb_phase2 + dfuel_cruise_phase2 + dfuel_descent_phase2 + 
                   dfuel_climb_phase1 + dfuel_cruise_phase1 + dfuel_descent_phase1_ft_10000 + 
                   dfuel_standby_phase1 + dfuel_descent_phase1_ft_0)

dist_estimated = (dist_climb_phase2 + dist_cruise_phase2 + dist_descent_phase2 + 
                  dist_climb_phase1 + dist_cruise_phase1 + dist_descent_phase1_ft10000 + 
                  dist_descent_phase1_ft0)

t_estimated = (t_climb_phase2 + t_cruise_phase2 + t_descent_phase2 + 
               t_climb_phase1 + t_cruise_phase1 + t_descent_phase1_ft10000 + 
               t_standby_phase1 + t_descent_phase1_ft0)

print(f'Combustível restante: {Wf-dfuel_estimated} [N]')
print(f'Distância percorrida total, somando alternativa: {dist_estimated} [km]')
print(f'Tempo total do percurso, incluindo espera na alternativa: {t_estimated} [min]')

####################################################################################################

# Etapa 3 - cálculo acoplando tudo

# Partida até chegada

h_disc_ph2 = 40
h1 = 0
h2 = 41000
delta_isa_ph2 = 15

# Subida

W = MTOW - 100*g # Desconta-se 100 Kg de combustível para decolagem

_, _, _, time_climb_ph2, DT_climb_ph2, FS_climb_ph2, h_climb_ph2 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph2, delta_isa_ph2, Thrust0, W, S, K, CD0, Emax, CL_Emax, TSFC
)

# Cruzeiro

h = 41000
h_real, _, _, _ = apf.atmosfera_real(apf.ft2m(h), delta_isa_ph2)
W -= np.cumsum(FS_climb_ph2)[-1]
dfuel_cruise_estimated = 52300 # [N]

leftover_fuel = (Wf - np.cumsum(FS_climb_ph2)[-1] - dfuel_cruise_estimated ) / Wf 
# Porcentagem de combustível de sobra desejado

time_cruise_ph2, DT_cruise_ph2, FS_cruise_ph2 = apf.long_range_cruise_fixed_W(
    h, delta_isa_ph2, S, K, CD0, TSFC, W, W_0fuel, leftover_fuel, Wf, rho
    )
h_cruise_ph2 = h_real * np.ones(len(time_cruise_ph2))

# Descida

h1 = 41000
h2 = 0
W -=  np.cumsum(fuel_spent_plot_cruise_ph2)[-1]
T_flight_iddle = 0.05 * Thrust0

_, _, _, time_descent_ph2, DT_descent_ph2, FS_descent_ph2, h_descent_ph2 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph2, delta_isa_ph2, T_flight_iddle, W, S, K, CD0, Emax, CL_Emax, TSFC
)

# Subida em direção à alternativa

delta_isa_ph1 = 20
h_disc_ph1 = 30
h1 = 0
h2 = 20000

W -= np.cumsum(FS_descent_ph2)[-1]


_, _, _, time_climb_ph1, DT_climb_ph1, FS_climb_ph1, h_climb_ph1 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph1, delta_isa_ph1, Thrust0, W, S, K, CD0, Emax, CL_Emax, TSFC 
)

# Cruzeiro em direção à alternativa

h = 20000
h_real, _, _, _ = apf.atmosfera_real(apf.ft2m(h), delta_isa_ph1)
W -= np.cumsum(FS_climb_ph1)[-1]
dfuel_cruise_estimated = 7300 # [N]

leftover_fuel = 0.16
# Sobra de combustível a partir da estimativa de combustível consumido

time_cruise_ph1, DT_cruise_ph1, FS_cruise_ph1 = apf.long_range_cruise_fixed_W(
    h, delta_isa_ph1, S, K, CD0, TSFC, W, W_0fuel, leftover_fuel, Wf, rho
    )
h_cruise_ph1 = h_real * np.ones(len(time_cruise_ph1))

# Descida até 1000 ft

h1 = 20000
h2 = 10000
W -= np.cumsum(fuel_spent_plot_cruise_ph1)[-1]

_, _, _, time_descent_ph1, DT_descent_ph1, FS_descent_ph1, h_descent_ph1 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph1, delta_isa_ph1, T_flight_iddle, W, S, K, CD0, Emax, CL_Emax, TSFC
)

# Espera por 30 min

h = apf.ft2m(10000)
_, _, rho, _ = apf.atmosfera_real(apf.ft2m(h), delta_isa_ph1)
W -= np.cumsum(FS_descent_ph1)[-1]

time_standby, FS_standby = apf.max_endurance_cruise(rho, S, CD0, K, W, TSFC, t_total=30*60)
DT_standby = np.cumsum(DT_descent_ph1)[-1]/1000 * np.ones(len(time_standby))
h_standby = h * np.ones(len(time_standby))

# Descida até alternativa

h1 = 10000
h2 = 2000
W -= np.cumsum(fuel_spent_plot_standby)[-1]

_, _, _, time_descent2_ph1, DT_descent2_ph1, FS_descent2_ph1, h_descent2_ph1 = apf.compute_climbing_vectors(
    h1, h2, h_disc_ph1, delta_isa_ph1, T_flight_iddle, W, S, K, CD0, Emax, CL_Emax, TSFC
)

# Compilando resultados 

distance = np.append(np.cumsum(DT_climb_ph2)/1000, np.cumsum(DT_climb_ph2)[-1]/1000 + DT_cruise_ph2)
distance = np.append(distance, distance[-1] + np.cumsum(DT_descent_ph2)/1000)
distance = np.append(distance, distance[-1] + np.cumsum(DT_climb_ph1)/1000)
distance = np.append(distance, distance[-1] + DT_cruise_ph1)
distance = np.append(distance, distance[-1] + np.cumsum(DT_descent_ph1)/1000)
distance = np.append(distance, distance[-1] + DT_standby)
distance = np.append(distance, distance[-1] + np.cumsum(DT_descent2_ph1)/1000)

time = np.append(np.cumsum(time_climb_ph2)/60, np.cumsum(time_climb_ph2)[-1]/60 + time_cruise_ph2)
time = np.append(time, time[-1] + np.cumsum(time_descent_ph2)/60)
time = np.append(time, time[-1] + np.cumsum(time_climb_ph1)/60)
time = np.append(time, time[-1] + time_cruise_ph1)
time = np.append(time, time[-1] + np.cumsum(time_descent_ph1)/60)
time = np.append(time, time[-1] + time_standby)
time = np.append(time, time[-1] + np.cumsum(time_descent2_ph1)/60)


dfuel = np.append(np.cumsum(FS_climb_ph2), np.cumsum(FS_climb_ph2)[-1] + np.cumsum(FS_cruise_ph2))
dfuel = np.append(dfuel, dfuel[-1] + np.cumsum(FS_descent_ph2))
dfuel = np.append(dfuel, dfuel[-1] + np.cumsum(FS_climb_ph1))
dfuel = np.append(dfuel, dfuel[-1] + np.cumsum(FS_cruise_ph1))
dfuel = np.append(dfuel, dfuel[-1] + np.cumsum(FS_descent_ph1))
dfuel = np.append(dfuel, dfuel[-1] + np.cumsum(FS_standby))
dfuel = np.append(dfuel, dfuel[-1] + np.cumsum(FS_descent2_ph1))

altitude = np.append(h_climb_ph2, h_cruise_ph2)
altitude = np.append(altitude, h_descent_ph2)
altitude = np.append(altitude, h_climb_ph1)
altitude = np.append(altitude, h_cruise_ph1)
altitude = np.append(altitude, h_descent_ph1)
altitude = np.append(altitude, h_standby)
altitude = np.append(altitude, h_descent2_ph1)


apf.show_results(time, altitude, 'Tempo [min]', 'Altitude [km]', 'Evolução temporal da altitude')
apf.show_results(time, distance, 'Tempo [min]', 'Distância horizontal percorrida [km]', 'Evolução temporal da \ndistância percorrida')
apf.show_results(distance, altitude, 'Distância [km]', 'Altitude [km]', 'Trajetória da aeronave')
apf.show_results(distance, MTOW - dfuel, 'Distância [Km]', 'Combustível [N]', 'Consumo de combustível\nao lngo da distância')

plt.show()
