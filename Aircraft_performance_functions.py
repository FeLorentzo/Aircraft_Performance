import numpy as np
import matplotlib.pyplot as plt

# Funções de conversão de unidades

def ft2m(ft):
    return 0.3048 * ft


def ftmin2ms(ftmin):
    return 0.00508 * ftmin


def nm2m(nm):
    return 1852 * nm


# Funções para atmosfera


def atmosisa(h):
    '''
    Return [T, P, rho, a]

    h = meters

    T = kelvin
    P = Pa
    rho = Kg/m³
    a = m/s
    '''
    T0 = 288.15
    P0 = 101325
    g = 9.81
    asl = -0.0065
    R = 287.1
    gamma = 1.4

    if h < 11000:
        T = T0 + asl * h
        P = P0 * (T/T0)**(-g/R/asl)
        rho = P / R / T
    else:
        T = 216.65
        P11 = P0 * (T/T0)**(-g/R/asl)
        rho11 = P11 / R / T
        P = P11 * np.exp(-g/R/T * (h-11000))
        rho = rho11 * np.exp(-g/R/T * (h-11000))

    a = np.sqrt(gamma * R * T)

    return T, P, rho, a


def atmosfera_real(h_isa, delta_isa):
    T0 = 288.15
    rho0 = 1.225
    g = 9.81
    asl = -0.0065
    R = 287.1

    
    Temp,_,_,_ = atmosisa(h_isa)
    T_isa = Temp
    delta = (T_isa/T0)**(-g/R/asl)
    T = T_isa+delta_isa
    theta = (T/T0)
    sigma = delta/theta
    rho = rho0*sigma
    h = T/T_isa * h_isa

    return h, T, rho, sigma


# Funções para subida e descida

def correct_thrust(h, thrust0, sigma):
    if h < ft2m(36890):
        thrust = thrust0 * sigma**(0.7)
    else:
        thrust = thrust0 * 1.439 * sigma
    return thrust

def compute_Emax(K, CD0):
    CL_Emax = np.sqrt(CD0/K)
    CD_Emax = 2*CD0
    return CL_Emax/CD_Emax, CL_Emax


def compute_ascending_flight(h, thrust0, sigma, Emax, S, phi, rho, CL_Emax, delta_h, TSFC, W):
    
    thrust = correct_thrust(ft2m(h), thrust0, sigma)
    
    step_climb = ftmin2ms(500)

    # Cálculo da velocidade

    v_cas = np.sqrt(2*W/S/rho/CL_Emax)

    # Cálculo do ângulo de máxima razão de subida
    gamma_max = thrust/W - 1/Emax

    # Cálculo da máxima razão de subida
    rate_of_climb_max = (v_cas * gamma_max)/(1+phi)

    # Verificação da necessidade de step climb
    if rate_of_climb_max < step_climb and rate_of_climb_max > 0:
        print('Precisa de step climb')

    # Cálculo do tempo de subida
    time_to_climb = delta_h / rate_of_climb_max

    # Cálculo da velocidade horizontal
    v_horizontal = v_cas * np.cos(gamma_max)

    # Cálculo da distância horizontal percorrida
    distance_travelled = v_horizontal * time_to_climb

    # Cálculo do consumo de combustível
    fuel_spent = thrust * TSFC * time_to_climb 

    return v_cas, gamma_max, rate_of_climb_max, time_to_climb, distance_travelled, fuel_spent

def get_acc_factor(T, S, W, rho, CD0, K):

    # Definindo parâmetros termodinâmicos

    gamma=1.4
    R = 287.14

    # Pré-alocando o vetor phi

    phi = []
    M = []

    # Definindo parâmetros aerodinâmicos

    CL_Emax = np.sqrt(CD0/K)
    v_cas_acc = []
    a = []


    # Cálculo do Mach, usando relação termodinâmica com temperatura real e
    # velocidade calibrada da aeronave

    a = np.sqrt(gamma*R*T)
    v_cas_acc = np.sqrt(2*W/S/rho/CL_Emax)
    M = v_cas_acc/a

    # Cálculo do fator de aceleração

    phi = ( ( (1+0.2*M**2)**3.5 - 1)/(1 + 0.2*M**2)**2.5 )

    return phi


def compute_climbing_vectors(h1, h2, h_disc, delta_isa, Thrust, W, S, K, CD0, Emax, CL_Emax, TSFC):
    
    # Cálculo de parâmetros atmosféricos
    
    h_isa = np.linspace(ft2m(h1), ft2m(h2), h_disc)
    delta_h_isa = h_isa[1]-h_isa[0]
    
    # Pré-alocando vetores

    v_plot = []
    rate_of_climb_plot = []
    time_to_climb_plot = []
    distance_travelled_plot = []
    fuel_spent_plot = []
    h_plot = []
    gamma_plot = []


    for h in h_isa:
        h_real, T, rho, sigma = atmosfera_real(h, delta_isa=delta_isa)
        h_real_next,_,_,_ = atmosfera_real(h+delta_h_isa, delta_isa=delta_isa)
        delta_h = h_real_next-h_real
        
        phi = get_acc_factor(T, S, W, rho, CD0, K)
        
        v, gamma_max, rate_of_climb_max, time_to_climb, distance_travelled, fuel_spent = compute_ascending_flight(
        h=h_real, thrust0=Thrust, sigma=sigma, Emax=Emax, S=S, phi=phi, rho=rho, CL_Emax=CL_Emax, 
        delta_h=delta_h, TSFC=TSFC, W=W)
        
        W -= fuel_spent 
        
        v_plot.append(v)
        gamma_plot.append(gamma_max*180/np.pi)
        rate_of_climb_plot.append(rate_of_climb_max)
        time_to_climb_plot.append(time_to_climb)
        distance_travelled_plot.append(distance_travelled)
        fuel_spent_plot.append(fuel_spent)
        h_plot.append(h_real)
    
    return v_plot, gamma_plot, rate_of_climb_plot, time_to_climb_plot, distance_travelled_plot, fuel_spent_plot, h_plot


# Funções para cruzeiro

def long_range_cruise_fixed_distance(h, delta_isa, S, K, CD0, TSFC, goal_distance, W0, rho):
    
    _, _, rho, _ = atmosfera_real(ft2m(h), delta_isa)
    
    # Pré alocando vetores

    t_plot = []
    x_plot = []
    fuel_spent_plot = []

    t_plot.append(0)
    x_plot.append(0)
    fuel_spent_plot.append(0)
    
    # Definindo variáveis

    V = np.linspace(10, 500, 100)

    dt = 1; #[s]
    t = 0; #[s]
    x = 0; #[m]

    W = W0

    while x < goal_distance: # Condição para distância a ser percorrida
    
        # Hipótese de pequenos ângulos para alpha
        L = W

        # Range de coeficientes aerodinâmicos
        CL = (2*L)/(rho*S*V**2)
        CD = CD0 + K*CL**2

        # Empuxo necessário para cruzeiro
        Thrust = 0.5*rho*V**2*S*CD

        # Combustível consumido
        delta_fuel = TSFC*Thrust*dt #[N]

        # Alcance específico
        SR = V/(delta_fuel/dt)

        # Velocidade de LRC - condição onde alcance específico é igual a 99% do
        # máximo alcance específico

        LRC, idx = get_long_range_cruise_velocity(V,SR)
        
        t += dt
        x += LRC*dt
        W -= delta_fuel[idx]
        
        t_plot.append(t/60) # min
        x_plot.append(x/1000) # km
        fuel_spent_plot.append(delta_fuel[idx])

    return t_plot, x_plot, fuel_spent_plot


def long_range_cruise_fixed_W(h, delta_isa, S, K, CD0, TSFC, W, W_0fuel, leftover_fuel, Wf, rho):
    
    _, _, rho, _ = atmosfera_real(ft2m(h), delta_isa)
    
    # Definindo variáveis

    V = np.linspace(10, 500, 100)

    dt = 1; #[s]
    t = 0; #[s]
    x = 0; #[m]

    t_plot = []
    x_plot = []
    delta_fuel_plot = []

    t_plot.append(0)
    x_plot.append(0)
    delta_fuel_plot.append(0)

    W = W

    while W > (W_0fuel+leftover_fuel*Wf): # Condição para reserva de combustível

        # Hipótese de pequenos ângulos para alpha
        L = W

        # Range de coeficientes aerodinâmicos
        CL = (2*L)/(rho*S*V**2)
        CD = CD0 + K*CL**2

        # Empuxo necessário para cruzeiro
        Thrust = 0.5*rho*V**2*S*CD

        # Combustível consumido
        delta_fuel = TSFC*Thrust*dt #[N]

        # Alcance específico
        SR = V/(delta_fuel/dt)

        # Velocidade de LRC - condição onde alcance específico é igual a 99% do
        # máximo alcance específico

        LRC, idx = get_long_range_cruise_velocity(V,SR)

        # Atualização do peso
        W = W - delta_fuel[idx]

        # Step no tempo e na distância
        t += dt
        x += LRC*dt

        # Informações desejadas
        t_plot.append(t/60)
        x_plot.append(x/1000)
        delta_fuel_plot.append(delta_fuel[idx])

    return t_plot, x_plot, delta_fuel_plot


def max_endurance_cruise(rho, S, CD0, K, W0, TSFC, t_total):
    # Range de velocidades
    V = np.linspace(10, 500, 100)

    # Inicializando variáveis

    t = 0
    dt = 1
    W = W0
    fuel_consumed = 0
    
    t_plot = []
    fuel_spent_plot = []
    
    t_plot.append(0)
    fuel_spent_plot.append(0)
    

    while t < t_total:
        # Coeficientes aerodinâmicos
        CL = (2*W)/(rho*S*V**2)
        CD = CD0 + K*CL**2

        # Condição de máxima autonomia
        idx = np.nanargmax(CL**(3/2)/CD)
        V_max_endurance = V[idx]

        # Empuxo necessário
        Thrust = 0.5*rho*V_max_endurance**2*S*CD[idx]

        # Combustível consumido
        dfuel = TSFC*Thrust*dt #[N]

        # Atualização do peso
        W -= dfuel

        # Step no tempo e no combustível
        t += dt
        fuel_consumed += dfuel

        # Informações desejadas
        t_plot.append(t/60) # min
        fuel_spent_plot.append(dfuel)
    
    return t_plot, fuel_spent_plot


def get_long_range_cruise_velocity(V,SR):

    SR_auxiliar = SR[:]
    error = abs(SR-0.99*max(SR))
    index_first_min = np.nanargmin(error)

    SR_auxiliar[index_first_min] = np.nan
    error = abs(SR_auxiliar-0.99*max(SR_auxiliar))
    index_second_min = np.nanargmin(error)

    V_min = V[index_first_min]
    V_second_min = V[index_second_min]
    vector_velocities = [V_min,V_second_min]
    LRC = max(vector_velocities)
    idx_vel = np.nanargmax(vector_velocities)
    idx = np.where(V==vector_velocities[idx_vel])[0][0]

    return LRC, idx

# Função para plot

def show_results(x,y,xlabel, ylabel, title):
    plt.figure()
    plt.plot(x,y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid()