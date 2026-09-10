'''
This script holds input dictionaries for several aircraft
which can be used as starting point for new designs.
'''

# GENERAL IMPORTS
import numpy as np
from .constants import ft2m, gravity, nm2m

#=================================================

def standard_airplane(name='fokker100'):
    '''
    The standard parameters refer to the Fokker 100, but they could be redefined for
    any new aircraft.
    '''

    if name == 'fokker100':

        # This model was taken from measuring the 3D view from
        # https://www.icas.org/ICAS_ARCHIVE/ICAS1988/ICAS-88-1.6.2.pdf
        # Obert, E. The Aerodynamic Development of the Fokker 100
        
        inputs = {'type': 'transport', # Can be 'transport', 'fighter', or 'general'
                    
                  'S_w' : 93.5, # Wing area [m2] - From Obert's paper
                  'AR_w' : 8.32,  # Wing aspect ratio
                  'taper_w' : 0.25, # Wing taper ratio
                  'sweep_w' : 15.76*np.pi/180, # Wing sweep [rad]
                  'dihedral_w' : 3*np.pi/180, # Wing dihedral [rad]
                  'xr_w' : 13.38, # Longitudinal position of the wing (with respect to the fuselage nose) [m]
                  'zr_w' : -0.99, # Vertical position of the wing (with respect to the fuselage nose) [m]
                  'tcr_w' : 0.123, # t/c of the root section of the wing
                  'tct_w' : 0.096, # t/c of the tip section of the wing
                    
                  'Cht' : 0.976, # Horizontal tail volume coefficient
                  'Lc_h' : 4.23, # Non-dimensional lever of the horizontal tail (lever/wing_mac)
                  'AR_h' : 4.52, # HT aspect ratio
                  'taper_h' : 0.44, # HT taper ratio
                  'sweep_h' : 26*np.pi/180, # HT sweep [rad]
                  'dihedral_h' : 2*np.pi/180, # HT dihedral [rad]
                  'zr_h' : 4.76, # Vertical position of the HT [m]
                  'tcr_h' : 0.1, # t/c of the root section of the HT
                  'tct_h' : 0.1, # t/c of the tip section of the HT
                  'eta_h' : 1.0, # Dynamic pressure factor of the HT (0.9 for conventional tail or 1.0 for T-tail)
                    
                  'Cvt' : 0.068, # Vertical tail volume coefficient
                  'Lb_v' : 0.469, # Non-dimensional lever of the vertical tail (lever/wing_span)
                  'AR_v' : 1.03, # VT aspect ratio
                  'taper_v' : 0.73, # VT taper ratio
                  'sweep_v' : 39*np.pi/180, # VT sweep [rad]
                  'zr_v' : 1.39, # Vertical position of the VT [m]
                  'tcr_v' : 0.1, # t/c of the root section of the VT
                  'tct_v' : 0.1, # t/c of the tip section of the VT
                    
                  'L_f' : 31.6, # Fuselage length [m]
                  'D_f' : 3.28, # Fuselage diameter [m]
                    
                  'x_n' : 21.4, # Longitudinal position of the nacelle frontal face [m]
                  'y_n' : 3.01, # Lateral position of the nacelle centerline [m]
                  'z_n' : 0.45, # Vertical position of the nacelle centerline [m]
                  'L_n' : 4.91, # Nacelle length [m]
                  'D_n' : 1.69, # Nacelle diameter [m]
                    
                  'n_engines' : 2, # Number of engines
                  'n_engines_under_wing' : 0, # Number of engines installed under the wing
                  'engine' : {'model' : 'Howe turbofan', # Check engineTSFC function for options
                              #'model' : 'Raymer turbofan', # Check engineTSFC function for options
                              'BPR' : 3.04, # Engine bypass ratio
                              #'weight' : 1500*gravity, # Single engine weight [N] (Can also be omitted to let designTool estimate it)
                              'C_ref' : 0.70/3600, # Reference thrust-specific fuel consumption [1/s] (Can also be omitted to let designTool estimate it)
                              'altitude_ref': 35000*ft2m, # Altitude that corresponds to the given TSFC [m]
                              'Mach_ref': 0.73, # Mach that corresponds to the given TSFC
                              },
                    
                  'x_nlg' : 3.7, # Longitudinal position of the nose landing gear [m]
                  'x_mlg' : 17.4, # Longitudinal position of the main landing gear [m]
                  'y_mlg' : 2.47, # Lateral position of the main landing gear [m]
                  'z_lg' : -4, # Vertical position of the landing gear [m]
                  'x_tailstrike' : 23.4, # Longitudinal position of critical tailstrike point [m]
                  'z_tailstrike' : -1.54, # Vertical position of critical tailstrike point [m]
                    
                  'c_tank_c_w' : 0.4, # Fraction of the wing chord occupied by the fuel tank
                  'x_tank_c_w' : 0.2, # Fraction of the wing chord where fuel tank starts
                  'b_tank_b_w_start' : 0.0, # Fraction of the wing semi-span where fuel tank starts
                  'b_tank_b_w_end' : 0.95, # Fraction of the wing semi-span where fuel tank ends
                    
                  'clmax_w' : 1.8, # Maximum lift coefficient of wing airfoil
                  'k_korn' : 0.91, # Airfoil technology factor for Korn equation (wave drag)
        
                  'flap_type' : 'double slotted',  # Flap type
                  'c_flap_c_wing' : 0.30, # Fraction of the wing chord occupied by flaps
                  'b_flap_b_wing' : 0.60, # Fraction of the wing span occupied by flaps (including fuselage portion)
                    
                  'slat_type' : None, # Slat type
                  'c_slat_c_wing' : 0.0, # Fraction of the wing chord occupied by slats
                  'b_slat_b_wing' : 0.0, # Fraction of the wing span occupied by slats

                  'c_ail_c_wing' : 0.27, # Fraction of the wing chord occupied by aileron
                  'b_ail_b_wing' : 0.34, # Fraction of the wing span occupied by aileron
                    
                  'h_ground' : 35.0*ft2m, # Distance to the ground for ground effect computation [m]
                  'k_exc_drag' : 0.06, # Excrescence drag factor

                  'winglet' : False, # Add winglet
                    
                  'altitude_takeoff' : 0.0, # Altitude for takeoff computation [m] - From Obert's paper
                  'distance_takeoff' : 2050.0, # Required takeoff distance [m] - From Obert's paper
                  'deltaISA_takeoff' : 15.0, # Variation from ISA standard temperature [ºC] - From Obert's paper
                    
                  'altitude_landing' : 0.0, # Altitude for landing computation [m]
                  'distance_landing' : 1340.0, # Required landing distance [m]
                  'deltaISA_landing' : 0.0, # Variation from ISA standard temperature [ºC]
                  'MLW_frac' : 40100/43090, # Max Landing Weight / Max Takeoff Weight - From Obert's paper
                    
                  'altitude_cruise' : 35000*ft2m, # Cruise altitude for design mission [m] - From Obert's paper
                  'Mach_cruise' : 0.73, # Cruise Mach number for design mission - From Obert's paper
                  'range_cruise' : 1310*nm2m, # Cruise range for design mission [m] - From Obert's paper

                  'altitude_maxcruise' : 35000*ft2m, # Altitude for high-speed cruise [m]
                  'Mach_maxcruise' : 0.77, # Mach for high-speed cruise [m]
                    
                  'time_loiter' : 45*60, # Loiter time [s]
                  'altitude_loiter' : 1500*ft2m, # Loiter time [s] - DOI: 10.2514/6.2025-3572
                    
                  'altitude_altcruise' : 15000*ft2m, # Alternative cruise altitude [m] - based on: DOI: 10.2514/6.2025-3572
                  'Mach_altcruise' : 0.4, # Alternative cruise Mach number (limitation of 250 kts)
                  'range_altcruise' : 200*nm2m, # Alternative cruise range [m]
                    
                  'W_payload' : 107*91*gravity, # Payload weight for design mission [N]
                  'xcg_payload' : 14.4, # Longitudinal position of the Payload center of gravity [m]
                  'W_maxpayload' : 11242*gravity, # Maximum payload weight [N]
                    
                  'W_crew' : 5*91*gravity, # Crew weight [N]
                  'xcg_crew' : 2.5, # Longitudinal position of the Crew center of gravity [m]

                  'block_range' : 400*nm2m, # Block range [m]
                  'block_time' : (1.0 + 2*40/60)*3600, # Block time [s]
                  'n_captains' : 1, # Number of captains in flight
                  'n_copilots' : 1, # Number of copilots in flight
                    
                  'rho_fuel' : 804, # Fuel density kg/m3 (This is Jet A-1)

                  'W0_guess' : 40000*gravity # Guess for MTOW
                  }

    elif name == 'my_airplane':

        # v3 (Lab 02): aeronave escolhida pela equipe -- ponto B, o joelho da
        # frente de Pareto W0 x Wf obtida com NSGA-II (lab2_opt_equipe_moga.py).
        # W0 = 291.291,1 kgf e Wf = 109.891,3 kgf, ou seja -5,04% de MTOW e
        # -10,47% de combustivel em relacao a v1 (configuracao final do PRJ-22).
        # Em relacao a v2 (otimo mono-objetivo de min W0, 290.117,1 kgf), B
        # aceita +0,4% de MTOW para queimar -1,0% de combustivel por missao.
        # A v1 continua disponivel em standard_airplane('my_airplane_v1'),
        # pois os scripts do Lab 02 partem dela.

        inputs = {'type': 'transport', # Can be 'transport', 'fighter', or 'general'

                  'S_w' : 368.833, # Wing area [m2] (v1: 390.0)
                  'AR_w' : 9.80998,  # Wing aspect ratio (v1: 8.0)
                  'taper_w' : 0.24, # Wing taper ratio
                  'sweep_w' : 0.60478, # Wing sweep [rad] (v1: 0.58)
                  'dihedral_w' : 6*np.pi/180, # Wing dihedral [rad]
                  'xr_w' : 15.9966, # Longitudinal position of the wing (with respect to the fuselage nose) [m] (v1: 17)
                  'zr_w' : -1.3, # Vertical position of the wing (with respect to the fuselage nose) [m]
                  'tcr_w' : 0.196226, # t/c of the root section of the wing (v1: 0.18)
                  'tct_w' : 0.08, # t/c of the tip section of the wing

                  'Cht' : 0.7, # Horizontal tail volume coefficient
                  'Lc_h' : 4.0, # Non-dimensional lever of the horizontal tail (lever/wing_mac)
                  'AR_h' : 4.6, # HT aspect ratio
                  'taper_h' : 0.45, # HT taper ratio
                  'sweep_h' : 0.66, # HT sweep [rad]
                  'dihedral_h' : 0.03490658503988659, # HT dihedral [rad]
                  'zr_h' : 2, # Vertical position of the HT [m]
                  'tcr_h' : 0.1, # t/c of the root section of the HT
                  'tct_h' : 0.1, # t/c of the tip section of the HT
                  'eta_h' : 0.9, # Dynamic pressure factor of the HT (0.9 for conventional tail or 1.0 for T-tail)

                  'Cvt' : 0.085, # Vertical tail volume coefficient
                  'Lb_v' : 0.51, # Non-dimensional lever of the vertical tail (lever/wing_span)
                  'AR_v' : 1.5, # VT aspect ratio
                  'taper_v' : 0.4, # VT taper ratio
                  'sweep_v' : 0.66, # VT sweep [rad]
                  'zr_v' : 3.0, # Vertical position of the VT [m]
                  'tcr_v' : 0.1, # t/c of the root section of the VT
                  'tct_v' : 0.1, # t/c of the tip section of the VT

                  'L_f' : 61.39, # Fuselage length [m]
                  'D_f' : 6.08, # Fuselage diameter [m]

                  'x_n' : 18, # Longitudinal position of the nacelle frontal face [m]
                  'y_n' : 11, # Lateral position of the nacelle centerline [m]
                  'z_n' : -3, # Vertical position of the nacelle centerline [m]
                  'L_n' : 7.3, # Nacelle length [m]
                  'D_n' : 4, # Nacelle diameter [m]

                  'n_engines' : 2, # Number of engines
                  'n_engines_under_wing' : 2, # Number of engines installed under the wing
                  'engine' : {'model' : 'Howe turbofan', # Check engineTSFC function for options
                              'BPR' : 9, # Engine bypass ratio
                              'C_ref' : 0.00013888888889, # Reference thrust-specific fuel consumption [1/s] (Can also be omitted to let designTool estimate it)
                              'altitude_ref': 10668.0, # Altitude that corresponds to the given TSFC [m]
                              'Mach_ref': 0.85, # Mach that corresponds to the given TSFC
                              'weight' : 86000, # Single engine weight [N] (Can also be omitted to let designTool estimate it)
                              'Tmax' : 513000, # Maximum thrust of a single engine [N]
                              },

                  'x_nlg' : 3.0, # Longitudinal position of the nose landing gear [m]
                  'x_mlg' : 29.3738, # Longitudinal position of the main landing gear [m] (v1: 31)
                  'y_mlg' : 6.93599, # Lateral position of the main landing gear [m] (v1: 5.5)
                  'z_lg' : -6.02342, # Vertical position of the landing gear [m] (v1: -5.75)
                  'x_tailstrike' : 46.19, # Longitudinal position of critical tailstrike point [m]
                  'z_tailstrike' : -3.04, # Vertical position of critical tailstrike point [m]

                  'c_tank_c_w' : 0.41, # Fraction of the wing chord occupied by the fuel tank
                  'x_tank_c_w' : 0.12, # Fraction of the wing chord where fuel tank starts
                  'b_tank_b_w_start' : 0.0, # Fraction of the wing semi-span where fuel tank starts
                  'b_tank_b_w_end' : 0.9, # Fraction of the wing semi-span where fuel tank ends

                  'clmax_w' : 1.8, # Maximum lift coefficient of wing airfoil
                  'k_korn' : 0.95, # Airfoil technology factor for Korn equation (wave drag)

                  'flap_type' : 'double slotted',  # Flap type
                  'c_flap_c_wing' : 0.3, # Fraction of the wing chord occupied by flaps
                  'b_flap_b_wing' : 0.6, # Fraction of the wing span occupied by flaps (including fuselage portion)

                  'slat_type' : 'slat', # Slat type
                  'c_slat_c_wing' : 0.15, # Fraction of the wing chord occupied by slats
                  'b_slat_b_wing' : 0.7, # Fraction of the wing span occupied by slats

                  'c_ail_c_wing' : 0.27, # Fraction of the wing chord occupied by aileron
                  'b_ail_b_wing' : 0.34, # Fraction of the wing span occupied by aileron

                  'h_ground' : 10.668000000000001, # Distance to the ground for ground effect computation [m]
                  'k_exc_drag' : 0.06, # Excrescence drag factor

                  'winglet' : True, # Add winglet

                  'altitude_takeoff' : 0.0, # Altitude for takeoff computation [m]
                  'distance_takeoff' : 2900.0, # Required takeoff distance [m]
                  'deltaISA_takeoff' : 0.0, # Variation from ISA standard temperature [ºC]

                  'altitude_landing' : 0.0, # Altitude for landing computation [m]
                  'distance_landing' : 2900.0, # Required landing distance [m]
                  'deltaISA_landing' : 0.0, # Variation from ISA standard temperature [ºC]
                  'MLW_frac' : 0.9306103504293339, # Max Landing Weight / Max Takeoff Weight

                  'altitude_cruise' : 10668.0, # Cruise altitude for design mission [m]
                  'Mach_cruise' : 0.85, # Cruise Mach number for design mission
                  'range_cruise' : 14816000.0, # Cruise range for design mission [m]

                  'altitude_maxcruise' : 10668.0, # Altitude for high-speed cruise [m]
                  'Mach_maxcruise' : 0.87, # Mach for high-speed cruise [m]

                  'time_loiter' : 2700, # Loiter time [s]
                  'altitude_loiter' : 457.20000000000005, # Loiter altitude [m]

                  'altitude_altcruise' : 4572, # Alternative cruise altitude [m]
                  'Mach_altcruise' : 0.4, # Alternative cruise Mach number (limitation of 250 kts)
                  'range_altcruise' : 370400.0, # Alternative cruise range [m]

                  'W_payload' : 313920.0, # Payload weight for design mission [N]
                  'xcg_payload' : 29.357, # Longitudinal position of the Payload center of gravity [m]

                  'W_crew' : 13390.650000000001, # Crew weight [N]
                  'xcg_crew' : 25, # Longitudinal position of the Crew center of gravity [m]

                  'block_range' : 740800.0, # Block range [m]
                  'block_time' : 8399.999999999998, # Block time [s]
                  'n_captains' : 1, # Number of captains in flight
                  'n_copilots' : 1, # Number of copilots in flight

                  'rho_fuel' : 804, # Fuel density kg/m3 (This is Jet A-1)

                  'W0_guess' : 2_857_566 # Guess for MTOW (= W0 da aeronave B, 291.291,1 kgf)
                  }

    elif name == 'my_airplane_v1':

        # v1: configuracao final do PRJ-22. E o ponto de partida das
        # otimizacoes do Lab 02 e a baseline das figuras de comparacao, entao
        # fica registrada aqui para que os scripts lab2_*.py continuem
        # reproduzindo o relatorio depois que 'my_airplane' passou a ser a v3.
        # W0 = 306.745,2 kgf e Wf = 122.743,5 kgf.
        #
        # A v1 difere da v3 apenas nas oito variaveis de projeto do Lab 02
        # (e no chute de MTOW), entao e escrita como um override da v3 para
        # nao duplicar o dicionario e nao permitir que as duas versoes se
        # desencontrem no resto dos parametros.

        inputs = standard_airplane('my_airplane')['inputs']

        inputs.update({'AR_w' : 8.0,   # Wing aspect ratio
                       'xr_w' : 17,    # Longitudinal position of the wing [m]
                       'S_w' : 390.0,  # Wing area [m2]
                       'sweep_w' : 0.58, # Wing sweep [rad] 0.5759586531581288
                       'x_mlg' : 31,   # Longitudinal position of the MLG [m]
                       'tcr_w' : 0.18, # t/c of the root section of the wing
                       'z_lg' : -5.75, # Vertical position of the landing gear [m]
                       'y_mlg' : 5.5,  # Lateral position of the MLG [m]
                       'W0_guess' : 2_844_900, # Guess for MTOW
                       })

    airplane = {'inputs':inputs}

    return airplane