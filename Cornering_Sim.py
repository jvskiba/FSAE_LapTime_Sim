import math
from Tire_Model import *
import numpy as np
import matplotlib.pyplot as plt

class Vehicle:
    def __init__(self, m, fwr, wr, p_mult, final_dr, shift_time, cg_h, fa, wheelbase, Cd, track_width):
        self.vMass = m
        self.Wf_ratio = fwr
        self.wheel_radius = wr
        self.power_mult = p_mult
        self.fdr = final_dr
        self.shift_delay_s = shift_time
        self.h_cg = cg_h
        self.frontal_area = fa
        self.wheelbase = wheelbase
        self.track_width = track_width
        self.Cd = Cd

        self.a = self.Wf_ratio * self.vMass * self.wheelbase / (self.vMass)   # front CG distance
        self.b = self.wheelbase - self.a                                 # rear CG distance
        self.Iz = self.vMass * (self.wheelbase**2 + self.track_width**2) / 12.0

    
    def __str__(self):
        return (
            f"Vehicle:\n"
            f"  Mass: {self.vMass:.1f} kg\n"
            f"  Weight dist (front): {self.Wf_ratio*100:.1f}%\n"
            f"  Wheel radius: {self.wheel_radius:.3f} m\n"
            f"  Power multiplier: {self.power_mult}\n"
            f"  Final drive ratio: {self.fdr}\n"
            f"  Shift delay: {self.shift_delay_s:.3f} s\n"
            f"  CG height: {self.h_cg:.3f} m\n"
            f"  Frontal area: {self.frontal_area:.2f} m^2\n"
            f"  Wheelbase: {self.wheelbase}\n"
            f"  Track width: {self.track_width}\n"
        )
#End


class Environment:
    def __init__(self, g=9.81, rho=1.225):
        self.g = g
        self.rho = rho


class Simulator:
    def __init__(self, vehicle, tire_model, environment, dt=0.01):
        self.car = vehicle
        self.tire_model = tire_model
        self.env = environment
        self.dt = dt

        # States
        self.u = 0.0  # longitudinal velocity [m/s]
        self.v = 0.0  # lateral velocity [m/s]
        self.r = 0.0  # yaw rate [rad/s]

    def weight_transfer(self, ax, ay):
        """
        Compute per-wheel loads from inertial weight transfer (mass only).
        """
        mass = self.car.vMass
        g = self.env.g
        h_cg = self.car.h_cg
        wheelbase = self.car.wheelbase
        track_width = self.car.track_width
        fwr = self.car.Wf_ratio

        # Static distribution (mass only, no aero)
        W_total = mass * g
        W_front = fwr * W_total
        W_rear  = (1 - fwr) * W_total
        
        W_fl = W_front / 2.0
        W_fr = W_front / 2.0
        W_rl = W_rear  / 2.0
        W_rr = W_rear  / 2.0

        # Lateral transfer (left tires gain, right tires lose)
        dFy = (mass * ay * h_cg) / track_width
        W_fl += dFy / 2.0
        W_rl += dFy / 2.0
        W_fr -= dFy / 2.0
        W_rr -= dFy / 2.0

        # Longitudinal transfer (front loses under accel, gains under braking)
        dFx = (mass * ax * h_cg) / wheelbase
        W_fl -= dFx / 2.0
        W_fr -= dFx / 2.0
        W_rl += dFx / 2.0
        W_rr += dFx / 2.0

        return {"FL": W_fl, "FR": W_fr, "RL": W_rl, "RR": W_rr}

    def downforce(self, v):
        return 0, 0

    def calc_tire_loads(self, v, ax, ay):
        tire_loads = self.weight_transfer(ax, ay)
        F_front, F_rear = self.downforce(v)

        tire_loads["FL"] += F_front / 2.0
        tire_loads["FR"] += F_front / 2.0
        tire_loads["RL"] += F_rear / 2.0
        tire_loads["RR"] += F_rear / 2.0

        return tire_loads
    
    def reset(self, u=10.0, v=0.0, r=0.0):
        self.u = u
        self.v = v
        self.r = r

if __name__ == "__main__":
    # Vehicle, environment, tire model
    env = Environment()
    tire = UGAMSTireModel()

    #VMass, FWeightRatio, Wheel Radius, Power Mult, FDR, Shift Time, CG Height, Frontal Area, wheelbase, Cd, track_width
    G25 = Vehicle(288, .47, 0.203, 1, 40/13, 0.18, 0.3, .95, 1.530, 0.7, 1.368)
    G24 = Vehicle(298, .53, 0.26, 0.9, 40/11, 0.2, 0.3, .95, 1.530, 0.7, 1.368) #Check Stacks

    # Simulator
    sim = Simulator(vehicle=G25, tire_model=tire, environment=env)


