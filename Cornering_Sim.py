class Vehicle:
    def __init__(self, m, fwr, wr, p_mult, final_dr, st, cg_h, fa, wheelbase, Cd, track_width=None):
        self.vMass = m
        self.Wf_ratio = fwr
        self.wheel_radius = wr
        self.power_mult = p_mult
        self.fdr = final_dr
        self.shift_delay_s = st
        self.h_cg = cg_h
        self.frontal_area = fa
        self.wheelbase = wheelbase
        self.track_width = track_width
        self.Cd = Cd
    
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


class TireModel:
    def max_lateral_force(self, Fz, **kwargs):
        """Return max lateral force for a given normal load."""
        raise NotImplementedError

class MuTireModel(TireModel):
    def __init__(self, mu):
        self.mu = mu
    
    def max_lateral_force(self, Fz, **kwargs):
        return self.mu * Fz


class Environment:
    def __init__(self, g=9.81, rho=1.225):
        self.g = g
        self.rho = rho

import math

class Simulator:
    def __init__(self, vehicle, tire_model, environment):
        self.car = vehicle
        self.tire_model = tire_model
        self.env = environment

    def max_corner_speed(self, radius, initial_guess=10.0):
        """Compute steady-state max corner speed for a given radius."""
        v = initial_guess
        for _ in range(30):  # simple fixed-point iteration
            # Aero downforce
            downforce = 0.5 * self.env.rho * self.car.frontal_area * self.car.Cd * v**2
            # Normal load (simplified, divided across 4 wheels)
            Fz_total = self.car.vMass * self.env.g + downforce
            F_lat_max = self.tire_model.max_lateral_force(Fz_total / 4.0) * 4  # sum all wheels

            v_new = math.sqrt(F_lat_max * radius / self.car.vMass)
            if abs(v_new - v) < 1e-3:
                break
            v = v_new
        return v





if __name__ == "__main__":
    G25 = Vehicle(288, .47, 0.203, 1, 40/13, 0.18, 0.3, .95, 1.530, 0.7)
    G24 = Vehicle(298, .53, 0.26, 0.9, 40/11, 0.2, 0.3, .95, 1.530, 0.7)
    print(G25)
    print(G25)

    # Create environment
    env = Environment()

    # Create tire model
    tire = MuTireModel(mu=1.6)

    # Create vehicle
    car = Vehicle(
        m=300, fwr=0.45, wr=0.25, p_mult=1.0,
        final_dr=4.0, st=0.05, cg_h=0.25, fa=1.2,
        wheelbase=2.0, Cd=1.2, track_width=1.5
    )

    # Create simulator
    sim = Simulator(vehicle=car, tire_model=tire, environment=env)

    # Compute max corner speed
    radius = 50  # meters
    v_max = sim.max_corner_speed(radius)
    print(f"Max corner speed for radius {radius} m: {v_max:.2f} m/s")
