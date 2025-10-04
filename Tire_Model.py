import numpy as np

class TireModel:
    def max_lateral_force(self, SA, FZ, P, IA):
        """Return max lateral force for a given normal load."""
        raise NotImplementedError

class MuTireModel(TireModel):
    def __init__(self, mu):
        self.mu = mu
    
    def max_lateral_force(self, SA, FZ, P=12, IA=0):
        return self.mu * FZ

class UGAMSTireModel(TireModel):
    """
    UGA Motorsports Tire Model Magic Formula 6.1.2
    Based on OptimumTire coefficients for Hoosier R20 43075 16 x 7.5
    Equations from Hans Pacejka - Tire and Vehicle Dynamics (3rd Ed.)
    """

    def __init__(self):
        # Lateral Coefficients
        self.PCY1 = 1.47938
        self.PDY1 = 2.30264
        self.PDY2 = -0.181517
        self.PDY3 = 12.2537
        self.PEY1 = 0.258174
        self.PEY2 = 0.134876
        self.PEY3 = 0.964615
        self.PEY4 = 0.0467532
        self.PEY5 = -156.583
        self.PKY1 = 33.3683
        self.PKY2 = 3.92018
        self.PKY3 = 0.963332
        self.PKY4 = 5.07208
        self.PKY5 = 49.6379
        self.PKY6 = -3.67812
        self.PKY7 = -1.32934
        self.PHY1 = -0.000546807
        self.PHY2 = 0.00152636
        self.PVY1 = -0.00246833
        self.PVY2 = -0.0198409
        self.PVY3 = -0.0718706
        self.PVY4 = 1.50852
        self.PPY1 = 0.290986
        self.PPY2 = 0.971999
        self.PPY3 = -0.185848
        self.PPY4 = 0.19433
        self.PPY5 = -0.351119

        # Nominal Values
        self.FNOMIN = -1110
        self.NOMPRES = 14 * 6894.76  # Pa

        # Scaling Coefficients
        self.LFZO = 1
        self.LCY = 1
        self.LMUY = 0.66
        self.LKY = 1
        self.LKYC = 1
        self.LVY = 1
        self.LHY = 1
        self.LEY = 1

        # Spin correction (set to 1 if neglected)
        self.ZETA0 = 1
        self.ZETA2 = 1
        self.ZETA3 = 1
        self.ZETA4 = 1

    def evaluate(self, SA, FZ, P, IA):
        """Evaluate the tire model for slip angle (SA), normal load (FZ),
        pressure (P in psi), and inclination angle (IA in rad)."""

        gammaAST = np.sin(IA)

        # (4.E1) Reference load
        Fz0_ = self.LFZO * self.FNOMIN

        # (4.E2a) Load difference
        dfz = (FZ - Fz0_) / Fz0_

        # (4.E2b) Pressure difference
        dpi = (P * 6894.76 - self.NOMPRES) / self.NOMPRES

        # (4.E21) Cornering stiffness factor
        Cy = self.PCY1 * self.LCY

        # (4.E23) Friction coefficient
        Muy = (self.PDY1 + self.PDY2 * dfz) * \
              (1 + self.PPY3 * dpi + self.PPY4 * dpi**2) * \
              (1 - self.PDY3 * gammaAST**2) * self.LMUY

        # (4.E22) Peak factor
        Dy = Muy * FZ * self.ZETA2

        # (4.E25) Cornering stiffness
        Kya = (self.PKY1 * Fz0_ * (1 + self.PPY1 * dpi) *
               (1 - self.PKY3 * np.abs(gammaAST)) *
               np.sin(self.PKY4 * np.arctan((FZ / Fz0_) /
               ((self.PKY2 + self.PKY5 * gammaAST**2) * (1 + self.PPY2 * dpi)))) *
               self.ZETA3 * self.LKY)

        # (4.E24) Curvature factor
        Ey = ((self.PEY1 + self.PEY2 * dfz) *
              (1 + self.PEY5 * gammaAST**2 -
               (self.PEY3 + self.PEY4 * gammaAST) * np.sign(SA)) * self.LEY)

        # (4.E39) Avoid division by zero
        signKya = np.sign(Kya)
        if signKya == 0:
            signKya = 1
        Kya_ = Kya + np.finfo(float).eps * signKya

        # (4.E28) Camber-induced side force
        SVyg = FZ * (self.PVY3 + self.PVY4 * dfz) * gammaAST * \
               self.ZETA2 * self.LKYC * self.LMUY

        # (4.E30) Camber stiffness
        Kyy0 = FZ * (self.PKY6 + self.PKY7 * dfz) * \
               (1 + self.PPY5 * dpi) * self.LKYC

        # (4.E29) Lateral shift
        SVy = FZ * (self.PVY1 + self.PVY2 * dfz) * \
              self.LVY * self.LMUY * self.ZETA2 + SVyg

        # (4.E27) Horizontal shift
        SHy = ((self.PHY1 + self.PHY2 * dfz) * self.LHY +
               (Kyy0 * gammaAST - SVyg) / Kya_ * self.ZETA0 +
               self.ZETA4 - 1)

        # (4.E20) Effective slip angle
        alphay = SA + SHy

        # (4.E26) Stiffness factor
        signCy = np.sign(Cy)
        if signCy == 0:
            signCy = 1
        By = Kya / (Cy * Dy + np.finfo(float).eps * signCy)

        # (4.E19) Final lateral force
        FY = Dy * np.sin(Cy * np.arctan(By * alphay -
              Ey * (By * alphay - np.arctan(By * alphay)))) + SVy

        return {
            "FY": FY,
            "Muy": Muy,
            "dfz": dfz,
            "Fz0_": Fz0_,
            "dpi": dpi,
            "By": By,
            "Cy": Cy,
            "Ey": Ey,
            "Dy": Dy,
            "Kya_": Kya_,
            "SVy": SVy,
            "SHy": SHy
        }
    
    def max_lateral_force(self, SA, FZ, P=12, IA=0):
        return self.evaluate(SA, FZ, P, IA)["FY"]
