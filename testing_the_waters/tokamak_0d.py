import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

# Physical Constants
E_CHARGE = 1.60217663e-19  # Joules per eV
KEV_TO_J = E_CHARGE * 1e3  # Joules per keV
MU_0 = 4 * np.pi * 1e-7    # Vacuum permeability (H/m)

@dataclass
class ReactorConfig:
    """Geometric and magnetic configuration of the toroidal device."""
    name: str = "SPARC-Class Tokamak"
    R: float = 1.85          # Major radius (m)
    a: float = 0.57          # Minor radius (m)
    B0: float = 12.2         # Toroidal magnetic field on axis (T)
    Ip: float = 8.7          # Plasma current (MA) - Set to 0 for pure Stellarator
    kappa: float = 1.97      # Elongation
    delta: float = 0.54      # Triangularity
    H_factor: float = 1.0    # Confinement enhancement factor (H98(y,2))
    Z_eff: float = 1.5       # Effective ion charge (impurity measure)
    A_mass: float = 2.5      # Average atomic mass (2.5 for 50:50 D-T)
    is_stellarator: bool = False
    iota_bar: float = 0.9    # Rotational transform (used if is_stellarator=True)

    @property
    def aspect_ratio(self) -> float:
        return self.R / self.a

    @property
    def volume(self) -> float:
        """Approximate plasma volume (m^3) accounting for elongation."""
        return 2.0 * (np.pi ** 2) * self.R * (self.a ** 2) * self.kappa

    @property
    def surface_area(self) -> float:
        """Approximate first-wall surface area (m^2)."""
        return 4.0 * (np.pi ** 2) * self.R * self.a * np.sqrt((1.0 + self.kappa**2) / 2.0)


class PlasmaEngine0D:
    """Zero-dimensional steady-state power balance and stability evaluator."""

    def __init__(self, config: ReactorConfig):
        self.cfg = config

    @staticmethod
    def bosch_hale_dt(T_keV: np.ndarray) -> np.ndarray:
        """
        Bosch-Hale parameterization for Maxwellian D-T reactivity <sigma*v> in m^3/s.
        Valid for 0.2 keV < T < 100 keV.
        """
        T = np.clip(T_keV, 0.2, 100.0)
        # Bosch-Hale coefficients for D(t,n)alpha
        BG = 34.3827
        mc2 = 1124656.0
        C1 = 1.17302e-9
        C2 = 1.51361e-2
        C3 = 7.51886e-2
        C4 = 4.60643e-3
        C5 = 1.35000e-2
        C6 = -1.06750e-4
        C7 = 1.36600e-5

        numer = T * (C2 + T * (C4 + T * C6))
        denom = 1.0 + T * (C3 + T * (C5 + T * C7))
        theta = T / (1.0 - (numer / denom))
        xi = ((BG ** 2) / (4.0 * theta)) ** (1.0 / 3.0)

        sigmav_cm3 = C1 * theta * np.sqrt(xi / (mc2 * (T ** 3))) * np.exp(-3.0 * xi)
        return sigmav_cm3 * 1e-6  # Convert cm^3/s to m^3/s

    def calc_safety_factor_q95(self) -> float:
        """
        Cylindrical/Shaped safety factor q_95 approximation for tokamaks,
        or inverse rotational transform for stellarators.
        """
        if self.cfg.is_stellarator or self.cfg.Ip == 0:
            return 1.0 / self.cfg.iota_bar
        
        eps = self.cfg.a / self.cfg.R
        # Shaping factor for elongation and triangularity
        f_shape = (1.0 + (self.cfg.kappa ** 2) * (1.0 + 2.0 * (self.cfg.delta ** 2) - 1.2 * (self.cfg.delta ** 3))) / 2.0
        q_star = (5.0 * (self.cfg.a ** 2) * self.cfg.B0 / (self.cfg.R * self.cfg.Ip)) * f_shape
        return q_star

    def calc_greenwald_limit(self) -> float:
        """Greenwald density limit n_G in units of 10^20 m^-3."""
        if self.cfg.is_stellarator:
            # Sudo density limit approximation for stellarators (assuming 10 MW absorbed power baseline)
            P_abs = 10.0 
            return 0.25 * np.sqrt((P_abs * self.cfg.B0) / (self.cfg.a**2 * self.cfg.R))
        return self.cfg.Ip / (np.pi * (self.cfg.a ** 2))

    def evaluate_state(self, n20: np.ndarray, T_keV: np.ndarray):
        """
        Evaluates 0D plasma parameters across a scalar or 2D meshgrid of:
        n20: Electron density in units of 10^20 m^-3
        T_keV: Volume-averaged temperature (assuming T_i = T_e) in keV
        """
        ne = n20 * 1e20  # m^-3
        # Fuel dilution slightly from Z_eff (assuming carbon/low-Z impurity Z=6)
        n_dt = ne * np.clip((6.0 - self.cfg.Z_eff) / 5.0, 0.5, 1.0)
        
        # 1. Fusion & Alpha Power (MW)
        sigmav = self.bosch_hale_dt(T_keV)
        E_fusion = 17.6e3 * KEV_TO_J  # 17.6 MeV per reaction
        E_alpha = 3.52e3 * KEV_TO_J   # 3.52 MeV to alpha heating
        
        # 50-50 D-T mixture: n_D = n_T = 0.5 * n_dt -> n_D * n_T = 0.25 * n_dt^2
        reaction_rate = 0.25 * (n_dt ** 2) * sigmav * self.cfg.volume
        P_fus_MW = (reaction_rate * E_fusion) * 1e-6
        P_alpha_MW = (reaction_rate * E_alpha) * 1e-6

        # 2. Bremsstrahlung Radiation Loss (MW)
        # Standard coefficient ~ 5.35e-37 W * m^3 / keV^0.5
        P_brem_MW = 5.35e-37 * (ne ** 2) * self.cfg.Z_eff * np.sqrt(T_keV) * self.cfg.volume * 1e-6

        # 3. Stored Thermal Energy W_th (MJ)
        # W = 3 * n_e * k_B * T * V (assuming n_i ~ n_e and T_i = T_e)
        W_th_MJ = 3.0 * ne * (T_keV * KEV_TO_J) * self.cfg.volume * 1e-6

        # 4. Energy Confinement Time tau_E (s) via ITER IPB98(y,2) Scaling
        # Note: IPB98(y,2) depends on loss power P_loss. We solve self-consistently:
        # tau_E = C * P_loss^-0.69, and P_loss = W_th / tau_E -> tau_E = C^(1 / (1 - 0.69)) * W_th^(-0.69 / 0.31)
        n19 = n20 * 10.0
        Ip_eff = self.cfg.Ip if not self.cfg.is_stellarator else (5.0 * self.cfg.a**2 * self.cfg.B0 * self.cfg.iota_bar / self.cfg.R)
        
        C_ipb = (0.0562 * self.cfg.H_factor * (Ip_eff ** 0.93) * (self.cfg.B0 ** 0.15) *
                 (n19 ** 0.41) * (self.cfg.R ** 1.97) * ((self.cfg.a / self.cfg.R) ** 0.58) *
                 (self.cfg.kappa ** 0.78) * (self.cfg.A_mass ** 0.19))
        
        # Self-consistent algebra eliminating P_loss circular dependency:
        # tau_E = C_ipb * P_loss^-0.69 = C_ipb * (W_th / tau_E)^-0.69
        # tau_E^(1 - 0.69) = C_ipb * W_th^-0.69
        tau_E = (C_ipb * (np.maximum(W_th_MJ, 1e-3) ** -0.69)) ** (1.0 / 0.31)

        # 5. Conductive Power Loss & Required Auxiliary Power (MW)
        P_cond_MW = W_th_MJ / tau_E
        P_loss_total_MW = P_cond_MW + P_brem_MW
        P_aux_MW = np.maximum(P_loss_total_MW - P_alpha_MW, 0.0)

        # 6. Fusion Gain Q
        Q = np.where(P_aux_MW > 0.01, P_fus_MW / P_aux_MW, 1000.0) # 1000 represents ignition

        # 7. MHD Beta & Operational Limits
        # Volume-averaged plasma pressure (Pa) = 2 * n_e * k_B * T
        pressure_Pa = 2.0 * ne * (T_keV * KEV_TO_J)
        B_pressure_Pa = (self.cfg.B0 ** 2) / (2.0 * MU_0)
        beta_toroidal_pct = (pressure_Pa / B_pressure_Pa) * 100.0
        
        # Normalized Beta: beta_N = beta(%) / (Ip / (a * B0))
        beta_N = beta_toroidal_pct / (Ip_eff / (self.cfg.a * self.cfg.B0))

        # Neutron Wall Loading (80% of fusion power is in 14.1 MeV neutrons)
        neutron_wall_load = (0.8 * P_fus_MW) / self.cfg.surface_area

        return {
            "P_fus_MW": P_fus_MW,
            "P_alpha_MW": P_alpha_MW,
            "P_brem_MW": P_brem_MW,
            "P_cond_MW": P_cond_MW,
            "P_aux_MW": P_aux_MW,
            "Q": Q,
            "tau_E": tau_E,
            "beta_pct": beta_toroidal_pct,
            "beta_N": beta_N,
            "q95": self.calc_safety_factor_q95(),
            "n_Greenwald": self.calc_greenwald_limit(),
            "neutron_wall_load_MW_m2": neutron_wall_load
        }


def plot_popcon(engine: PlasmaEngine0D):
    """Generates a 2D Plasma Operation Contour (POPCON) plot."""
    n_grid = np.linspace(0.2, 5.0, 200)   # 10^20 m^-3
    T_grid = np.linspace(1.0, 30.0, 200)  # keV
    N, T = np.meshgrid(n_grid, T_grid)

    res = engine.evaluate_state(N, T)

    fig, ax = plt.subplots(figsize=(10, 7), dpi=120)
    plt.style.use('dark_background')
    ax.set_facecolor('#0d1117')
    fig.patch.set_facecolor('#0d1117')

    # 1. Shaded heatmap of Fusion Power
    cp = ax.contourf(T, N, res["P_fus_MW"], levels=30, cmap="inferno", alpha=0.85)
    cbar = fig.colorbar(cp, ax=ax)
    cbar.set_label("Total Fusion Power $P_{fus}$ (MW)", fontsize=11)

    # 2. Auxiliary Heating Power Contours (P_aux)
    aux_levels = [5, 15, 30, 60, 100]
    c_aux = ax.contour(T, N, res["P_aux_MW"], levels=aux_levels, colors='cyan', linestyles='dashed', linewidths=1.2)
    ax.clabel(c_aux, inline=True, fontsize=9, fmt="$P_{aux}=%d$ MW")

    # 3. Fusion Gain Q Contours (Breakeven Q=1, Burning Q=5, Ignition)
    q_levels = [1.0, 5.0, 10.0, 50.0]
    c_q = ax.contour(T, N, res["Q"], levels=q_levels, colors='lime', linewidths=2.0)
    ax.clabel(c_q, inline=True, fontsize=10, fmt="$Q=%g$")

    # 4. Operational Limits (Greenwald Density & Troyon Beta Limit)
    n_G = res["n_Greenwald"]
    ax.axhline(n_G, color='red', linestyle='-', linewidth=2.5, label=f"Greenwald Limit ($n_G = {n_G:.2f} \\times 10^{{20}} m^{{-3}}$)")
    ax.fill_between(T_grid, n_G, 5.0, color='red', alpha=0.2, hatch='//')

    # Troyon Limit: beta_N = 2.8
    c_beta = ax.contour(T, N, res["beta_N"], levels=[2.8], colors='magenta', linewidths=2.5)
    ax.clabel(c_beta, inline=True, fontsize=10, fmt="Troyon Limit ($\\beta_N = 2.8$)")

    ax.set_title(f"POPCON Operational Space: {engine.cfg.name}\n($R={engine.cfg.R}$m, $a={engine.cfg.a}$m, $B_0={engine.cfg.B0}$T, $q_{{95}}={res['q95']:.2f}$)", fontsize=13, pad=12)
    ax.set_xlabel("Volume-Averaged Temperature $\\langle T \\rangle$ (keV)", fontsize=11)
    ax.set_ylabel("Electron Density $\\langle n_e \\rangle$ ($10^{20}$ m$^{-3}$)", fontsize=11)
    ax.set_ylim(0.2, 5.0)
    ax.set_xlim(1.0, 30.0)
    ax.grid(True, linestyle=':', alpha=0.3)
    ax.legend(loc="upper left", facecolor="#161b22", edgecolor="white")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Initialize a SPARC-class high-field compact tokamak
    sparc_cfg = ReactorConfig()
    engine = PlasmaEngine0D(sparc_cfg)

    # Evaluate a single design point: n_e = 3.0e20 m^-3, T = 10.0 keV
    point = engine.evaluate_state(n20=np.array(3.0), T_keV=np.array(10.0))
    
    print(f"=== {sparc_cfg.name} Design Point Evaluation ===")
    print(f"Safety Factor (q95):       {point['q95']:.2f} (Stable > 2.0)")
    print(f"Greenwald Limit (n_G):     {point['n_Greenwald']:.2f} x 10^20 m^-3")
    print(f"Confinement Time (tau_E):  {point['tau_E']:.3f} s")
    print(f"Fusion Power (P_fus):      {point['P_fus_MW']:.1f} MW")
    print(f"Alpha Heating (P_alpha):   {point['P_alpha_MW']:.1f} MW")
    print(f"Bremsstrahlung (P_brem):   {point['P_brem_MW']:.1f} MW")
    print(f"Required Aux Power:        {point['P_aux_MW']:.1f} MW")
    print(f"Fusion Gain (Q):           {point['Q']:.2f}")
    print(f"Normalized Beta (beta_N):  {point['beta_N']:.2f} (Limit ~ 2.8)")
    print(f"Neutron Wall Loading:      {point['neutron_wall_load_MW_m2']:.2f} MW/m^2")

    # Generate POPCON Visualization
    plot_popcon(engine)