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
    name: str = "HELIAS-5B (W7-X Power Reactor Concept)"
    R: float = 22.0          # Major radius (m) [W7-X experiment is 5.5m]
    a: float = 1.8           # Effective minor radius (m) [W7-X experiment is 0.53m]
    B0: float = 5.6          # Magnetic field on axis (T) [W7-X experiment is 2.5T]
    Ip: float = 0.0          # Zero net plasma current in a Stellarator!
    kappa: float = 1.0       # Effective elongation (absorbed into effective minor radius 'a' for 3D stellarators)
    delta: float = 0.0       # Triangularity
    H_factor: float = 1.35   # Confinement renormalization factor (f_ren * H_ISS04 for optimized QI stellarator)
    Z_eff: float = 1.3       # Effective ion charge (impurity measure)
    A_mass: float = 2.5      # Average atomic mass (2.5 for 50:50 D-T)
    is_stellarator: bool = True
    iota_bar: float = 0.91   # Rotational transform at 2/3 radius (W7-X is ~ 5/5 = 1.0 at edge, ~0.86-0.91 core)
    beta_max_pct: float = 4.5 # W7-X / HELIAS equilibrium & ballooning stability limit (%)

    @property
    def aspect_ratio(self) -> float:
        return self.R / self.a

    @property
    def volume(self) -> float:
        """Approximate plasma volume (m^3)."""
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
        """Inverse rotational transform q = 1 / iota for stellarators."""
        if self.cfg.is_stellarator or self.cfg.Ip == 0:
            return 1.0 / self.cfg.iota_bar
        
        f_shape = (1.0 + (self.cfg.kappa ** 2) * (1.0 + 2.0 * (self.cfg.delta ** 2) - 1.2 * (self.cfg.delta ** 3))) / 2.0
        return (5.0 * (self.cfg.a ** 2) * self.cfg.B0 / (self.cfg.R * self.cfg.Ip)) * f_shape

    def evaluate_state(self, n20: np.ndarray, T_keV: np.ndarray):
        """
        Evaluates 0D plasma parameters across a scalar or 2D meshgrid of:
        n20: Electron density in units of 10^20 m^-3
        T_keV: Volume-averaged temperature (assuming T_i = T_e) in keV
        """
        ne = n20 * 1e20  # m^-3
        n_dt = ne * np.clip((6.0 - self.cfg.Z_eff) / 5.0, 0.5, 1.0)
        
        # 1. Fusion & Alpha Power (MW)
        sigmav = self.bosch_hale_dt(T_keV)
        E_fusion = 17.6e3 * KEV_TO_J
        E_alpha = 3.52e3 * KEV_TO_J
        
        reaction_rate = 0.25 * (n_dt ** 2) * sigmav * self.cfg.volume
        P_fus_MW = (reaction_rate * E_fusion) * 1e-6
        P_alpha_MW = (reaction_rate * E_alpha) * 1e-6

        # 2. Bremsstrahlung Radiation Loss (MW)
        P_brem_MW = 5.35e-37 * (ne ** 2) * self.cfg.Z_eff * np.sqrt(T_keV) * self.cfg.volume * 1e-6

        # 3. Stored Thermal Energy W_th (MJ)
        W_th_MJ = 3.0 * ne * (T_keV * KEV_TO_J) * self.cfg.volume * 1e-6

        # 4. Energy Confinement Time tau_E (s)
        n19 = n20 * 10.0
        if self.cfg.is_stellarator:
            # International Stellarator Scaling 2004 (ISS04):
            # tau_E = 0.134 * H * a^2.28 * R^0.64 * P_loss^-0.61 * n19^0.54 * B0^0.84 * iota^0.41
            C_iss04 = (0.134 * self.cfg.H_factor * (self.cfg.a ** 2.28) * (self.cfg.R ** 0.64) *
                       (n19 ** 0.54) * (self.cfg.B0 ** 0.84) * (self.cfg.iota_bar ** 0.41))
            # Self-consistent elimination of P_loss = W_th / tau_E:
            # tau_E = C_iss04 * (W_th / tau_E)^-0.61 -> tau_E^(1 - 0.61) = C_iss04 * W_th^-0.61
            tau_E = (C_iss04 * (np.maximum(W_th_MJ, 1e-3) ** -0.61)) ** (1.0 / 0.39)
        else:
            # Tokamak IPB98(y,2) scaling
            C_ipb = (0.0562 * self.cfg.H_factor * (self.cfg.Ip ** 0.93) * (self.cfg.B0 ** 0.15) *
                     (n19 ** 0.41) * (self.cfg.R ** 1.97) * ((self.cfg.a / self.cfg.R) ** 0.58) *
                     (self.cfg.kappa ** 0.78) * (self.cfg.A_mass ** 0.19))
            tau_E = (C_ipb * (np.maximum(W_th_MJ, 1e-3) ** -0.69)) ** (1.0 / 0.31)

        # 5. Conductive Power Loss & Required Auxiliary Power (MW)
        P_cond_MW = W_th_MJ / tau_E
        P_loss_total_MW = P_cond_MW + P_brem_MW
        P_aux_MW = np.maximum(P_loss_total_MW - P_alpha_MW, 0.0)

        # Total heating power sustaining the plasma (Alpha + Aux)
        P_heat_total_MW = P_alpha_MW + P_aux_MW

        # 6. Density Limit: Sudo Limit for Stellarators vs Greenwald for Tokamaks
        if self.cfg.is_stellarator:
            # Sudo density limit (10^20 m^-3) scales with sqrt(P_heat)
            # We compute n_Sudo across the grid; n > n_Sudo means radiative collapse
            n_limit = 0.25 * np.sqrt((np.maximum(P_heat_total_MW, 1.0) * self.cfg.B0) / ((self.cfg.a ** 2) * self.cfg.R))
        else:
            n_limit = np.full_like(n20, self.cfg.Ip / (np.pi * (self.cfg.a ** 2)))

        # 7. Fusion Gain Q
        Q = np.divide(P_fus_MW, P_aux_MW, out=np.full_like(P_fus_MW, 1000.0), where=(P_aux_MW > 0.01))

        # 8. MHD Beta (%)
        pressure_Pa = 2.0 * ne * (T_keV * KEV_TO_J)
        B_pressure_Pa = (self.cfg.B0 ** 2) / (2.0 * MU_0)
        beta_pct = (pressure_Pa / B_pressure_Pa) * 100.0

        neutron_wall_load = (0.8 * P_fus_MW) / self.cfg.surface_area

        return {
            "P_fus_MW": P_fus_MW,
            "P_alpha_MW": P_alpha_MW,
            "P_brem_MW": P_brem_MW,
            "P_cond_MW": P_cond_MW,
            "P_aux_MW": P_aux_MW,
            "Q": Q,
            "tau_E": tau_E,
            "beta_pct": beta_pct,
            "q95": self.calc_safety_factor_q95(),
            "n_limit": n_limit,
            "neutron_wall_load_MW_m2": neutron_wall_load
        }


def plot_popcon(engine: PlasmaEngine0D, design_point=None, save_path=None):
    """Generates a 2D Plasma Operation Contour (POPCON) plot with high-contrast dark theme."""
    # Apply dark background BEFORE creating figure so all labels default to crisp white!
    plt.style.use('dark_background')
    
    n_grid = np.linspace(0.2, 3.5, 200)   # 10^20 m^-3
    T_grid = np.linspace(1.0, 25.0, 200)  # keV
    N, T = np.meshgrid(n_grid, T_grid)

    res = engine.evaluate_state(N, T)

    fig, ax = plt.subplots(figsize=(10.5, 7.2), dpi=120)
    ax.set_facecolor('#0d1117')
    fig.patch.set_facecolor('#0d1117')

    # 1. Shaded heatmap of Fusion Power
    cp = ax.contourf(T, N, res["P_fus_MW"], levels=30, cmap="inferno", alpha=0.85)
    cbar = fig.colorbar(cp, ax=ax)
    cbar.set_label("Total Fusion Power $P_{fus}$ (MW)", fontsize=11, color="white")
    cbar.ax.yaxis.set_tick_params(color="white", labelcolor="white")

    # 2. Auxiliary Heating Power Contours (P_aux)
    aux_levels = [10, 30, 60, 120, 250]
    c_aux = ax.contour(T, N, res["P_aux_MW"], levels=aux_levels, colors='cyan', linestyles='dashed', linewidths=1.3)
    ax.clabel(c_aux, inline=True, fontsize=9, fmt="$P_{aux}=%d$ MW")

    # 3. Fusion Gain Q Contours
    q_levels = [1.0, 5.0, 20.0]
    c_q = ax.contour(T, N, res["Q"], levels=q_levels, colors='#00ff66', linewidths=2.2)
    ax.clabel(c_q, inline=True, fontsize=10, fmt="$Q=%g$")

    # 4. Stellarator Sudo Density Limit (where N - n_Sudo = 0)
    sudo_margin = N - res["n_limit"]
    # Draw the limit line only where the plasma still needs auxiliary heating (outside the ignited zone)
    sudo_line = np.where(res["P_aux_MW"] > 0.0, sudo_margin, np.nan)
    c_sudo = ax.contour(T, N, sudo_line, levels=[0.0], colors='#ff3333', linewidths=2.5)
    ax.clabel(c_sudo, inline=True, fontsize=10, fmt="Sudo Density Limit ($n_{Sudo}$)")
    # Shade unstable high-density radiative collapse zone
    ax.contourf(T, N, sudo_margin, levels=[0.0, 10.0], colors=['#ff3333'], alpha=0.18, hatches=['//'])

    # 5. Stellarator Equilibrium Beta Limit (<beta> = 4.5%)
    lvl = engine.cfg.beta_max_pct
    c_beta = ax.contour(T, N, res["beta_pct"], levels=[lvl], colors='#ff44ff', linewidths=2.5)
    ax.clabel(c_beta, inline=True, fontsize=10, fmt={lvl: f"Equilibrium Limit ($\\langle\\beta\\rangle = {lvl}\\%$)"})

    # 6. Ignited region (alpha heating alone exceeds losses, P_aux = 0).
    # A 0D steady state is not self-consistent here, so it is hatched rather than trusted.
    ignited = (res["P_aux_MW"] <= 0.0).astype(float)
    ax.contourf(T, N, ignited, levels=[0.5, 1.5], colors='none', hatches=['xx'])
    ax.contour(T, N, ignited, levels=[0.5], colors='white', linewidths=1.0, linestyles='dotted')

    # 7. Design point marker
    if design_point is not None:
        n_dp, T_dp = design_point
        ax.plot(T_dp, n_dp, marker='*', markersize=18, color='gold', markeredgecolor='black', zorder=10)
        ax.annotate("Design point", (T_dp, n_dp), xytext=(10, -16), textcoords='offset points',
                    color='gold', fontsize=10, fontweight='bold')

    # Crisp white labels and styling
    ax.set_title(
        f"Stellarator POPCON Space: {engine.cfg.name}\n"
        f"($R={engine.cfg.R}$m, $a={engine.cfg.a}$m, $B_0={engine.cfg.B0}$T, $\\iota_{{2/3}}={engine.cfg.iota_bar:.2f}$, $I_p=0$ MA)",
        fontsize=13, color="white", fontweight="bold", pad=14
    )
    ax.set_xlabel("Volume-Averaged Temperature $\\langle T \\rangle$ (keV)", fontsize=11, color="white", labelpad=8)
    ax.set_ylabel("Electron Density $\\langle n_e \\rangle$ ($10^{20}$ m$^{-3}$)", fontsize=11, color="white", labelpad=8)
    ax.tick_params(axis='both', colors='white', labelsize=10)
    
    ax.set_ylim(0.2, 3.5)
    ax.set_xlim(1.0, 25.0)
    ax.grid(True, linestyle=':', color='white', alpha=0.25)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, facecolor=fig.get_facecolor())
    plt.show()


if __name__ == "__main__":
    # Initialize HELIAS-5B Stellarator Reactor
    stellarator_cfg = ReactorConfig()
    engine = PlasmaEngine0D(stellarator_cfg)

    # Design point chosen to sit inside both operational limits with margin:
    # n_e = 1.6e20 m^-3 (~90% of n_Sudo), <T> = 10 keV (<beta> ~ 4.1%, ~91% of the 4.5% limit)
    n_dp, T_dp = 1.6, 10.0
    point = engine.evaluate_state(n20=np.array(n_dp), T_keV=np.array(T_dp))
    
    print(f"=== {stellarator_cfg.name} Design Point ===")
    print(f"Rotational Transform (iota): {stellarator_cfg.iota_bar:.2f} (q = {point['q95']:.2f})")
    print(f"Sudo Density Limit (n_Sudo): {point['n_limit']:.2f} x 10^20 m^-3 (Operating at {n_dp:.2f})")
    print(f"Confinement Time (tau_E):    {point['tau_E']:.2f} s (ISS04 Scaling)")
    print(f"Fusion Power (P_fus):        {point['P_fus_MW']:.1f} MW")
    print(f"Alpha Heating (P_alpha):     {point['P_alpha_MW']:.1f} MW")
    print(f"Bremsstrahlung (P_brem):     {point['P_brem_MW']:.1f} MW")
    print(f"Required Aux Power (P_aux):  {point['P_aux_MW']:.1f} MW")
    print(f"Fusion Gain (Q):             {point['Q']:.2f}")
    print(f"Volume-Avg Beta (<beta>):    {point['beta_pct']:.2f}% (Limit: {stellarator_cfg.beta_max_pct}%)")
    print(f"Neutron Wall Loading:        {point['neutron_wall_load_MW_m2']:.2f} MW/m^2")

    # Generate POPCON Visualization
    plot_popcon(engine, design_point=(n_dp, T_dp))