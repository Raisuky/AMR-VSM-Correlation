import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Compatibility for trapezoidal integration across NumPy versions
trapz_func = getattr(np, 'trapezoid', getattr(np, 'trapz', None))

# ==========================================
# 1. RELATIVE PATH CONFIGURATION
# ==========================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

file_agfm_long = os.path.join(DATA_DIR, "NiCuO_SP090626_4koe_3.txt")      # Longitudinal (0°)
file_agfm_trans = os.path.join(DATA_DIR, "NiCuO_SP090626_4Koe_plano.txt") # Transverse (90°)
file_mr = os.path.join(DATA_DIR, "Ni10CuO15.txt")                         # 4-Probe AMR

# ==========================================
# 2. DATA PROCESSING FUNCTIONS
# ==========================================
def load_and_correct_agfm(filepath, h_sat_threshold=2500.0):
    """Loads AGFM data and subtracts the linear diamagnetic glass substrate slope."""
    df = pd.read_csv(filepath, header=None, skiprows=1, sep=r'\s+')
    H_raw = df[0].values
    M_raw = df[1].values
    
    mask_pos = H_raw > h_sat_threshold
    mask_neg = H_raw < -h_sat_threshold
    
    if np.sum(mask_pos) > 5 and np.sum(mask_neg) > 5:
        slope_pos, intercept_pos = np.polyfit(H_raw[mask_pos], M_raw[mask_pos], 1)
        slope_neg, intercept_neg = np.polyfit(H_raw[mask_neg], M_raw[mask_neg], 1)
        chi_glass = 0.5 * (slope_pos + slope_neg)
        offset_m = 0.5 * (intercept_pos + intercept_neg)
        M_corr = M_raw - (chi_glass * H_raw) - offset_m
    else:
        M_corr = M_raw - np.mean(M_raw)
    return H_raw, M_corr

def get_coercive_fields(H, M):
    """Finds coercive fields (-Hc and +Hc) via linear interpolation around M = 0."""
    idx_min_H = np.argmin(H)
    H_down, M_down = H[:idx_min_H+1], M[:idx_min_H+1]
    mask_d = (H_down > -1000) & (H_down < 1000)
    Hc_down = np.interp(0.0, M_down[mask_d][::-1], H_down[mask_d][::-1])
    
    H_up, M_up = H[idx_min_H:], M[idx_min_H:]
    mask_u = (H_up > -1000) & (H_up < 1000)
    Hc_up = np.interp(0.0, M_up[mask_u], H_up[mask_u])
    return Hc_down, Hc_up

def get_anhysteretic_and_loop_area(H, m_norm, h_grid):
    """Interpolates down/up branches to compute enclosed hysteresis area and anhysteretic curve."""
    idx_min = np.argmin(H)
    H_down, m_down = H[:idx_min+1], m_norm[:idx_min+1]
    H_up, m_up = H[idx_min:], m_norm[idx_min:]
    
    sort_d = np.argsort(H_down)
    sort_u = np.argsort(H_up)
    m_down_interp = np.interp(h_grid, H_down[sort_d], m_down[sort_d])
    m_up_interp = np.interp(h_grid, H_up[sort_u], m_up[sort_u])
    
    loop_area_Oe = trapz_func(np.abs(m_down_interp - m_up_interp), h_grid)
    m_anhyst = 0.5 * (m_down_interp + m_up_interp)
    return m_down_interp, m_up_interp, m_anhyst, loop_area_Oe

# ==========================================
# 3. MAIN ANALYSIS PIPELINE
# ==========================================
if __name__ == "__main__":
    H_long, M_long_corr = load_and_correct_agfm(file_agfm_long, h_sat_threshold=2600.0)
    H_trans, M_trans_corr = load_and_correct_agfm(file_agfm_trans, h_sat_threshold=2400.0)

    m_long_norm = M_long_corr / np.max(np.abs(M_long_corr))
    m_trans_norm = M_trans_corr / np.max(np.abs(M_trans_corr))

    Hc_down_long, Hc_up_long = get_coercive_fields(H_long, m_long_norm)
    Hc_down_trans, Hc_up_trans = get_coercive_fields(H_trans, m_trans_norm)
    Hc_mean_long = 0.5 * (abs(Hc_down_long) + abs(Hc_up_long))
    Hc_mean_trans = 0.5 * (abs(Hc_down_trans) + abs(Hc_up_trans))
    H_eb = 0.5 * (Hc_down_long + Hc_up_long)

    # Thermodynamic Energy Density Integration
    M_s_Ni = 485.0  # Bulk Ni saturation magnetization at 300 K (emu/cm^3)
    H_grid_full = np.linspace(-3500, 3500, 1400)
    _, _, m_long_anhyst, area_loop_long_Oe = get_anhysteretic_and_loop_area(H_long, m_long_norm, H_grid_full)
    _, _, m_trans_anhyst, area_loop_trans_Oe = get_anhysteretic_and_loop_area(H_trans, m_trans_norm, H_grid_full)

    E_diss_long = area_loop_long_Oe * M_s_Ni
    E_diss_trans = area_loop_trans_Oe * M_s_Ni

    mask_q1 = (H_grid_full >= 0) & (H_grid_full <= 2500)
    H_q1 = H_grid_full[mask_q1]
    m_long_q1 = np.clip(m_long_anhyst[mask_q1], 0, 1.0)
    m_trans_q1 = np.clip(m_trans_anhyst[mask_q1], 0, 1.0)
    K_eff = trapz_func(np.maximum(m_trans_q1 - m_long_q1, 0), H_q1) * M_s_Ni

    # 4-Probe Magnetoresistance
    df_mr = pd.read_csv(file_mr, header=None, sep=r'\s+')
    H_mr = df_mr[1].values[1:-1]
    R_mr = df_mr[2].values[1:-1]

    idx_mr_minH = np.argmin(H_mr)
    H_mr_down, R_mr_down = H_mr[:idx_mr_minH+1], R_mr[:idx_mr_minH+1]
    H_mr_up, R_mr_up = H_mr[idx_mr_minH:], R_mr[idx_mr_minH:]

    R_perp = np.percentile(R_mr, 5)
    R_par = np.max(R_mr)
    Delta_R = R_par - R_perp
    AMR_ratio = (Delta_R / R_perp) * 100.0

    # Stoner-Wohlfarth Coherent Rotation Model
    H_k = 2200.0
    H_sim = np.linspace(-3500, 3500, 700)
    sin_theta_sw = np.clip((H_sim - H_eb) / H_k, -1.0, 1.0)
    R_sw_teorico = R_perp + Delta_R * (1.0 - sin_theta_sw**2)

    print("=== EXTRACTED PHYSICAL PARAMETERS ===")
    print(f"Longitudinal (0°):  -Hc = {Hc_down_long:.1f} Oe | +Hc = {Hc_up_long:.1f} Oe | <|Hc|> = {Hc_mean_long:.1f} Oe")
    print(f"Transverse (90°):   -Hc = {Hc_down_trans:.1f} Oe | +Hc = {Hc_up_trans:.1f} Oe | <|Hc|> = {Hc_mean_trans:.1f} Oe")
    print(f"Exchange Bias (Heb): {H_eb:.1f} Oe")
    print(f"K_eff:               {K_eff:.2e} erg/cm^3")
    print(f"E_diss (0° / 90°):   {E_diss_long:.2e} / {E_diss_trans:.2e} erg/cm^3")
    print(f"AMR States:          R_perp = {R_perp:.3f} Ohm | R_par = {R_par:.3f} Ohm | AMR = {AMR_ratio:.2f}%")

    # ==========================================
    # 4. PUBLICATION-GRADE FIGURES
    # ==========================================
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    # Figure 1: In-Plane Magnetic Anisotropy
    plt.figure(figsize=(7.5, 5.2))
    plt.fill_between(
        H_q1, m_long_q1, m_trans_q1, color='#1f77b4', alpha=0.18,
        label=fr'Anisotropy Energy ($K_{{\mathrm{{eff}}}} = {K_eff/1e5:.2f}\times 10^5$ erg/cm$^3$)'
    )
    plt.plot(H_long, m_long_norm, color='black', lw=2, label=r'Longitudinal $0^\circ$ (Hard Axis)')
    plt.plot(H_trans, m_trans_norm, color='#1f77b4', lw=2, label=r'Transverse $90^\circ$ (Easy Axis)')
    plt.axhline(0, color='gray', linestyle='--', alpha=0.6)
    plt.axvline(0, color='gray', linestyle='--', alpha=0.6)
    plt.title('In-Plane Magnetic Anisotropy', fontsize=13, fontweight='bold')
    plt.xlabel('Applied Magnetic Field $H$ (Oe)', fontsize=12)
    plt.ylabel('Normalized Magnetization $M / M_s$', fontsize=12)
    plt.xlim(-3500, 3500)
    plt.ylim(-1.15, 1.15)
    plt.legend(loc='lower right', fontsize=10, framealpha=0.95)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Plot1_Poster_Clean.png'), dpi=300)
    plt.close()

    # Figure 2: Magnetization Reversal and AMR Transport Correlation
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.0, 8.5), sharex=True)

    ax1.plot(H_long, m_long_norm, color='black', lw=2, label=r'AGFM Longitudinal ($0^\circ$)')
    ax1.axvline(Hc_down_long, color='#d62728', linestyle='--', lw=1.5, alpha=0.85, label=f'$-H_c = {Hc_down_long:.1f}$ Oe')
    ax1.axvline(Hc_up_long, color='#2ca02c', linestyle='--', lw=1.5, alpha=0.85, label=f'$+H_c = +{Hc_up_long:.1f}$ Oe')
    ax1.axvspan(Hc_down_long, Hc_up_long, color='#ff7f0e', alpha=0.15)
    ax1.axhline(0, color='gray', linestyle='-', alpha=0.5)
    ax1.set_ylabel(r'Normalized $M / M_s$', fontsize=12)
    ax1.set_title('Magnetization Reversal and AMR Transport Correlation', fontsize=13, fontweight='bold')
    ax1.legend(loc='lower right', fontsize=10, framealpha=0.95)
    ax1.set_ylim(-1.15, 1.15)

    ax2.plot(H_mr_down, R_mr_down, 'o-', color='#d62728', markersize=4, lw=1.2, alpha=0.85, label='Exp. AMR (Down-sweep)')
    ax2.plot(H_mr_up, R_mr_up, 's-', color='#2ca02c', markersize=4, lw=1.2, alpha=0.85, label='Exp. AMR (Up-sweep)')
    ax2.plot(H_sim, R_sw_teorico, color='#1f77b4', lw=2.5, label=f'Stoner-Wohlfarth Model ($H_k={H_k:.0f}$ Oe)')
    ax2.axvline(Hc_down_long, color='#d62728', linestyle='--', lw=1.5, alpha=0.85)
    ax2.axvline(Hc_up_long, color='#2ca02c', linestyle='--', lw=1.5, alpha=0.85)
    ax2.axvspan(
        Hc_down_long, Hc_up_long, color='#ff7f0e', alpha=0.18,
        label=fr'Switching Window ($\langle |H_c| \rangle = {Hc_mean_long:.1f}$ Oe)'
    )
    ax2.set_xlabel('Applied Magnetic Field $H$ (Oe)', fontsize=12)
    ax2.set_ylabel(r'Longitudinal Resistance $R_{xx}$ ($\Omega$)', fontsize=12)
    ax2.set_xlim(-3500, 3500)
    ax2.set_ylim(27.045, 27.175)
    ax2.legend(loc='upper right', fontsize=9.5, framealpha=0.95)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Plot2_Poster_Clean.png'), dpi=300)
    plt.close()
    print(f"Figures saved to {FIG_DIR}")
