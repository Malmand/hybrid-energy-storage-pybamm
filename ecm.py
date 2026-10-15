"""Equivalent-circuit cell and pack models (one RC branch) built from step 2 parameters."""
import numpy as np


class ECMCell:
    """Cell with OCV(SOC), R0(SOC) and one RC branch. Current positive = discharge."""

    def __init__(self, m, soc0=0.5, capacity_Ah=None):
        self.m = m
        self.Q = capacity_Ah or m["capacity_Ah"]
        self.soc, self.v1 = soc0, 0.0

    def params(self):
        m, s = self.m, self.soc
        return (np.interp(s, m["ocv_soc"], m["ocv_v"]),
                np.interp(s, m["soc"], m["R0_ohm_Ah"]) / self.Q,  # ohm*Ah -> ohm for this capacity
                np.interp(s, m["soc"], m["R1_ohm_Ah"]) / self.Q,
                np.interp(s, m["soc"], m["tau_s"]))

    def step_current(self, i, dt):
        ocv, r0, r1, tau = self.params()
        v = ocv - self.v1 - r0 * i
        a = np.exp(-dt / tau)  # exact RC update for constant current
        self.v1 = self.v1 * a + r1 * i * (1 - a)
        self.soc -= i * dt / 3600 / self.Q  # coulomb counting
        return v


class ECMPack:
    """n_series identical cells; capacity set by the pack energy. Power positive = discharge."""

    def __init__(self, m, energy_kWh, n_series, soc0=0.5, soc_min=0.1, soc_max=0.95):
        self.m, self.ns = m, n_series
        v_nom = float(np.interp(0.5, m["ocv_soc"], m["ocv_v"]))
        self.Q = energy_kWh * 1000 / (n_series * v_nom) if energy_kWh > 0 else 0.0  # Ah per cell string
        self.cell = ECMCell(m, soc0, self.Q) if self.Q > 0 else None
        self.soc_min, self.soc_max = soc_min, soc_max

    @property
    def soc(self):
        return self.cell.soc if self.cell else np.nan

    def max_power(self, discharge=True):
        if not self.cell:
            return 0.0
        ocv, r0, _, _ = self.cell.params()
        e = (ocv - self.cell.v1) * self.ns
        r = r0 * self.ns
        v_lim = (self.m["v_min"] if discharge else self.m["v_max"]) * self.ns
        p = v_lim * (e - v_lim) / r               # power at the voltage limit
        if discharge and self.soc <= self.soc_min:
            return 0.0
        if not discharge and self.soc >= self.soc_max:
            return 0.0
        return max(p, 0.0) if discharge else max(-p, 0.0)

    def step_power(self, p_w, dt):
        """Apply a power request (clipped to limits). Returns delivered power, current, voltage, loss."""
        if not self.cell or p_w == 0:
            if self.cell:
                self.cell.step_current(0.0, dt)
            return 0.0, 0.0, 0.0, 0.0
        p = float(np.clip(p_w, -self.max_power(False), self.max_power(True)))
        ocv, r0, r1, _ = self.cell.params()
        e = (ocv - self.cell.v1) * self.ns
        r = r0 * self.ns
        disc = e * e - 4 * r * p  # solve p = (e - r i) i for i
        i = (e - np.sqrt(max(disc, 0.0))) / (2 * r)  # smaller root: physical branch
        v_cell = self.cell.step_current(i, dt)
        loss = (i * i * r0 + (self.cell.v1 ** 2) / r1) * self.ns  # ohmic + RC branch losses
        return p, i, v_cell * self.ns, loss


class SuperCapModule:
    """Supercapacitor bank: n_series x n_parallel cells, each an ideal capacitor with ESR.

    Cell values are typical of a 3000 F / 2.7 V cell (ESR 0.29 mOhm). The usable window is
    half to full rated voltage, which holds 75% of the stored energy. Power positive = discharge.
    """
    C_CELL, ESR_CELL, V_CELL = 3000.0, 0.29e-3, 2.7

    def __init__(self, n_series, n_parallel, v0_frac=0.8):
        self.ns, self.np_ = n_series, n_parallel
        self.C = self.C_CELL * n_parallel / n_series if n_parallel else 0.0  # series divides C
        self.R = self.ESR_CELL * n_series / n_parallel if n_parallel else 0.0
        self.v_max = self.V_CELL * n_series
        self.v_min = 0.5 * self.v_max
        self.v = v0_frac * self.v_max

    @property
    def energy_usable_kWh(self):
        return 0.5 * self.C * (self.v_max ** 2 - self.v_min ** 2) / 3.6e6

    @property
    def soc(self):
        """Share of the usable energy window that is filled (0 at half voltage, 1 at rated)."""
        return (self.v ** 2 - self.v_min ** 2) / (self.v_max ** 2 - self.v_min ** 2) if self.C else np.nan

    def step_power(self, p_w, dt):
        if not self.C:
            return 0.0, 0.0, 0.0, 0.0
        # energy-based limits for this step
        e_avail = 0.5 * self.C * (self.v ** 2 - self.v_min ** 2)
        e_room = 0.5 * self.C * (self.v_max ** 2 - self.v ** 2)
        p = float(np.clip(p_w, -e_room / dt, e_avail / dt))
        disc = self.v ** 2 - 4 * self.R * p  # p = (v - R i) i
        i = (self.v - np.sqrt(max(disc, 0.0))) / (2 * self.R)
        loss = i * i * self.R
        # energy balance: d(C v^2 / 2) = -(p + loss) dt
        self.v = float(np.sqrt(max(self.v ** 2 - 2 * (p + loss) * dt / self.C, 0.0)))
        return p, i, self.v, loss
