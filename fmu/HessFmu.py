"""FMI 2.0 co-simulation wrapper of the HESS plant, built with pythonfmu.

Inputs  : P_bat_set, P_sc_set (W, positive = discharge)
Outputs : P_bat, P_sc (delivered, W), SOC_bat, SOC_sc, V_bat, V_sc (V), P_loss (W)
Params  : bat_kWh, bat_series, sc_series, sc_parallel, soc0

Build:  pythonfmu build -f fmu/HessFmu.py hess_plant.py ecm.py models --dest fmu
The FMU runs the Python plant inside the FMI interface, so the host needs Python with numpy.
"""
import os
import sys

from pythonfmu import Fmi2Causality, Fmi2Slave, Fmi2Variability, Real


class HessFmu(Fmi2Slave):
    author = "Mohamed Al Mandhari"
    description = "Li-ion battery (reduced from PyBaMM) + supercapacitor bank"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        sys.path.insert(0, self.resources)
        self.P_bat_set = self.P_sc_set = 0.0
        self.P_bat = self.P_sc = self.P_loss = 0.0
        self.SOC_bat = self.SOC_sc = self.V_bat = self.V_sc = 0.0
        self.bat_kWh, self.bat_series, self.sc_series, self.sc_parallel, self.soc0 = 1000.0, 200.0, 150.0, 40.0, 0.5
        self.plant = None
        for n in ("P_bat_set", "P_sc_set"):
            self.register_variable(Real(n, causality=Fmi2Causality.input))
        for n in ("P_bat", "P_sc", "SOC_bat", "SOC_sc", "V_bat", "V_sc", "P_loss"):
            self.register_variable(Real(n, causality=Fmi2Causality.output))
        for n in ("bat_kWh", "bat_series", "sc_series", "sc_parallel", "soc0"):
            self.register_variable(Real(n, causality=Fmi2Causality.parameter, variability=Fmi2Variability.tunable))

    def exit_initialization_mode(self):
        from hess_plant import HessPlant
        self.plant = HessPlant(self.bat_kWh, int(self.bat_series), int(self.sc_series), int(self.sc_parallel),
                               self.soc0, model_folder=os.path.join(self.resources, "models"))

    def do_step(self, current_time, step_size):
        r = self.plant.step(self.P_bat_set, self.P_sc_set, step_size)
        self.P_bat, self.P_sc = float(r["p_bat"]), float(r["p_sc"])
        self.SOC_bat, self.SOC_sc = float(r["soc_bat"]), float(r["soc_sc"])
        self.V_bat, self.V_sc = float(r["v_bat"]), float(r["v_sc"])
        self.P_loss = float(r["loss_bat"] + r["loss_sc"])
        return True
