"""HESS plant: a Li-ion battery pack (reduced ECM from step 2) and a supercapacitor bank.

Used directly by step 4 and wrapped as an FMU (fmu/HessFmu.py) for MATLAB/Simulink or any
FMI 2.0 co-simulation master.
"""
import json
import os

from ecm import ECMPack, SuperCapModule

HERE = os.path.dirname(os.path.abspath(__file__))


def load_ecm(name="liion", folder=None):
    with open(os.path.join(folder or os.path.join(HERE, "models"), f"ecm_{name}.json")) as fh:
        return json.load(fh)


class HessPlant:
    def __init__(self, bat_kWh=1000.0, bat_series=200, sc_series=150, sc_parallel=0, soc0=0.5, model_folder=None):
        self.bat = ECMPack(load_ecm("liion", model_folder), bat_kWh, bat_series, soc0=soc0)
        self.sc = SuperCapModule(sc_series, sc_parallel)

    def step(self, p_bat_w, p_sc_w, dt):
        pb, ib, vb, lb = self.bat.step_power(p_bat_w, dt)
        ps_, is_, vs, ls = self.sc.step_power(p_sc_w, dt)
        return dict(p_bat=pb, i_bat=ib, v_bat=vb, loss_bat=lb, soc_bat=self.bat.soc,
                    p_sc=ps_, i_sc=is_, v_sc=vs, loss_sc=ls, soc_sc=self.sc.soc)
