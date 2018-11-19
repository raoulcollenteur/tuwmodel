import pandas as pd
from numpy import array
import matplotlib.pyplot as plt

prec = pd.read_csv("data/prec.csv", header=0, squeeze=True,
                   infer_datetime_format=True, index_col="Date")
evap = pd.read_csv("data/evap.csv", header=0, squeeze=True,
                   infer_datetime_format=True, index_col="Date")
river = pd.read_csv("data/river.csv", header=0, squeeze=True,
                    infer_datetime_format=True, index_col="Date")
airt = pd.read_csv("data/temp.csv", infer_datetime_format=True,
                   squeeze=True, index_col=0,
                   usecols=["Time [YYYY-MM-DD HH:MM:SS]",
                            "Air temperature [degC]"])
#"SCF","DDF","Tr","Ts","Tm","LPrat","FC","BETA","k0","k1","k2","lsuz",
# "cperc","bmax","croute"


params = array([1.2, # snow correction factor [-] (e.g., 0.9-1.5);
                1.2, # DDF degree day factor [mm/degC/timestep] (e.g., 0.0-5.0 mm/degC/day);
                2, #Tr threshold temperature above which precipitation is rain [degC] (e.g., 1.0-3.0 degC);
                -2, # Ts threshold temperature below which precipitation is snow [degC] (e.g., -3.0-1.0 degC);
                0, #Tm threshold temperature above which melt starts [degC] (e.g., -2.0-2.0 degC);
                0.9, # LPrat parameter related to the limit for potential evaporation [-] (e.g., 0.0-1.0);
                60, # FC field capacity, i.e., max soil moisture storage [mm](e.g., 0-600 mm);
                3.3, # BETA the non linear parameter for runoff production [-] (e.g., 0.0-20.0);
                100, # k0 storage coefficient for very fast response [
                # timestep] (e.g., 0.0-2.0 days);
                1, # k1 storage coefficient for fast response [timestep] (
                # e.g., 2.0-30.0 days);
                5, # k2 storage coefficient for slow response [timestep] (
                # e.g., 30.0-250.0 days);
                1, # lsuz threshold storage state, i.e., the very fast response start if exceeded [mm] (e.g., 1.0-100.0 mm);
                2, # cperc constant percolation rate [mm/timestep] (e.g.,
                # 0.0-8.0 mm/day);
                1, # bmax maximum base at low flows [timestep] (e.g.,
                # 0.0-30.0 days);
                26.5 #croute free scaling parameter [timestep^2/mm] (e.g., 0.0-50.0 days^2/mm);
                ]).reshape(15, 1)
incon = array([50, 0, 2.5, 2.5]).reshape(4, 1)


import tuwmodel as tuw

data = tuw.simulate(area=1, param=params, incon=incon, prec=prec,
                    airt=airt, evap=evap, tmin="2007", tmax="2012")


(data["qzones"]*15000).plot()
river.loc["2007":"2012"].plot()
plt.ylim(0, 120)