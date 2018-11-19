from pandas import DataFrame
from numpy import array

from .hbvmodel import hbvmodel


def simulate(param=None, incon=None, prec=None, evap=None, airt=None, area=1,
             tmin=None, tmax=None):
    """Method to simulate the waterbalance using the TUWmodel.

    Parameters
    ----------
    param
    incon
    prec
    airt
    ep
    area

    Returns
    -------
    data: pandas.DataFrame
        DataFrame containing the different fluxes and system states.


    """
    prec = prec.loc[tmin:tmax]
    evap = evap.loc[tmin:tmax]
    airt = airt.loc[tmin:tmax]

    index = prec.index

    a = hbvmodel(area=area, param=param, incon=incon, prec=prec, airt=airt,
                 ep=evap)

    names = ["qzones", "swe", "moist", "rain", "snow", "melt", "q0", "q1",
             "q2", "eta", "suz", "slz"]

    data = DataFrame(a[0][:][:12].transpose(), columns=names, index=index)

    return data
