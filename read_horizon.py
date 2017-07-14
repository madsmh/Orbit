import numpy as np

def readdata(planet):
    path = 'data/' + planet + '.csv'
    data = np.loadtxt(path, delimiter=',', dtype=np.float64, skiprows=1, usecols=range(2, 8))
    data = data*1e3
    return data


def readdiagnosticdata(planet):
    path = 'diagnostic_data/' + planet + '.csv'
    data = np.loadtxt(path, delimiter=',', dtype=np.float64, skiprows=1, usecols=range(2, 8))
    data = data * 1e3
    return data
