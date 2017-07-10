import csv

def readdata(planet):
    data = []
    path = 'data/' + planet + '.csv'
    with open(path) as csvfile:
        datareader = csv.reader(csvfile, delimiter=',', skipinitialspace=True)
        next(datareader)
        for row in datareader:
            # Convert srings to floats and to m/s and m.
            floatrow = [float(i)*1000 for i in row[2:8]]
            data.append(floatrow)
    return data


