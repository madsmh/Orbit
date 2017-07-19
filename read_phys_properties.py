import csv
def read_phys_properties():
    names = []
    radii = []
    gms = []

    with open('physical_properties/properties.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            names.append(row['NAME'])
            radii.append(float(row['RADIUS']))
            gms.append(float(row['GM']))
    return names, radii, gms


