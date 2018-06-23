def calc_nis_percent(nis_file_name, nis_value):
    nis_file = open(nis_file_name, "r")
    above = 0
    below = 0
    for line in nis_file:
        if (float(line) >= nis_value):
            above += 1
        else:
            below += 1
    nis_file.close()
    
    return(above * 100.0 / (above + below))

nis_radar = 7.815
nis_lidar = 5.991

print("Radar: %.2f %%" % calc_nis_percent("../build/nis-radar.txt", nis_radar))
print("Lidar: %.2f %%" % calc_nis_percent("../build/nis-lidar.txt", nis_lidar))
