import re
import sys

root = str(sys.argv[1]) + "/"
cluster = str(sys.argv[2])


for method in ['glmnet', 'enet', 'oracle']:
    for i in range(1200):
        try:
            filename = root + "output/MiscParams." + cluster + "." + str(i) + "." + method + ".csv"
            f = open(filename, 'r')
            contents = f.read()
            f.close()

            contents = re.sub("c\([\d.,\s-]*,\s*(?P<final>-?[0-9.]+)\n?\)", "\g<final>", contents)

            filename = root + "output/MiscParams." + cluster + "." + str(i) + "." + method + ".csv"
            f = open(filename, 'w')
            f.write(contents)
            f.close()
        except: pass