#! /usr/bin/python

import re
import sys

root = str(sys.argv[1]) + "/"

for method in ['glmnet', 'enet', 'oracle']:
    for i in range(900):
        try:
            filename = root + "/" + str(i+1) + "/MiscParams.poisson." + str(i) + "." + method + ".csv"
            f = open(filename, 'r')
            contents = f.read()
            f.close()

            contents = re.sub("c\([\d.,\s-]*,\s*(?P<final>-?[0-9.]+)\n?\)", "\g<final>", contents)

            filename = root + "/" + str(i+1) + "/MiscParams.poisson." + str(i) + "." + method + ".csv"
            f = open(filename, 'w')
            f.write(contents)
            f.close()
        except: pass