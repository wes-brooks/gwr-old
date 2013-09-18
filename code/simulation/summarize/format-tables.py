#! /Library/Frameworks/Python.framework/Versions/Current/bin/python
import argparse
import sys
import re

#Parse command line arguments
parser = argparse.ArgumentParser(description='Format a table for publication.')
parser.add_argument('--infile', type=str, nargs='?', default="")
parser.add_argument('--outfile', type=str, nargs='?', default="")
parser.add_argument('--boxsize', type=str, nargs='?', default="12cm")

args = vars(parser.parse_args())
infile = args['infile']
outfile = args['outfile']
boxsize = args['boxsize']

#Extract the content of the table to a list
lines = []
with open(infile) as input_file:
    for i, line in enumerate(input_file):
        if re.search("\\\\begin{table}", line):
            lines.append(line)
            lines.append("\\thispagestyle{empty}\n")
            lines.append("\\begin{center}\n")
            lines.append("\\resizebox{!}{" + boxsize + "}{\n")
        elif re.search("\\\\caption{.*}", line):
            line = re.sub("\\caption{(?P<caption>.*)}", "\\\\caption{\g<caption>}}\n", line)
            lines.append(line)
            lines.append("\\end{center}\n")
        else: lines.append(line)

#Write the altered table to the file.
out = open(outfile, "w")
for line in lines:
    out.write(line)

out.close()


