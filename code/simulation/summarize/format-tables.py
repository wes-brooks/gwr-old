#! /Library/Frameworks/Python.framework/Versions/Current/bin/python
import argparse
import sys
import re

#Parse command line arguments
parser = argparse.ArgumentParser(description='Format a table for publication.')
parser.add_argument('--infile', type=str, nargs='?', default="")
parser.add_argument('--outfile', type=str, nargs='?', default="")
parser.add_argument('--sizex', type=str, nargs='?', default="")
parser.add_argument('--sizey', type=str, nargs='?', default="")
parser.add_argument('--tabcolsep', type=str, nargs='?', default="")
parser.add_argument('--header', type=str, nargs='?', default="")

#Interpret the command line arguments:
args = vars(parser.parse_args())
infile = args['infile']
outfile = args['outfile']
tabcolsep = str(args['tabcolsep'])
header = str(args['header'])
if args['sizex']: sizex = str(args['sizex'])
else: sizex = "!"
if args['sizey']: sizey = str(args['sizey'])
else: sizey = "!"


#Import the header file if one is specified
if header:
    with open(header) as header_file:
        header = header_file.read()
else: header = ""


#Extract the content of the table to a list
lines = []
with open(infile) as input_file:
    for i, l in enumerate(input_file):
        if re.search(r'\\begin{table}', l) or re.search(r'\\begin{sidewaystable}', l):
            lines.append(l)
            lines.append("\\thispagestyle{empty}\n")
            lines.append("\\begin{center}\n")
            lines.append("\\resizebox{" + sizex + "}{" + sizey + "}{\n")
            if (tabcolsep): lines.append("\\setlength\\tabcolsep{" + tabcolsep + "}\n")
        elif re.search(r'\\centering', l):
            pass
        elif re.search(r'\\begin{tabular}', l):
            lines.append(l)
            lines.append(header)
        elif re.search(r'\\end{tabular}', l):
            l = re.sub(r'\\end{tabular}', r'\end{tabular}}', l)
            lines.append(l)
        elif re.search(r'\\caption{.*}', l):
            lines.append(l)
            lines.append("\\end{center}\n")
        else: lines.append(l)

#Write the altered table to the file.
out = open(outfile, "w")
for line in lines:
    out.write(line)

#Done!
out.close()


