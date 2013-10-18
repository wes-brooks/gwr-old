#! /Library/Frameworks/Python.framework/Versions/Current/bin/python

import argparse
import sys
import re
import operator

#Parse command line arguments
parser = argparse.ArgumentParser(description='Format a table for publication.')
parser.add_argument('order', type=str, nargs=3)
parser.add_argument('--infile', type=str, nargs='?', default="")
parser.add_argument('--outfile', type=str, nargs='?', default="")

args = vars(parser.parse_args())
order = args['order']
infile = args['infile']
outfile = args['outfile']

#Establish the simulation settings
taux = []
[taux.extend([x]*9) for x in [0.03, 0.1]]
rho = []
[rho.extend([x]*3) for x in [0, 0.5, 0.8]]
rho = rho * 2
taus = [0, 0.03, 0.1] * 6

#Extract the content of the table to a list
intable = False
outtable = False
header = []
rows = []
footer = []
with open(infile) as input_file:
    for i, line in enumerate(input_file):
        if not intable:
            if re.search("\\\\begin{tabular}", line):
                intable = True
                #line = re.sub("{ccccccccc}", "{ccc | ccc | ccc | ccc }", line)
            header.append(line)
        elif not outtable:
            if not re.search("\\\\end{tabular}", line): rows.append(line)
            else: outtable = True
        
        if outtable:
            footer.append(line)


multihead = ["&&&\\multicolumn{3}{c}{GWEN}&\\multicolumn{3}{c}{GWEN-LLE}&\\multicolumn{3}{c}{Oracle}\\\\", "MSE & bias & var &  MSE & bias & var &  MSE & bias & var\\\\ "]
rows[0] = multihead[1]

#Deal with the header: First, convert the column headers to latex code
code = ['taux', 'taus', 'rho']
latex = ['$\\tau_x$', '$\\tau_\\sigma$', '$\\rho$']
var_order = [latex[code.index(x)] for x in order]
var_header = " & ".join(var_order)
params = dict(zip(code, [taux, taus, rho]))

#Now prepend the variable names:
rows[0] = var_header + " & " + rows[0]
rowdiv = [len(set(params[order[x]])) for x in range(len(order))]

#Now layout the rest of the rows with the variable headers:
roworder = []
rownums = range(len(rows[2:]))

j = 0
vardivs = list(set(params[order[j]]))
vardivs.sort()
for vd in vardivs:
    rownums1 = [x for x in rownums if params[order[j]][x]==vd]
    
    k = 1
    vardivs2 = list(set(params[order[k]]))
    vardivs2.sort()
    for vd2 in vardivs2:
        rownums2 = [x for x in rownums1 if params[order[k]][x]==vd2]
        
        l = 2
        vardivs3 = list(set(params[order[l]]))
        vardivs3.sort()
        for vd3 in vardivs3:
            rownums3 = [x for x in rownums2 if params[order[l]][x]==vd3]
            
            roworder.extend(rownums3)


rr = rows[0:2]
rr.extend([rows[i+2] for i in roworder])
rows = rr

vardivs = [list(set(params[order[j]])) for j in range(len(order))]
[x.sort() for x in vardivs]

width = len(rows[0].split("&"))
rowrange = range(len(rows[2:]))
for r in rowrange:
    postpend = ""
    row = str(params[order[-1]][roworder[r]]) + " & " + rows[r+2]
    for j in range(len(order))[-2::-1]:
        divlength = reduce(operator.mul, [len(x) for x in vardivs[j+1:]])
        if r % divlength == 0:
            row = "\multirow{" + str(divlength) + "}{*}{" + str(params[order[j]][roworder[r]]) + "} & " + row
        elif r % divlength == divlength-1:
            row = " & " + row
            postpend = "\\cline{" + str(j+1) + "-" + str(width) + "}"
        else:
            row = " & " + row
    
    if r != rowrange[-1]: row = row[:-1] + postpend + row[-1:]
    rows[r+2] = row


out = open(outfile, "w")
for line in header:
    out.write(line)

out.write(multihead[0])

for line in rows:
    out.write(line)

for line in footer:
    out.write(line)

out.close()


