import consumer_resource_data_analysis as cf

import re

#### permet de transformer une log file en un txt lisible par numpy

log_file = "/Users/Shared/Master/Master_Thesis/results/optimize_matrices_core_0.log"
new_text = ""
with open(log_file) as file:
    i=0

    for line in file:
        if i>5:
            to_add = (re.sub("[^0-9.]"," ", line[51:len(line):1]))
            to_add = re.sub(' +', ';', to_add)
            new_text+= (to_add[0:len(to_add)-1:1] + "\n")
        i+=1

fout = open(log_file[:len(log_file)-3]+"txt","w")
fout.write(new_text)
fout.close
