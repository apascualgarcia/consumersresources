import numpy as np
import re
# finds the asymptote value of an asymptotic function y, which is assumed to be constant
# or constantly oscillating for its last Npoints
def asymptote(y, Npoints):
    return np.mean(y[-Npoints:])
def remove_strings_from_file(matrices_folder,filename):
    file = open(filename + '.out', "r")
    metadata = []
    for x in file:
        name = x.replace(matrices_folder + '/RandTrix_Nr', '')
        name = name.replace('_Nc', ' ')
        name = name.replace('_Nest', ' ')
        name = name.replace('_Conn', ' ')
        name = name.replace('.txt', '')
        metadata.append(name)
    file.close()
    f = open(filename + '_filtered.out', 'w')
    for a in metadata:  # python will convert \n to os.linesep
        f.write(a + '\n')
    f.close()
    return

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def extract_numbers_from_string(filename):
    # printing original string
    # using re.findall()
    # getting numbers from string
    temp = re.findall(r"[-+]?\d*\.\d+|\d+", filename)
    res = list(map(float, temp))
    return res
