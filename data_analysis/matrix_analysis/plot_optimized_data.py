import numpy as np
import pandas as pd
import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt

data_file_path = "data_output/optimized_matrices_characteristics.csv"


data = pd.read_csv(data_file_path)
for connG in [0.13]:
    to_plot = data[data['connG']==connG]
    to_plot.plot(x = 'connG', y = 'connA')
