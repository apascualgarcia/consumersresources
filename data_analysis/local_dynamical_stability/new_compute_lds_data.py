import consumer_resource_data_analysis as cf
import pandas as pd


data_location = "../results/all_data_NR25_NS25_verbose-level=1.out"

col_names = ["G matrix", "A matrix", "G connectance", "G nestedness", "A connectance", "A nestedness", "A-mode","alpha0","number of sims","percentage of feasible systems", "percentage of dynamically stable systems", "proportion of dynamically unstable systems","proportion of marginally stable systems","average largest eigenvalue"]
df = pd.read_csv(data_location,names=col_names,comment="#",sep=" ")

for alpha0 in cf.alpha0:
    print("alpha0 = ", alpha0)
    print(df[df["alpha0"]==alpha0][df["A-mode"]=='random_structure']['average largest eigenvalue'].to_list())
