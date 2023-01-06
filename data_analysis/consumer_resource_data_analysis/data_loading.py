import pandas as pd

def load_data_frame(file_location):
    col_names = ["G matrix", "A matrix", "G connectance", "G nestedness", "A connectance", "A nestedness", "A-mode","alpha0","number of sims","feasible volume", "dynamically stable volume", "dynamically unstable volume","marginally stable volume","av. dominant eigenvalue"]
    df = pd.read_csv(file_location,names=col_names,comment="#",sep=" ")
    return df
