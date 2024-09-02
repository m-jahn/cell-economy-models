import pandas as pd

# define class result where model results are collected
class result:
    def __init__(self, name, model, v, a, c, u):
        self.name = name
        self.model = model
        self.table = pd.DataFrame(model.load_results())
        self.v = v
        self.a = a
        self.c = c
        self.u = u
