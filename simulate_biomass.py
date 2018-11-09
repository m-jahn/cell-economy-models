# LIBRARIES ------------------------------------------------------------

from gekko import GEKKO


# determine biomass increase over time from previously executed 
# gekko models. Now directly integrated in other models.


b = GEKKO()

# time and growth rate mu are reused from previous gekko model
b.time = time
mu = b.Param(value=mu.value, name='mu')

# biomass accumulated over time
bm = b.Var(value=1, name='bm')

# biomass accumulation over time with initial biomass concentration bm
b.Equation(bm.dt() == mu*bm)

b.options.IMODE=4
b.solve()
