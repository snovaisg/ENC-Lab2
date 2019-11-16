# read in data
x <- c(-1.16225338, -0.01236762, 0.23144468, -2.08263805, 2.64870304, -0.52868938,
       -1.52280636, 0.03085357, -1.54249244, -0.23164183, 2.23750583, -0.64678326,
       -0.82817693, -1.50508448, 0.87125402, -1.67430203, 1.63702338, -0.85729792,
       -0.67855079, -1.36315013, -1.70412209, -0.05967576, -0.88579241, -0.33469221,
       -0.74940615)

# (a) Use the observed data and the method of Newton-Raphson to approximate the ML 
# estimate of ??. Justify your choice of the initial estimate and all the intermediate steps.

# (b) Consider the following reparametrization b = log(??) in the pdf above. Implement the 
# Newton-Rapshon algorithm to approximate theML estimate of ?? using this reparametrization

# (c) Which approach, (a) or (b) is more sensitive to the initial values?