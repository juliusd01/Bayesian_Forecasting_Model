import pyreadr

# Load the .RData file
result = pyreadr.read_r('Bayesian_Modelling/Homework/Homework#6/Test_2.RData')

# The result is a dictionary where keys are the names of the objects in the .RData file
print(result.keys()) # let's check what objects we got