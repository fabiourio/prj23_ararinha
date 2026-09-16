'''
This script executes the Fokker 100 example
'''

#IMPORTS
from designTool.standard_airplane import standard_airplane
from designTool.analyze import analyze
from designTool.moment_of_inertia import moment_of_inertia

#=========================================

# SETUP

# Define airplane dictionary
airplane = standard_airplane('fokker100')

#=========================================

# EXECUTION

# Run standard test case and show results
analyze(airplane, plot=True, print_log=True)

# Execute the moment of inertia calculation
# for given factors of fuel and payload weight.
moment_of_inertia(airplane, fuel_frac=0.75, payload_frac=1.0)

# Print results
print('Computed Moments of Inertia')
print(airplane['moment_of_inertia'])