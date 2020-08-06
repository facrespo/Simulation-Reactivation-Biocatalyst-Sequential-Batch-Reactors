# Simulation-of-the-Reactivation-of-Partially-Inactivated-Biocatalyst-in-Sequential-Batch-Reactors
Simulation of the Reactivation of Partially Inactivated Biocatalyst in Sequential Batch Reactors

This repository has the files to reply the simulation of reactors. The files using Python for simulate and R for analysis.

The reactors are define in the file reactorbioquim.py.

The function reaccionactivada make the case with reactivation, the function reaccion only make a simple reaction.

The example for reactivation is simulacion_experiment_with_activation.py.

The example for simple case is simulacion_experiment_without_activation.py.

The curve is calculated with simulation_curve_without_activation.py in the case without activation, y the curve is calculated with simulation_curve_with_activation.py in the case with activation.

The case with random parameters use the file: simulacion_reactores_reactivada_random.py, this file generate the file to make the analysis.

The result change the number of batches appear in result_example. 

The file of analysis is results_analysis.R.
