### DSO Market Clearing Formulation

This tool is companion software to the paper 'A DSO-Level Contract Market for Conditional Demand Response' submitted to the IEEE PES PowerTech Conference in Milan, 2019.
Using the parameters from the case study in the paper, determines the optimal combination of DR services, paying aggregators accordingly.

### INPUT
Sets - Contains the sets used in the DSO Market clearing tool

Params - Contains the parametes used in the DSO Market clearing tool

### PROCESSES
EcoGridDSO.gpr - The Project file for the Market Clearing Tool

EcoGrid2.gms - The "DSO Market clearing tool". GAMS file containing the set+parameter definition, formulation, and payment calculation

### OUTPUT
DSO_Market_Output.gdx - Output file storing sets, parameters, DR decisions, payments and profit outcomes
