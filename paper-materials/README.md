To reproduce the environment within which these notebooks were ran, type: 
``` 
conda env create --file=pb-environment.yaml 
```

This directority contains the code for the main text figures:
- M1-sub-domain-map.ipynb --- map of the domain of the Pb model
- M2-boundary-conditions.ipynb --- temperature-salinity diagrams and transects with observed and estimated dissolved Pb concentrations used for the boundary conditions in the dissolved Pb model
- R3-evaluation-fields.ipynb --- depth slices of simulated dissolved Pb concentrations alongside observations from GEOTRACES cruises
- R4-evaluation-transect.ipynb  --- a transect connecting GEOTRACES cruise sample stations, with observed and simulated dissolved Pb concentrations
- R5-residence-time.ipynb --- comparison of dissolved Pb residence time in the model and from observational studies in Arctic and sub-Arctic oceans
- R6-budget.ipynb --- calculations and figure of the Pb budget estimated from these simulations
- R7-atlantic.ipynb --- depth slices and transect of regions most impacted by Atlantic water using dissolved Pb as a tracer
- R8-isolation.ipynb --- similar to above but looking at how isolated Baffin Bay is from contributions from the shelves / Atlantic water
- R9,S10,S11-jet.ipynb --- visualization of a seasonal jet of Atlantic water that crosses through Davis Strait as identified with dissolved Pb, looking at different depth levels and months

And, supplementary material figures:
- S1-boundary-condition-method.ipynb --- explanation of the method used to develop the model dissolved Pb boundary conditions
- S2-initial-conditions.ipynb --- model dissolved Pb initial conditions
- S3-monthly-boundary-condition-transect.ipynb --- transects of monthly dissolved Pb boundary conditions in the Labrador Sea
- S4-hudson-bay-boundary-sensitivity.ipynb --- visualization of the impact of the Hudson Bay boundary on dissolved Pb concentrations in Baffin Bay
- S5-Pb-black-carbon-ratio.ipynb --- derivation of the ratio of Pb to black carbon
- S6-sediment-resuspension.ipynb --- visualization of the sediment resuspension forcing field
- S7-sediment-in-ice.ipynb --- visualization of the sediment in sea ice forcing field
- S8-spin-up.ipynb --- assessment of spin up of model experiments based on month-to-month change in Pb concentrations at evaluation stations
- S9-region-definitions.ipynb --- regions as defined for the residence time evaluation
