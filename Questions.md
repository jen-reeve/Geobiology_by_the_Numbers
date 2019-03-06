# Questions that arose
A resource for finding the answers to the questions that stymied us (at least briefly).

## How to control the range of the simulation?
```r.simulate(0,10,100)``` Simulates from time 0 to time 10 with 100 steps but these aren't evenly spaced steps.

To get what you would expect: ```r.simulate(start = 0, end = 10e5, steps = 100)```

```r.simulate(0,10,100,["substrate","rate law"])``` Still runs over the defined time frame, but now the results being produced are just (substrate,rate law), the results aren't being listed as controlled by time, but that is still the control.

Remember that this is fundamentally a kinetic package, so time will always be the controlling variable.

Handy trick: ```r.getSimulationData()``` will print the data coming out of the simulation and make it more obvious what is going on.

## Reversible and irreversible reactions
```->``` indicates irreversible

```=>``` indicates reversible

But Antimony doesn't actually care, this is just to help us keep track.

## Why can't I move cells around???
Ask the developer.

## What the frick are compartments and how can I actually use them?
Ask the developer.

## Can I tell the steady state solver to ignore a constantly accumulating species or to look for a steady state rate?
Can you pull rates out of an Antimony model? Yes!

This will print the steady state values for the things defined in the model as S1 and J12. Basically with the first command you are telling it to give you the steady state values for the things you request.
```
r.steadyStateSelections = ['[S1]','J12']
values = r.getSteadyStateValues()
selections = r.steadyStateSelections
print(values)
print(selections)
```
