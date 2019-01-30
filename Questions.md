# Questions that arose
A resource for finding the answers to the questions that stymied us (at least briefly).

## How to control the range of the simulation?
```r.simulate=(0,10,100)``` Simulates from time 0 to time 10 with 100 steps

```r.simulate(0,10,100,["substrate","rate law"])``` Still runs over the defined time frame, but now the results being produced are just (substrate,rate law), the results aren't being listed as controlled by time, but that is still the control.

Remember that this is fundamentally a kinetic package, so time will always be the controlling variable.

Handy trick: ```r.getSimulationData()``` will print the data coming out of the simulation and make it more obvious what is going on.

## Question title (use keywords so someone can search for this)
Question!

Answer! Once we have it.
