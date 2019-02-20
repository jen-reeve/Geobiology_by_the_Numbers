# Class notes
## Jan 16 2019
using [Tellurium](tellurium.analogmachine.org) which is built for biological networks (aka metabolic networks)

If you want get Systems Biology: Introduction to Pathway Modeling (on the [analogmachine website](http://books.analogmachine.org/))

Boz wants to be able to go to sites like [Biomodels](https://www.ebi.ac.uk/biomodels/) (where can browse 10k models) and mess around. Tellurium makes it very easy to import models from Biomodels.

Today let's install and play with some simple models.

Boz's goals for the class:
* get familiar with kinetic modeling platform
* use to reproduce some figures from paper from Halevy and Wing 2014 (sulfate reduction modeling)
* to use platform to calculate isotopic fraction in networks
* to enter a proposed early network and then test it
* another early metabolic network model

### concentration of a single molecule in an E. coli cell
Guess radius to be 1um, so the volume (assuming a cube) is 1um^3. There are 1000L in 1m^3 and 1mL is 1cm3. The volume of the cell in L is 1x10^(-15)L

1 molecule is 1x10^(-24)moles

So the concentration of 1 molecule in an E. coli cell is 1x10^(-9)M aka 1nM so this is the lowest possible concentration of a metabolite/substance of interest in a cell. And while we are in the nM range we can't have halves.

The total concentration of stuff in an E. coli cell is ~100mM.

For the future, use [Bionumbers](https://bionumbers.hms.harvard.edu/search.aspx) or [Cell Biology by the Numbers](book.bionumbers.org)

### Running code?!
Successfully ran the simple and slightly more complex examples in the Tellurium Notebook interface

We also tested defining models using the usage example method of defining the model. This also introduces the ability to define compartments.

Book for Earth system model: Mathematical modeling of Earth's Dynamic Systems

For next week, think about explicit examples that you would like to reproduce using this platform (id papers or models). If you have time play with it a little bit. Come in with one new thing that you have discovered. And buy the Systems Bio book!

## Jan 23 2019
Things to not do: biochemistry

Try to find pathways where the steps are known and the kinetics are constrained.

Maybe present models as final SuperLab meeting.

Build a simple graphical representation of the model we want to tackle. We don't want to create new biology, want to take a conceptual model from the literature and build a Tellurium/Antimony model of it.

Looking at the Log Plotting (Calzone) example notebook.

What do the compartment units mean? Probably volume units to change the concentration values into total units.

Often good to start at the end of your metabolic network as the final steps are often considered to be irreversible. A useful fiction.

Quick derivation of this type of transformation that is enzymatically catalyzed (from Noor et. al. 2013 paper):

"Quantum" (unit) process of biochemistry can be represented graphically: Substrate + Enzyme <=> Enzyme-Substrate (Complex)=> Product + Enzyme

This is Michaelis-Menten model (but we will call it the Menten model because she deserves all the credit).

Three kinetic parameters for the model: k_(+1) (forward direction of first step), k_(-1) (reverse direction of the first step), k_(+2) (forward direction of the second step).

We know from elementary chemical theory that the rate of this reaction is equal to the rate of accumulation of the product.

dP/dt = V(S) (where V(S) is the reaction velocity which is dependent on the concentration of the substrate)

This doesn't depend on the concentration of the product, but definitely depends on the concentration of substrate

V(S) = k_(+2) * C (the complex)

Menten assumed that the concentration of the substrate was much much greater than the concentration of the enzyme. so: d(ES)/dt = dC/dt = 0

You end up with a pair of rate laws that reflect these transformations.

1. (k_(+2)+k_(-1)) C = k_(+1) S * E (consumption of the enzyme-substrate complex is equal to the production of the enzyme-substrate complex).

   These rate constants have different units. Consumption has units per time. Production has units of per time*something

2. E_(total) = E + C

   The total enzyme present is equal to the free enzyme times the amount bound up in a complex

E = E_(total) - C

\[
(k_(+2)+k_(-1)) C + k_(+1) S * C = k_(+1) E_(total) * S
\]

$$C = (k_(+1) E_(total) * S) / ((k_(+2)+k_(-1)) + k_(+1) * S) = (E_(total) * S) / ((k_(+2) + k_(-1))/k_(+1) + S)$$

Endmember 1: If the concentration of the substrate is much much bigger than the kinetic ratio (k_(+2)+k_(-1))/k_(+1) then V = k_(+2) * E_(total) which defines V_(max) = k_(cat) * E_(total)

Endmember 2: if the concentration of the substrate is equal to the kinetic ratio (see above) then V = (k_(+2) * E_(total) * S) / (2 * S) = (k_(+2) * E_(total)) / 2 = V_(max) / 2

Now getting back to our lactose model:

Use [Brenda](https://www.brenda-enzymes.org/) to get the enzyme kinetic information.

## Feb 6 2019
Set up:

Diffusion across the membrane: 0 = area * permeability * (S1-S2)

Consumption of the internal substance, with a silent 3rd substance implied: 0 = (Vmax * S2)/(km + S2)

Steady state: dS/dt = 0

At equilibrium all of the net rates are 0, at steady state all the net rates are constant and all the same.

What here defines our net rate? S1 is constant, we start with no S2 and no S3. Once we start we begin to accumulate S2 and S3 and after a while the rate of accumulation of S3 is linear. This linear accumulation rate is the net rate of the system.

You can't have a steady state without the system having "openness" b/c you are accumulating your final product constantly. Defining your first substrate as constant implies/forces an open system. This is a nuance we should consider.

Perhaps clearer to initially define our equations as:

J_ss = area * permeability * (S1 - S2)

J_ss = (Vmax * S2)/(km + S2)

Final form for quadratic: 0 = A*P*S2^2 + S2 * (Km*A*P - A*P*S1 + Vmax) - Km*A*P*S1

a = A*P

b = km*A*P - A*P*S1 + Vmax

c = -km*A*P*S1

### Units

S1/S2 has units of mass/length^3

km has units of mass/length^3

Area has units of length^2

permeability has units of length/time

so Vmax has units of mass/time

Remember that we determined that Vmax = kcat * Etotal

kcat had units of turnover time aka 1/time but really it was mass of substrate/(mass of enzyme * time)

So Vmax has units of mass substrate per time per length^3

So now we have incorporated our volume into the Vmax term, this is the total Vmax, aka on a per cell basis. So to compare these two Vmax terms (the Vmax per enzyme and the Vmax per cell): Vmax_cell = constant * Vmax_enyzme where the constant must have units of length^3 so it must be the volume of the cell.

Vmax_cell is essentially the amount of reactions in the cell

Use our minimum and maximum concentrations in cells as a way to detect units problems.

Our concentrations are given in mass per length cubed given as molarity which is moles per liter. This means that the length unit we are using is the cube root of a liter. Which is (10^3 cm^3)^(1/3) or 10 cm or 1 dm.

We have taken our cell, called it a sphere and then unfurled it and put walls at each end of the membrane so that we can consider things in terms of the surface area. This implies that there is not differential spatial distribution of things in the cells.

For next class:

1. take your antimony code and compare it at steady state to the values you just got. Going to need to know some commands about steady state values (already in code stolen from Boz).
2. explore how the steady state values relate to the external lactose concentration. Make plots of internal vs. external lactose concentration and use the python code.

### Calculations
Calculating what we expect the steady state concentration of S3 (intracellular lactose) to be:
```
import math
import numpy as np
import tellurium as te
# this is our calculation to check our numerical calculations
# 1 liter = 1 dm^3
area = 1.0e-10; # dm^2
permeability = 1.0e-8; # dm s-1
km = 5.5e-3; # M
S1 = 1; # M
kcat = 6.42e2; # s-1
etot = 1.0e-6; # M
vol_cell = 10e-15; # L or dm^3
a = area * permeability;
b = km*area*permeability-area*permeability*S1+kcat*etot*vol_cell;
c = -1*km*area*permeability*S1;
S3_ss=(-b+(b**2-4*a*c)**0.5)/(2*a);
print(a,b,c)
print('S3_ss',S3_ss)
```
Gives 86uM which is nicely in between our limits of physiological concentrations.

Checking that a steady state concentration gives net rates that are equivalent for our two J_ss equations:
```
J_ss1=area*permeability*(S1-S3_ss);
print('J_ss1',J_ss1)
J_ss2=(kcat*etot*vol_cell*S3_ss)/(km+S3_ss);
print('J_ss2',J_ss2)
```
Woo these give the same number. 10^-21 moles per second (per cell) which converts to 10^-16 moles per day (per cell) aka one femtomole per day (per cell) we will call this unit a Jorgensen aka a Jorg. This is a reality check, we put in numbers out of nowhere and got a number that was hard to understand, so we can convert into Jorgs and that is a reasonable process rate for cells in the environment. In culture these rates will always be faster.

Model:
```
model LacPath
    compartment cytoplasm;
    cytoplasm = 1.0;
    compartment environment;
    environment = 1000.0;
    species lactose, glucgalac;
    const species ex_lactose;
    unit M = mole / liter
    unit inv_sec = 1. / seconds
    unit um = 10.e-6 meters
    unit vol_conv = 1000 liter/1 meter^3

    lactose = 0.0 M; # M; free parameter
    glucgalac = 0.0 M;
    ex_lactose = 0.1 M;

    lactose in cytoplasm; glucgalac in cytoplasm;
    ex_lactose in environment

    J34: lactose -> glucgalac; kcat*Etot*lactose/(Km+lactose);

    kcat = 6.42e2 inv_sec; #s^-1; from Juers et al (2012)
    Etot = 1.e-6 M; # M; free parameter
    Km = 5.5e-3 M; # M; from BRENDA

    #J13: ex_lactose => lactose; area*permeability*(ex_lactose-lactose)*conversion;
    #J31: lactose => ex_lactose; area*permeability*(lactose-ex_lactose)*conversion;

    area = 5 um^2;
    permeability = 3 um * 1 inv_sec;
    conversion = 1 vol_conv;
end
```

Running the model
```
LacPath.reset()
LacPath.simulate(0., 100., 1000) #'S1','J1'
LacPath.plot()

values = LacPath.getSteadyStateValues()
selections = LacPath.steadyStateSelections
print(values)
print(selections)
```

## Homework from Feb 5 2019
not functional code
```
import math
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
# this is our calculation to check our numerical calculations
# 1 liter = 1 dm^3
area = 1.0e-10; # dm^2
permeability = 1.0e-8; # dm s-1
km = 5.5e-3; # M
S1 = np.linspace(0.1,10,100); # M
kcat = 6.42e2; # s-1
etot = 1.0e-6; # M
vol_cell = 10e-15; # L or dm^3
a = area * permeability;
b = km*area*permeability-area*permeability*S1+kcat*etot*vol_cell;
c = -1*km*area*permeability*S1;
S3_ss=(-b+(b**2-4*a*c)**0.5)/(2*a);
#print(a,b,c)
#print('S3_ss',S3_ss)
plt(S1,S3_ss)


```

## Feb 13 2019 class
Diffusion:

Looking at the movement of lactose from the environment (external lactose)
 into the cell (internal lactose). There is resistance to diffusion given
  by the permeability multiplied by the surface area. We end up with a
  diffusive flux that is equal to the permeability times the surface area
   times the difference in concentrations between the external and
   internal concentrations. This is a net flux. If you wanted to
   calculate the flux in each direction separately then the flux in would
    be the permeability times the area times the external concentration
    and the flux out would be the permeability times the area times the
    internal concentration.

In this case the net flux is easy and equivalent. If we get into
isotopic stuff then we need to look at the ratio of the flux in
each direction.

Next steps: do we want to work on changing our units into biochemists'
units or add another step? By putting things in biochemists' units
(per mg cell protein) we might be able to make adding the periplasm a
 lot easier.

The periplasm is approximately 10nm thick. So the relative volume of the
 periplasm is 10^-6 times the volume of the cell aka a million times smaller.
 This would make the concentration of one molecule in the periplasm
 approximately 1mM. And the maximum concentration is still 100mM, so
 we have a very narrow range of concentrations possible in the periplasm.

For next week: add a second rate law that we are just considering passive
transport across the inner membrane. And then next week we can translate
this into an active transport law. And bring Boz something that he doesn't
know about Antimony.

### Code
```
import math
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
# this is our calculation to check our numerical calculations
# 1 liter = 1 dm^3
area = 1.0e-10; # dm^2
permeability = 1.0e-8; # dm s-1
km = 5.5e-3; # M
S1 = np.linspace(0.001,1,100); # M
kcat = 6.42e2; # s-1
etot = 1.0e-6; # M
vol_cell = 1.0e-15; # L or dm^3
a = area * permeability;
b = km*area*permeability-area*permeability*S1+kcat*etot*vol_cell;
c = -1*km*area*permeability*S1;
S3_ss=(-b+(b**2-4*a*c)**0.5)/(2*a);
#print(a,b,c)
#print('S3_ss',S3_ss)
plt.plot(S1,S3_ss)
plt.ylabel('Internal lactose (M)')
plt.xlabel('External lactose (M)')
```

```
J_ss1=area*permeability*(S1-S3_ss)*(60*60*24)*10e15;#fmol/day
#print('J_ss1',J_ss1)
J_ss2=(kcat*etot*vol_cell*S3_ss)/(km+S3_ss)*(60*60*24)*10e15; #fmol/day
#print('J_ss2',J_ss2)
plt.plot(S1,J_ss1)
plt.xlabel('External lactose (M)')
plt.ylabel('Diffusion (fmol/day)')
```

```
plt.plot(S1,J_ss2)
plt.xlabel('External lactose (M)')
plt.ylabel('Lactose conversion (fmol/day)')
```

Model:
```
model LacPath
    compartment cytoplasm;
    cytoplasm = 1.0e-6;
    compartment environment;
    environment = 1.0;
    species lactose, glucgalac;
    const species ex_lactose;
    unit M = mole / liter
    unit inv_sec = 1. / seconds
    unit um = 10.e-6 meters
    unit vol_conv = 1000 liter/1 meter^3
    unit dm = 0.01 meters

    lactose = 0.0; # M; free parameter
    glucgalac = 0.0;
    ex_lactose = 0.1;

    lactose in cytoplasm; glucgalac in cytoplasm;
    ex_lactose in environment

    J34: lactose -> ; kcat*Etot*vol_cell*lactose/(Km+lactose);

    kcat = 6.42e2; #s^-1; from Juers et al (2012)
    Etot = 1e-6; # M; free parameter
    Km = 5.5e-3; # M; from BRENDA
    vol_cell = 1.0e-15; # L or dm^3

    J13: ex_lactose => lactose; surf_area*permeability*(ex_lactose-lactose);
    #J31: lactose => ex_lactose; surf_area*permeability*(lactose-ex_lactose)*conversion;

    surf_area = 1.0e-10; # dm^2
    permeability = 1.0e-8; # dm s-1
    #conversion = 1 vol_conv;
end
```

```
LacPath.reset()
LacPath.simulate(0., 10e10, 1000,['time','[lactose]']) #'S1','J1'
LacPath.plot()

LacPath.steadyStateSelections = ['[lactose]','J34','J13']
values = LacPath.getSteadyStateValues()
selections = LacPath.steadyStateSelections
print(values)
print(selections)
```

## Homework for Feb 20th 2019
Something new about Antimony: the ability to have events (where passing a threshold changes other things)
and to have DNA strands so that production of enzymes is tied to a regulatory event.

Plotting of steady state rates against external lactose concentration is in
the previous code.

Model edits:
```
model LacPath
    compartment cytoplasm;
    cytoplasm = 1.0e-6;
    compartment environment;
    environment = 1.0;
    compartment periplasm;
    periplasm = 1.0e-6;
    species lactose, glucgalac, peri_lactose;
    const species ex_lactose;
    unit M = mole / liter
    unit inv_sec = 1. / seconds
    unit um = 10.e-6 meters
    unit vol_conv = 1000 liter/1 meter^3
    unit dm = 0.01 meters
    #pi = 3.1415;

    lactose = 0.0; # M; free parameter
    glucgalac = 0.0;
    peri_lactose = 0.0;
    ex_lactose = 0.1;

    lactose in cytoplasm; glucgalac in cytoplasm;
    peri_lactose in periplasm;
    ex_lactose in environment

    J34: lactose -> ; kcat*Etot*vol_cell*lactose/(Km+lactose);

    kcat = 6.42e2; #s^-1; from Juers et al (2012)
    Etot = 1e-6; # M; free parameter
    Km = 5.5e-3; # M; from BRENDA

    J12: ex_lactose => peri_lactose; surf_area_peri*permeability*(ex_lactose-peri_lactose);
    J23: peri_lactose => lactose; surf_area*permeability*(peri_lactose-lactose);
    #J31: lactose => ex_lactose; surf_area*permeability*(lactose-ex_lactose)*conversion;

    # cell area/volume calculations

    cell_radius = 400.0e-9 # E coli radius from Bionumbers
    vol_cell = (4/3)*pi*cell_radius^3
    surf_area = 4*pi*cell_radius^2
    #surf_area = 1.0e-10; # dm^2
    peri_thickness = 10.0e-9 # 10nm
    cell_peri_radius = cell_radius+peri_thickness
    vol_peri = ((4/3)*pi*cell_peri_radius^3)-vol_cell
    surf_area_peri = 4*pi*cell_peri_radius^2

    # membrane properties
    permeability = 1.0e-8; # dm s-1

end
```

```
LacPath.reset()
LacPath.simulate(0., 10e14, 1000,['time','[lactose]','[peri_lactose]']) #'S1','J1'
LacPath.plot()

LacPath.steadyStateSelections = ['[lactose]','J34','J12','J23','[peri_lactose]']
values = LacPath.getSteadyStateValues()
selections = LacPath.steadyStateSelections
print(values)
print(selections)
```

## February 20th 2019
### New things:
* don't close the Tellurium notebook until the "save successful" message appears, double check the "last save" time in the lower right hand corner to make sure the save happened
* DNA! You can pass rates onto downstream elements and use promoters to control gene product production
* Events: can have threshold values that lead to other changes in the system (changing parameters)
* Units: it is mostly just for record keeping, but it can be very helpful
* Arrays: can make an array ```np.arange(start,stop,step size)``` or ```np.linspace(start,stop,step number)``` in a python cell with numpy (hence np.)
* Subplots: load in at the top ```te.newTiledFigure(rows=2,cols=3)```, and then define each plot
* Use Shift+Enter to run scripts without using the button
* membrane transporter protein: e is the number of moles of pore per area (# of proteins/membrane area)
* transport saturation: J_a=(e x kf x S1 - kr x S2)/(1 + (S1 / km1) + (S2 / km2))
* when simulaing the tellurium model the third option isn't actually steps divided evenly, but if you explicitly define these terms then it will do what you think it should do: ```LacPath.simulate(start = 0, end =5.64e4, steps = 100)```

Plotting:

Can use matplotlib ```import matplotlib.pyplot as plt``` and then use commands like ```plt.plot(x,y)``` or  ```plt.scatter(x,y)```

Tellurium includes the plotly engine and as long as you have imported Tellurium you can use ```te.plot(x,y)```

Even better is that we can change our axes to logarithmic via ```te.plot(x,y,logx=True,logy=False)```

Time of culturing is on the order of hours to days. Growth stops because you run out of something essential.

We have an idea about the steady state rate of a cell and can hopefully scale this to a culture and then calculate how long it would take for the culture to run out of lactose.

An OD of 1.0 is approximately 10^9 or 10^10 cells per mL.

But we are in the GEOLOGY DEPARTMENT we don't care about cell time scales, we care about radioactive decay time scales. The rate law that defines our world is: dNp/dt = - lambda Np, Np = Np(0)e^(-lambda t). Lambda is a rate constant with units per time.

Np/Np(0) = 0.5 = e^(-lambda t1/2) so ln(0.5) = -lambda t1/2 so t1/2 = ln(0.5)/-lambda

We can calculate a "half life" of a culture. To use half the stuff in a culture takes about a day. This gives us a way to check if the rate we are getting is reasonable given the real world.

We need to find our rate constant, what are the units? J_ss is in moles per cell per time. How do we convert? We need to know something about moles and something about cells. Something about moles: S1 = 10^-3 moles per liter so multiply by the volume of the culture (in L). Something about cells: 10^9 cells per mL so multiply by the volume of the culture (in mL) to get total number of cells.

How are we going to take our steady state rate and convert it into a decay constant? Jss * (1/S1) * (1/volume) * cell numbers

The actual half life should be longer because the rate decreases as the nutrient concentration decreases. And we have assumed constant cell numbers which has issues... b/c we used 10^9 which is essentially a final cell density, the half life should be longer yet again.

French press detour.

What do we know? We know that our steady state rate gives us a time constant that is sensible.

How can we use this to evaluate our Tellurium model? Take the steady state rate that we calculate from the model and calculate a half life that is on the order of hours to days. Basically do the same conversions as with the analytical model.

### for next time
go try to get the things that boz has working working for you too!
* figure out compartments (Boz did this by messing with radioactive decay)
  * compartments seem to modify the rate constant
  * change the rate law to include the compartment size to avoid the rate constant depending on the compartment size
  * redefine our compartments in terms of things we actually care about (culture volume, cell volumes etc)



### Code
```
model LacPath
    compartment cytoplasm;
    cytoplasm = 1.0e-6;
    compartment environment;
    environment = 1.0;
    compartment periplasm;
    periplasm = 1.0e-6;
    species lactose, glucgalac, peri_lactose;
    const species ex_lactose;
    unit M = mole / liter
    unit inv_sec = 1. / seconds
    unit um = 10.e-6 meters
    unit vol_conv = 1000 liter/1 meter^3
    unit dm = 0.01 meters
    #pi = 3.1415;

    lactose = 0.0; # M; free parameter
    glucgalac = 0.0;
    peri_lactose = 0.0;
    ex_lactose = 0.1;

    lactose in cytoplasm; glucgalac in cytoplasm;
    peri_lactose in periplasm;
    ex_lactose in environment

    J34: lactose -> ; kcat*Etot*vol_cell*lactose/(Km+lactose);

    kcat = 6.42e2; #s^-1; from Juers et al (2012)
    Etot = 1e-6; # M; free parameter
    Km = 5.5e-3; # M; from BRENDA

    J12: ex_lactose => peri_lactose; surf_area_peri*permeability*(ex_lactose-peri_lactose);
    J23: peri_lactose => lactose; surf_area*permeability*(peri_lactose-lactose);
    #J31: lactose => ex_lactose; surf_area*permeability*(lactose-ex_lactose)*conversion;

    # cell area/volume calculations

    cell_radius = 400.0e-9 # E coli radius from Bionumbers
    vol_cell = (4/3)*pi*cell_radius^3
    surf_area = 4*pi*cell_radius^2
    #surf_area = 1.0e-10; # dm^2
    peri_thickness = 10.0e-9 # 10nm
    cell_peri_radius = cell_radius+peri_thickness
    vol_peri = ((4/3)*pi*cell_peri_radius^3)-vol_cell
    surf_area_peri = 4*pi*cell_peri_radius^2

    # membrane properties
    permeability = 1.0e-8; # dm s-1

end
```

```
LacPath.reset()
LacPath.simulate(0., 10e14, 1000,['time','[lactose]','[peri_lactose]']) #'S1','J1'
LacPath.plot()

LacPath.steadyStateSelections = ['[lactose]','J34','J12','J23','[peri_lactose]']
values = LacPath.getSteadyStateValues()
selections = LacPath.steadyStateSelections
print(values)
print(selections)
```

```
area = 1.0e-10; # dm^2
permeability = 1.0e-8; # dm s-1
km = 5.5e-3; # M
S1 = np.linspace(0.001,1,100); # M
kcat = 6.42e2; # s-1
etot = 1.0e-6; # M
vol_cell = 1.0e-15; # L or dm^3
a = area * permeability;
b = km*area*permeability-area*permeability*S1+kcat*etot*vol_cell;
c = -1*km*area*permeability*S1;
S3_ss=(-b+(b**2-4*a*c)**0.5)/(2*a);
#print(a,b,c)
#print('S3_ss',S3_ss)
#plt.plot(S1,S3_ss)
#plt.ylabel('Internal lactose (M)')
#plt.xlabel('External lactose (M)')
te.plot(S1,S3_ss,logx=True,logy=False) # te.plot seems to ignore requests to set logy=False

# or try this by taking the log of the values and then plot
logS1 = np.log10(S1) # log() is the natural log!
logS3ss = np.log10(S3_ss)

te.plot(logS1,logS3ss)
```

```
J_ss1=area*permeability*(S1-S3_ss)
culture_volume_mL = 10 # mL
culture_volume_L = culture_volume_mL * 1e-3
cell_density = 1e9 # cells per mL
cell_number = cell_density * culture_volume_mL
k_ss = J_ss1 * (1/S1) * (1/culture_volume_L) * cell_number
half_life_ss = -1 * np.log(0.5) / k_ss
half_life_days = half_life_ss /(60*60*24)
te.plot(S1,half_life_ss)
te.plot(S1,half_life_days)
```
