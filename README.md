# README
The scripts in this folder implement the preventive self-triggered control mechanism described in Gleizer and Mazo Jr. (2020). 

This work is part of the [SENTIENT project](https://mmazojr.3me.tudelft.nl/sentient/).

##  Requirements
Optimization Toolbox, Ellipsoidal Toolbox.

The Ellipsoidal Toolbox can be found [here](https://nl.mathworks.com/matlabcentral/fileexchange/21936-ellipsoidal-toolbox-et).

## How to reproduce the results of the paper

1. Install the Ellipsoidal Toolbox.
2. Replace the files `distance.m` and `intersection_ea.m` into `$MATLAB_DIR$/ellipsoids_112_lite/ellipsoids/@ellipsoid/`<sup>1</sup>.
3. Run `main.m`.

## How to test this STC with other linear systems

The script that simulates the STC is `mainnoise.m`. Its code is not very user-friendly, but to change the plant and controller, 
you should edit the block starting with `%% Load data` and, for defining the bounds on disturbances and noise, edit the block starting 
with `%% Assumption 4: Bound on disturbance`.

## References

[1] G. A. Gleizer and M. Mazo Jr. Self-triggered output-feedback control of LTI systems subject to disturbances and noise. Submitted to Automatica, 2020.

---
<sup>1</sup> There are compatibility issues between the Ellipsoidal Toolbox and recent versions of YALMIP. The scripts in this folder do not need the part of the Toolbox that uses YALMIP. However, some functions within the toolbox call YALMIP even without needing it. The provided replacements solve this issue. </small>
