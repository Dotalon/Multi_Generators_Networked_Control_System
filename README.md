# SISOgna_NW
SISOgna il Network

We consider a Power Network System (PNS) composed of several power generation areas coupled
through tie-lines. Tie lines allow for power exchange between areas. In this way, if the load profile of
an area grows too much with respect to the power generated in that area, neighbouring areas can
contribute to satisfy the request.
We consider thermal power stations with single-stage turbines. In the figure below the typical
structure of a simple network of three areas is shown as an example along with measurable input
and output variables. The term 𝑃𝑖𝑗 is function of the coupling between areas i and j.

<img width="616" alt="Screenshot 2024-03-22 alle 10 39 39" src="https://github.com/Dotalon/SISOgna_NW/assets/117781101/c922455c-efa4-4a68-b34f-b351095e7af3">


Here we consider the dynamics of a PNS with five areas, where each area is endowed with primary
control. The models are linearized around an equilibrium value. The state of area i consists of the
vector 𝑥𝑖 = [Δ𝜃𝑖, Δ𝜔𝑖, Δ𝑃𝑚,𝑖, Δ𝑃𝑣,𝑖]^𝑇 (all variables are assumed measurable for simplicity), where:
 Δ𝜃𝑖: Deviation of the angular displacement of the rotor with respect to the stationary
reference axis on the stator.
 Δ𝜔𝑖: Speed deviation of rotating mass from nominal value.
 Δ𝑃𝑚,𝑖: Deviation of the mechanical power from nominal value (p.u.).
 Δ𝑃𝑣,𝑖: Deviation of the steam valve position from nominal value (p.u.).
The input 𝑢𝑖 of each area is the deviation of the power set-point from the nominal value. A system
composed of 5 areas, neglecting the disturbance term (i.e., the deviation of the requested load from
the nominal value), can be described by the following model:
𝑥̇ = 𝐴𝑥 + 𝐵𝑢
where 𝑥 = [𝑥1', … , 𝑥5']'
, 𝑢 = [𝑢1, … , 𝑢5]'
, and where the system matrices are specified in the file.



Problem:
1. Decompose the state and input vectors into subvectors, consistently with the physical
description of the system. Obtain the corresponding decomposed model.
2. Generate the system matrices (both continuous-time and discrete-time, the latter with a
sampling time ℎ = 1 s). Perform the following analysis:
a. Compute the eigenvalues and the spectral abscissa of the (continuous-time) system.
Is it open-loop asymptotically stable?
b. Compute the eigenvalues and the spectral radius of the (discrete-time) system. Is it
open-loop asymptotically stable?
3. For different state-feedback control structures (i.e., centralized, decentralized, and different
distributed schemes) perform the following actions
a. Compute the continuous-time fixed modes
b. Compute the discrete-time fixed modes
c. Compute, if possible, the CONTINUOUS-TIME control gains using LMIs to achieve the
desired performances. Apply, for better comparison, different criteria for computing
the control laws.
d. Compute, if possible, the DISCRETE-TIME control gains using LMIs to achieve the
desired performances. Apply, for better comparison, different criteria for computing
the control laws.
e. Analyze the properties of the so-obtained closed-loop systems (e.g., stability,
eigenvalues) and compute the closed-loop system trajectories (generated both in
continuous-time and in discrete-time) of the speed deviations Δ𝜔𝑖 starting from a
common random initial condition.
