# Overview

- the idea here is to be able to properly distinguish pigment concentrations from each other
- the flow cytometer excites pigments using one of two lasers (red or blue, channels suffixes -R or -B respectively) and detects flourescence due to that excitation with several channels (GRN-B, NIR-B, NIR_R, RED-B, RED-R, and YEL-B)
- these detection bands encompass a range of wavelengths which will include excitation from several pigments which "contaminates" the output a bit
  - i.e., measures of chlorophyll using the RED-B channel will include some amount of phycocyanin
- to this end, we will:
  1. integrate the total amount of absorption across all pigments for each channel
  2. integrate the amount of absorption for each pigment across each channel
  3. turn the into *relative* quantities of each pigment per channel
  4. scale these by the amount of excitation for each pigment per channel
  5. create a function that calculates, based on the input cytometer channels, the relative concentrations of the differnt pigments
- in theory we should arrive at some function like `pigments(GRN.B, NIR.B, NIR.R, RED.B, RED.R, YEL.B)` returning the expected relative concentrations of each pigment and the "total pigmentation" so that these relative quantities can be turned into absolute-ish quantities

## Example

Say that some hypothetical cytometer channel has a laser that excites pigments at 100 nm wavelength ($\lambda = 100$) and detects absorption between 200 and 300 nm (bandwidth $B = \{200, 300 \}$). With two pigments $A$ and $B$, with excitation values at 100 nm being $E_{A,100}(\lambda)$ and $E_{B,100}(\lambda)$ respectively and flourescence values, $F$, within the bandwidth given by $\int^{300}_{200} F_{A}(\lambda) \, \mathrm{d}\lambda$ and $\int^{300}_{200} F_{B}(\lambda) \, \mathrm{d}\lambda$ respectively.

From existing data, we already know (approximately) the values $E_A(\lambda)$ & $E_A(\lambda)$ and distributions $F_A(\lambda)$ & $F_B(\lambda)$ because other people have done assays on the absorption and flourescence of pigments.

The relative concentration of pigment $A$, $\rho_A$, should be computed as follows:

$$
\rho_A = \frac{E_{A,100}\displaystyle\int^{300}_{200} F_{A}(\lambda) \, \mathrm{d}\lambda}{\left(E_{A,100}(\lambda) \displaystyle\int^{300}_{200} F_{A}(\lambda) \, \mathrm{d}\lambda \right) + \left(E_{B,100}(\lambda)\displaystyle\int^{300}_{200} F_{B}(\lambda) \, \mathrm{d}\lambda \right)}
$$

The cytometer can provide us with a measurements of the total absorption only: $\mathrm{abs}_{tot} = \left(E_{A,100}(\lambda) \displaystyle\int^{300}_{200} F_{A}(\lambda) \, \mathrm{d}\lambda \right) + \left(E_{B,100}(\lambda)\displaystyle\int^{300}_{200} F_{B}(\lambda) \, \mathrm{d}\lambda \right)$. 

Therefore we calculate $\rho_A$ using the existing data and substitute in the value of $\mathrm{abs}_tot$ to get the relative concentration of pigment $A$.

---

*Note:* this makes the assumption that we measure all absorbed light equally across the bandwidth which I doubt is true but it's not the end of the world and is still better than just using the individual values that the band is centred around.
