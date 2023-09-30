# Overview

- the idea here is to be able to properly distinguish pigment concentrations from each other
- the flow cytometer excites pigments using one of two lasers (red or blue, channels suffixes -R or -B respectively) and detects flourescence due to that excitation with several channels (GRN-B, NIR-B, NIR_R, RED-B, RED-R, and YEL-B)
- these detection bands encompass a range of wavelengths which will include excitation from several pigments which "contaminates" the output a bit
  - i.e., measures of chlorophyll using the RED-B channel will include some amount of phycocyanin
- to this end, we will:
  1. integrate the amount of absorption for each pigment across each channel
  2. turn the into *relative* quantities of each pigment per channel using the total amount of absorption across all pigments for each channel
  3. scale these by the amount of excitation for each pigment per channel
  4. calculate the weighted mean of the pigment using the relative excitation-flourescence as weights
- in theory we should arrive at some function like `pigments(GRN.B, NIR.B, NIR.R, RED.B, RED.R, YEL.B)` returning the expected relative concentrations of each pigment and the "total pigmentation" so that these relative quantities can be turned into absolute-ish quantities for but now its just a csv and you can use the `weighted.mean` function

*Note:* this makes the assumption that we measure all absorbed light equally across the bandwidth which I doubt is true but it's not the end of the world and is still better than just using the individual values that the band is centred around.
