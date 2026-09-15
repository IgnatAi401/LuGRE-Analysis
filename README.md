# LuGRE

This repository contains **MATLAB analysis scripts** for the [LuGRE](https://etd.gsfc.nasa.gov/our-work/lunar-gnss-receiver-experiment-lugre/) (Lunar GNSS Receiver Experiment) [Mission Data](https://zenodo.org/records/16411687).

## Overview

- **Mission:** LuGRE — Lunar GNSS Receiver Experiment
- **Data:** Public LuGRE/BGM1 mission data, including raw GPS/Galileo L1/E1 and L5/E5 observables and signal sample batches
- **Language:** MATLAB
- **Purpose:** GNSS-based positioning, navigation, and orbit determination in cislunar and lunar scenarios

## Data Sources and Ephemerides

The analysis uses:

- **NASA hybrid broadcast ephemeris** for GNSS satellite orbits
- **JPL Horizons lunar ephemeris** for lunar position and geometry

These inputs support the reconstruction and analysis of LuGRE observations in the Moon-transfer and lunar environment.

## Algorithms and Methods

The implemented workflow includes:

- **Observation-level Direct Position Estimation (DPE)**  
  Direct estimation of position from GNSS observables at the observation level.

- **Least-squares estimation**  
  Used for parameter estimation and navigation solution computation.

- **Restricted three-body orbital dynamics optimization**  
  Used for improved orbit determination and dynamical consistency in the Earth–Moon system.

## Key Highlights

- Integrates **NASA hybrid broadcast ephemeris** and **JPL Horizons lunar ephemeris**
- Combines **observation-level DPE**, **least-squares estimation**, and **restricted three-body orbital dynamics optimization**
- Supports **GNSS-based cislunar/lunar navigation and orbit-determination analysis**
- Built on **public LuGRE mission data** with **MATLAB** analysis scripts