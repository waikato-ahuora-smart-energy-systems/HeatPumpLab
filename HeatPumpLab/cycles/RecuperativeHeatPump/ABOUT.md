# Recuperative heat pump model for high temperature glide heating

Written by Tim Walmsley and Keri Adamson, University of Waikato, 2021-2025 (tim.walmsley@waikato.ac.nz)

## Introduction
The idea of heat recuperation within vapour-compression and transcritical heat pump cycles is not new. Internal heat exchangers (IHXs) are widely used to reduce flash gas generation after expansion, preheat suction gas, provide greater superheat after compression, and enable slightly lower high-side pressures (e.g. when producing hot water). In conventional cycles, the IHX transfers heat from the warm refrigerant immediately before expansion to the evaporator pressure. Pinch analysis states that heat is in excess below the pinch temperature and in deficit above it. A heat pump should ideally absorb heat below the pinch and reject upgraded heat above it. Applying this logic to the heat pump cycle itself suggests that repositioned IHXs could deliver thermodynamic benefits.

This model calculates heat pump perfromance for alternate IHX placements. To distinguish a conventional IHX from a repositioned one, this study uses the term recuperative heat exchanger (RHX). Below is a comparsion of four cycles that can be analysed using this model. 

<img width="452" height="514" alt="image" src="https://github.com/user-attachments/assets/cff96d29-63b3-4fed-af2a-3db739343726" />

### Standard heat pump (HP)
The baseline cycle is a conventional single-stage vapour compression heat pump. It comprises a compressor, a gas cooler, an expansion valve, internal heat exchanger, and an evaporator. No extra internal recuperation is included. This configuration serves as the performance reference for evaluating the benefits of the recuperative modifications.

### Recuperative heat pump 1 (RecuperativeHP1)
RecuperativeHP1 introduces a single internal recuperator between the outlet of the gas cooler and the compressor suction line. In this configuration, partially cooled high-pressure refrigerant transfers heat to the low-pressure suction stream before compression. The aim is to superheat the suction gas, thereby reducing the maximum pressure required while maintaining similar COP. This cycle represents the simplest form of recuperative modification.

### Recuperative heat pump 2 (RecuperativeHP2)
RecuperativeHP2 extends the concept by placing the source side of the recuperator between two gas coolers. This positioning allows the recuperator to operate at a higher temperature level than in RecuperativeHP1, increasing the degree of suction gas superheating. The expectation is further reduction of high-side pressure, though the COP trade-off is less predictable.

### Recuperative heat pump 3 (RecuperativeHP3)
RecuperativeHP3 has two compression stages with two recuperators. It combines elements of RecuperativeHP1 and RecuperativeHP2. This configuration may unlock higher suction superheat entering the second compression stage while respecting the pinch point(s) that arise. While the potential reduction in maximum pressure may be greater, this cycle is mechanically more complex and requires greater capital investment.

Contact: Tim Walmsley, tim.walmsley@waikato.ac.nz
