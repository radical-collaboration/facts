# FittedISMIP/icesheet

This module implemented the parametric fit to ISMIP6 Greenland Ice Sheet projections described in IPCC AR6 WG1 9.SM.4.4. As described therein:

A polynomial fit to the ISMIP6 results is employed to calculate rates of change. The parametric fit is a cubic fit to temperature and a quadratic fit over time:

$$ 𝜕𝑠/𝜕𝑡 =𝛽_0 +𝛽_1𝑇+𝛽_2𝑇^2+𝛽_3𝑇^3+𝛽_4𝑡+𝛽_5𝑡^2 $$

Where $s$ indicates the sea-level equivalent contribution in mm, $T$ is GSAT in °C, and $t$ is time in years. For the purposes of fitting this function, T and t are anomalies to their respective values in year 2015. Fitting is done using maximum a posteriori estimation.
