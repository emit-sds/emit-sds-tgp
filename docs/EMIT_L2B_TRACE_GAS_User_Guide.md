# Earth Surface Mineral Dust Source Investigation (EMIT)

## Level 2B Greenhouse Gas Data Product User Guide

**Version:** 2.0 </br>
**Release Date:** TBD </br>
**JPL- D-107862** </br>

Jet Propulsion Laboratory 
California Institute of Technology 
Pasadena, California 91109-8099 


**Change Log**
| Version | Date       | Comments |
|---------|------------|----------|
| 0.0     | August 30, 2023 | Initial Draft |
| 0.1     | June 25, 2024 | CO2 Product Additions |
| 0.2     | March 10, 2025 | Per-pixel Uncertainty and Sensitivity Estimates (V002 Enhancements) |
| 2.0     | November 20, 2025 | V002 Plume Estimates |

<div style="page-break-after: always;"></div>

## 1	Introduction

### 1.1	Identification
This document describes information about the file structure and datasets provided in the EMIT Level 2B Greenhouse Gas data products. The algorithms and data content of the Level 2B CH4ENH, CH4PLM. CO2ENH, CO2PLMdata products are described briefly in this guide, with the purpose of providing the user with sufficient information about the content and structure of the data files to enable the user to access and use the data, in addition to understanding the uncertainties involved in the products.  More detail on the detection methods used are available in the [L2B GHG ATBD](EMIT_L2B_TRACE_GAS_ATBD.md).

### 1.2	Overview
Mineral dust aerosols originate as soil particles lifted into the atmosphere by wind erosion.  Mineral dust created by human activity makes a large contribution to the uncertainty of direct radiative forcing (RF) by anthropogenic aerosols (USGCRP and IPCC). Mineral dust is a prominent aerosol constituent around the globe. However, we have poor understanding of its direct radiative effect, partly due to uncertainties in the dust mineral composition. Dust radiative forcing is highly dependent on its mineral-specific absorption properties. The current range of iron oxide abundance in dust source models translates into a large range of values, even changing the sign of the forcing (-0.15 to 0.21 W/m2) predicted by Earth System Models (ESMs) (Li et al., 2020). The National Aeronautics and Space Administration (NASA) recently selected the Earth Mineral Dust Source Investigation (EMIT) to close this knowledge gap. EMIT was launched on July 14, 2022 to the International Space Station (ISS) to directly measure and map the soil mineral composition of critical dust-forming regions worldwide.

In addition to its primary objective described above, EMIT has demonstrated the capacity to characterize methane (CH4) and carbon dioxide (CO2) point source emissions by measuring gas absorption features in the shortwave infrared.  This document breaks from the other mission Algorithm Theoretical Basis Documents (ATBDs), as the CH4 and CO2 products are not part of the standard series of interconnected science products.  Readers should consult the L1B ATBD for the precursor products to what are used here.

The EMIT Project is part of the Earth Venture-Instrument (EV-I) Program directed by the Program Director of the NASA Earth Science Division (ESD). EMIT is comprised of a VSWIR Infrared Dyson imaging spectrometer adapted for installation on the International Space Station (ISS). EMIT measures radiance between 380 and 2500 nanometers, with an approximate 7 nm bandpass.  Data are collected in a swath that is approximately 75 km wide at the equator, with an approximate ground sampling distance of 60 m. 


### 1.3 Product Overview

Shortly after completing initial data validation, it became evident that EMIT was a particularly useful tool for mapping out greenhouse gases, including methane, carbon dioxide, and water vapor.  This is consistent with previous findings from airborne data, but global nature, revisit frequency and wide swath of EMIT provided an unprecedented opportunity to investigate greenhouse gas retrievals.  Early prototypes of this data appeared on the VSWIR Imaging Spectroscopy Interface for Open Science (VISIONS; https://earth.jpl.nasa.gov/emit/data/data-portal/Greenhouse-Gases/), and based on high demand the EMIT GHG product suite was developed.

There are currently four greenhouse gas data products being produced by EMIT, two for methane and two for carbon dioxide.  The methane products are the EMIT L2B Methane Enhancement Data product (EMIT L2B CH4 ENH) and the EMIT L2B Estimated Methane Plume Complexes (EMIT L2B CH4 PLM). The carbon dioxide products are the EMIT L2B Carbon Dioxide Enhancement Data product (EMIT L2B CO2 ENH) and the EMIT L2B Estimated Carbon Dioxide Plume Complexes (EMIT L2B CO2 PLM).  All of these products are provided as Cloud Optimized GeoTIFFs (COGs), with accompanying GeoJSON metadata.

More products, include per-pixel uncertainties and emission rate estimates are planned for the near future.


### 1.4 File Formats
#### 1.4.1 Metadata Structure

EMIT operates from the ISS, orbiting Earth approximately 16 times in a 24-hour period. EMIT starts and stops data recording based on a surface coverage acquisition mask. The top-level metadata identifier for EMIT data is an orbit, representing a single rotation of the ISS around Earth. A period of continuous data acquisition within an orbit is called an orbit segment, where each orbit segment can cover up to thousands of kilometers down-track depending on the acquisition mask map. Each orbit segment is broken into scenes of 1280 down-track lines for convenience, though scenes may be seamlessly reassembled into orbit segments. To prevent a very small number of lines in any scene, the last scene can extend up to 2559 lines.

Because plume complexes are not constrained to a single scene, the CH4PLM product does not utilize scene or orbit numbers – instead, a unique, global plume complex identifier is utilized.  The plume complex identifier does not monotonically increase nor should it be used to infer information about when a plume complex was observed or identified.  In the case of CH4PLM, the timestamp used is the timestamp associated with the earliest scene that a plume intersects.

The EMIT Greenhouse Gas collections (EMITL2BCH4ENH, EMITL2BCH4PLM, EMITL2BCO2ENH, and EMITL2BCO2PLM) contain a combination of COG, PNG, and GeoJSON files, as described in Table 1-1. 

Table 1-1: EMIT Greenhouse Gas collection file list and naming convention
Collection: EMIT L2B Methane Enhancement Data – EMITL2BCH4ENH
Methane Enhancements: 

<table>
<tr>
    <td> File </td>  <td>  Description  </td>
</tr>
<tr>
<tr><td colspan="2"><b>Collection: EMIT L2B Methane Enhancement Data – EMITL2BCH4ENH</b></td></tr>
<td>EMIT_L2B_CH4ENH_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif</td> <td>Methane Enhancement</td>
</tr>
<tr>
<td>EMIT_L2B_CH4UNCERT_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif</td> <td>Methane Uncertainty</td>
</tr>
<tr>
<td>EMIT_L2B_CH4SENS_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif</td> <td>Methane Sensitivity</td>
</tr>
<tr>
<td>EMIT_L2B_CH4ENH_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.png</td> <td>Browse</td>
</tr>
<tr><td colspan="2"><b>Collection: EMIT L2B Methane Plume Complexes – EMITL2BCH4PLM</b></td></tr>
<tr>
<td>EMIT_L2B_CH4PLM_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;IIIIII&gt;.tif</td> <td> Methane Plume Complex Data</td>
</tr>
<tr>
<td>EMIT_L2B_CH4PLM_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;IIIIII&gt;.json</td> <td> Methane Plume Complexe Metadata, Emissions, and Uncertainty</td>
</tr>
<td>EMIT_L2B_CH4PLM_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;IIIIII&gt;.png</td> <td> Methane Plume Complex Browse</td>



<tr><td colspan="2"><b>Collection: EMIT L2B Carbon Dioxide Enhancement Data – EMITL2BCO2ENH</b></td></tr>
<td>EMIT_L2B_CO2ENH_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif</td> <td>Carbon Dioxide Enhancement</td>
</tr>
<tr>
<td>EMIT_L2B_CO2UNCERT_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif</td> <td>Carbon Dioxide Uncertainty</td>
</tr>
<tr>
<td>EMIT_L2B_CO2SENS_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.tif</td> <td>Carbon Dioxide Sensitivity</td>
</tr>
<tr>
<td>EMIT_L2B_CO2ENH_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;OOOOO&gt;_&lt;SSS&gt;.png</td> <td>Browse</td>
</tr>
<tr><td colspan="2"><b>Collection: EMIT L2B Carbon Dioxide Plume Complexes – EMITL2BCO2PLM</b></td></tr>
<tr>
<td>EMIT_L2B_CO2PLM_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;IIIIII&gt;.tif</td> <td> Carbon Dioxide Plume Complex Data</td>
</tr>
<tr>
<td>EMIT_L2B_CO2PLM_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;IIIIII&gt;.json</td> <td> Carbon Dioxide Plume Complexe Metadata, Emissions, and Uncertainty</td>
</tr>
<tr>
<td>EMIT_L2B_CO2PLM_&lt;VVV&gt;_&lt;YYYYMMDDTHHMMSS&gt;_&lt;IIIIII&gt;.png</td> <td> Carbon Dioxide Plume Complex Browse</td>
</tr>
</table>

\<VVV> gives the software version number, e.g., 001
\<YYYYMMDDTHHMMSS> is a time stamp, e.g., 20220101T083015
\<OOOOO> is the orbit identification number, e.g., 12345
\<SSS> is the scene identification number, e.g., 007
\<IIIIII> gives a unique global plume identifier for the plume complex, e.g. 000120 


#### 1.4.2 L2B GHG Data Products

The EMIT L2B Greenhouse Gas collections contain estimated greenhouse gas enhancements, GHG sensitivities, and GHG pixel uncertainties.  Files are provided either as COGs, which are orthorectified (latitude/longitude, projected using WGS 84, EPSG:4326) using nearest neighbor resampling, consistent with geometric lookup tables provided in other EMIT products, or as GeoJSONs (for vector layers). 


Table 1-2: EMIT L2B Data Products Summary

<table>
<tr><td>Earth Science Data Type</td><td>Product Level</td><td>Data Dimension</td><td>Spatial Resolution</td><td>Swath Width</td><td>Map Projection</td></tr>

<tr><td colspan="6"><b>Collection: EMIT L2B Methane Enhancement Data – EMITL2BCH4ENH</b></td></tr>
<tr><td>Methane Enhancement Data</td><td>L2B</td><td>x, y, 1</td><td>60 m*</td><td>72 km*</td><td>WGS-84, EPSG:4326</td></tr>
<tr><td>Methane Pixel Uncertainty</td><td>L2B</td><td>x, y, 1</td><td>60 m*</td><td>72 km*</td><td>WGS-84, EPSG:4326<td></tr>
<tr><td>Methane Pixel Sensitivity </td ><td> L2B </td><td> x, y, 1 </td><td> 60 m* </td ><td> 72 km* </td><td> WGS-84, EPSG:4326 </td></tr>

<tr><td colspan="6"><b>>Collection: EMIT L2B Methane Plume Complexes – EMITL2BCH4PLM</b></td></tr>
<tr><td>Methane Plume Complexes</td><td>L2B</td><td>x, y, 1</td><td>60 m*</td><td>72 km*</td><td>WGS-84, EPSG:4326</td></tr>

<tr><td>Earth Science Data Type</td><td>Product Level</td><td>Data Dimension</td><td>Spatial Resolution</td><td>Swath Width</td><td>Map Projection</td></tr>

<tr><td colspan="6"><b>Collection: EMIT L2B Carbon Dioxide Enhancement Data – EMITL2BCO2ENH</b></td></tr>
<tr><td>Carbon Dioxide Enhancement Data</td><td>L2B</td><td>x, y, 1</td><td>60 m*</td><td>72 km*</td><td>WGS-84, EPSG:4326</td></tr>
<tr><td>Carbon Dioxide Pixel Uncertainty</td><td>L2B</td><td>x, y, 1</td><td>60 m*</td><td>72 km*</td><td>WGS-84, EPSG:4326<td></tr>
<tr><td>Carbon Dioxide Pixel Sensitivity </td ><td> L2B </td><td> x, y, 1 </td><td> 60 m* </td ><td> 72 km* </td><td> WGS-84, EPSG:4326 </td></tr>

<tr><td colspan="6"><b>>Collection: EMIT L2B Carbon Dioxide Plume Complexes – EMITL2BCO2PLM</b></td></tr>
<tr><td>Carbon Dioxide Plume Complexes</td><td>L2B</td><td>x, y, 1</td><td>60 m*</td><td>72 km*</td><td>WGS-84, EPSG:4326</td></tr>



<tr><tr>
</table>
* Nominal at equator

### 1.5 Product Availability

The EMIT L2B Greenhouse Gas products will be available at the NASA Land Processes Distributed Active Archive Center (LP DAAC, https://lpdaac.usgs.gov/) and through NASA Earthdata (https://earthdata.nasa.gov/).

## 2. Greenhouse Gas Products

The EMIT Level 2B Greenhouse Gas products are a series of products that can be used to identify and quantify point source emissions. The first step is an enhancement estimate, which is fundamentally based on an adaptive matched filter approach. This yields an estimate of GHG enhancement in parts per million meter (ppm m), a total column enhancement estimate. This is provided at the EMITL2BCH4ENH product for methane and EMITL2BCO2ENH for carbon dioxide. Next, individual plumes are identified and vetted by multiple scientists, and high confidence plume complexes are provided for methane in the EMITL2BCH4PLM product and for carbon dioxide in the EMITL2BCO2PLM product. All datasets are provided as cloud optimized GeoTIFFs (COGs), with some supporting metadata provided in the EMITL2BCH4PLM and EMITL2BCO2PLM GeoJSON files.  

### 2.1 Delivery Frequency

Recognizing the value of both high confidence and low latency, EMIT products look to strike a balance on delivery.  Matched filter results are computed on the EMIT SDS on a daily basis, following the generation of L1B and L2A products. The ENH products for all scenes are then delivered to the LP DAAC and ingested on a continuous basis, similar to the L1B and L2A products.  However, a manual review of each scene is required before plumes are identified.  This can introduce additional latency from the time the time the initial processing is complete.  In general, most scenes are available to the EMIT SDS within a week of observation.  Manual plume complex identification can introduce another week or so of latency, after which plumes are delivered.  Once a plume complex is identified and confirmed by the requisite three scientists (see ATBD for details), the plume complex is sent to the LP DAAC as soon as operationally viable. 

### 2.2 File Structure

The L2BCH4ENH, L2BCH4PLM, L2BCO2ENH, and L2BCO2PLM are provided as COGs, and uncertainty estimates are provided along with the CH4PLM and CO2PLM products in the GeoJSON metadata. Products are single band and vary in size; L2BCH4ENH and L2BCO2ENH products follow EMIT scene sizes but are projected and so vary in size.  L2BCH4PLM and L2BCO2PLM data depend on the size of an individual plume complex and vary substantially more.  L2BCH4PLM and L2BCO2PLM data also include GeoJSON files, which include an outline of each plume complex along with other metadata.  Each GeoTIFF includes a variety of metadata fields to specify additional parameters.



### 3 References

* L. Li, N.M. Mahowald, R.L. Miller, C. Pérez García-Pando, M. Klose, D.S. Hamilton, M. Gonçalves Ageitos, P. Ginoux, Y. Balkanski, R.O. Green, O. Kalashnikova, J.F. Kok, V. Obiso, D. Paynter, and D.R. Thompson: Quantifying the range of the dust direct radiative effect due to source mineralogy uncertainty, Atmos. Chem. Phys., 21, 3973–4005, https://doi.org/10.5194/acp-21-3973-2021 (2021).

### 4	Acronyms
| Acronym | Definition |
|---------|------------|
| ADC | Analog to Digital Converter |
| APID | Application Identifier |
| ASCII | American Standard Code for Information Interchange |
| BIL | Band Interleaved by Line |
| CCSDS | Consultative Committee for Space Data Systems |
| DAAC | Distributed Active Archive Center |
| DCID | Data Collection Identifier |
| DN | Digital Number |
| EMIT | Earth Mineral dust source InvesTigation |
| ENVI | Environment for Visualizing Images |
| ESDIS | Earth Science Data and Information System |
| ESM | Earth System Model |
| FPA | Focal Plane Array |
| FPGA | Field Programmable Gate Array |
| FPIE | Focal Plane Interface Electronics |
| FPIE-A | Focal Plane Interface Electronics - Analog |
| FSW | Flight Software |
| Gbps | Gigabits per second |
| GLT | Geometry Lookup Table |
| HOSC | Huntsville Operations and Support Center |
| ICD | Interface Control Document |
| IOS | Instrument Operations System |
| ISS | International Space Station |
| JPL | Jet Propulsion Laboratory |
| kHz | Kilohertz |
| L0 | Level 0 (compressed, raw packets) |
| L1A | Level 1A (reconstructed, uncompressed data reassembled into scenes) |
| L1B | Level 1B (calibrated radiances with geolocation parameters) |
| L2A | Level 2A (atmospherically-corrected surface reflectance) |
| L2B | Level 2B (mineral feature depth maps) |
| L3 | Level 3 (gridded global map of mineral composition and abundances) |
| L4 | Level 4 (model runs of GISS ModelE2 and NCAR CESM) |
| LP DAAC | Land Processes Distributed Active Archive Center |
| LSB | Least Significant Bit |
| MSB | Most Significant Bit |
| NASA | National Aeronautics and Space Administration |
| NetCDF | Network Common Data Format |
| PGE | Product Generation Executable |
| PLRA | Program Level Requirements Appendix |
| ROIC | Readout Integrated Circuit |
| SDS | Science Data System |
| SIS | Software Interface Specification |
| SSDR | Solid State Data Recorder |
| UTC | Universal Time Coordinated |
