
# Solar Flares Display

Easy access to relevant images and information of any solar flare the past 12 years using python.

## Purpose

A solar flare display that easily accesses data and images using a common programming language. Most scientists use IDL— a more expensive method— to get solar data.

*PLEASE WAIT 15-20 SECONDS. THERE SHOULD BE 3 VIDEOS APPEARING BELOW*
![page1](https://user-images.githubusercontent.com/100271213/183282760-9c4211a5-a795-45d3-9bab-3f54a295ea9a.gif)
![page2](https://user-images.githubusercontent.com/100271213/183282807-a2605831-83a0-426b-be25-714270332dbe.gif)
![page3](https://user-images.githubusercontent.com/100271213/183282840-e1ad0e9d-ca34-4c5d-a4fc-00198797982c.gif)

## Control Panel
<img width="1398" alt="Screen Shot 2022-08-07 at 1 46 24 AM" src="https://user-images.githubusercontent.com/100271213/183282874-e1bcabf3-d62f-4333-a5bc-43ac83e1f39a.png">

## Features
- High speed solar movies (60 FPS)
- Priortizes automatic exposure control (AEC) images and unsaturated images for movie
- Displays flare contours of 1600 angstrom images
- Calculates area and location of HMI magnetic regions
- Running difference, base difference, running ratio and base difference movies are availiable for each wavelength
- Complex movie and graph controls
- Interchangable widgets with smooth zoom in/out and pan

## Data Types
- 5 AIA Wavelengths (131, 171, 304, 335, 1600)
- HMI Magnetogram
- Reconstructed RHESSI Images
- GOES X-ray Flux (XRSA, XRSB)
- RHESSI TimeSeries (9 wavelengths)
- EUV Flux for each AIA wavelength
- Night/SAA Time Range Data
- CME Info
- Flare Class Info
- Hale Class Info

## Databases Used
- NASA GSFC/RHESSI
- Stanford JSOC
- LMSAL SolarSoft
- NOAA SWPC

## New Flare Detection System
- Developed a numerical differentiation algorithm to more accurately identify each flare's start, peak and end time
- Recalibrated GOES flux curves to better estimate each flare's class and magnitude
<img width="827" alt="Screen Shot 2022-02-24 at 2 01 43 AM" src="https://user-images.githubusercontent.com/100271213/155502444-af072834-c333-4b64-b8c2-0d7ef491dda0.png">
 
## What's inside the Zip File?
- Cleaned GOES TimeSeries data per year from 2010-2018 that eliminates all gaps by combining GOES 13, 14, 15 data
- RHESSI TimeSeries from 2010-2018
- CME Info

## How do I run this program and access a specific flare?
- Download entire repository
- Click on the google drive link and unzip the file(it contains 8 years worth of data, so it's a big file)
- Run create_flare_master_list.py
- Input starting date, end date and click "Confirm Selections"
- Run widget_controller.py
- Click on any flare, then click "Show Flare"

## Documentation - In Progress

## License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
