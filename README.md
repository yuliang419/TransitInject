# TransitInject

This is a set of tools that allow users to download and select M dwarf light curves from K2, inject transits and probe their recovery rate with BLS or other methods.

## Pre-Requisites
    numpy 1.11.1 
    scipy 0.18.1
    joblib 0.9.4
    matplotlib 1.5.3
    pandas 0.19.1
    ldtk
    batman 0.9.0
    random
    urllib 1.17
    lxml
    FFITools
    vartools 1.33


## Usage

These are examples of how to use TransitInject functions.

### Downloading light curves from K2
Download csv file of campaign targets from https://keplerscience.arc.nasa.gov/k2-approved-programs.html

Run getfiles.py after editing the script to contain the desired proposal numbers. Currently the file is configured to download all M dwarfs from campaigns 1-3.

### Injecting transits into downloaded K2 light curves
Create list of targets to be injected with transits.
    
    cd k2/k2mdwarfs
    ls [desired files] >mdwarfs.ls

Then run InjectionTools.py

### Transit recovery with BLS
Not yet implemented.

## Copyright and License

Copyright(c) Liang Yu, all rights reserved

Distributed under the MIT License.