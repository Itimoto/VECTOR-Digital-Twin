# VECTOR Digital Twin
 Digital Twin for Vector Estimation for Continuous Tracking of Radio Signals

## Scripts of Interest
### `base_doa_pos.m`
The first iteration of the Digital Twin. Generates an 802.11az Null Data Packet, transmits it, estimates the Channel State Information, then determines the linear distance to a station.
Does not account for multipath.

## Master Documentation
This is the primary source of documentation, and is found under `/docs/README.md`. `/docs` also contains some useful official MATLAB examples that were used, then modified to fit the Digital Twin.

## System Requirements
- Developed on MATLAB R2023a
- Toolboxes:
    - Signal Processing Toolbox
    - Communications Toolbox
    - WLAN Toolbox