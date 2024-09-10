# PGW_ERA5_cmip6_WRF

procedure

1. remapping the cmip6 data to the grid of ERA5
2. adding the warming signal (cmip6_future-cmip6_present) into ERA5 data
3. applying function creat_FILT_for_WPS to write the warming signal into WRF simulation after the step ./ungrib.exe
