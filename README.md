# pitd_insert

insert PIOMAS ITD into GEOS5 sea ice restart
only works with GEOS-5 style sea ice (thermo) restart files 

## there are 2 versions:

* piomas_itd_insert.py: unoptimized, runs very slow

* piomas_itd_insert_fast.py: vectorized most loops, runs fast
