@echo off

set "folder_name=COMAP_to_healpix"
mkdir "%folder_name%"

for %%f in (./Galactic_2021-10-19\*) do (
    python COMAP_to_healpix.py "%%f" "%%~nxf"

)
