set -x


folder_name='COMAP_to_healpix'
mkdir "$folder_name"

mkdir "./COMAP_to_healpix/Fits_files"
mkdir "./COMAP_to_healpix/Plots"
mkdir "./COMAP_to_healpix/TXT_files"

for file in ./Galactic_2021-10-19/*; do
    python COMAP_to_healpix.py "$file" "$(basename $file)"
done
