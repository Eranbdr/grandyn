#!/bin/bash

run=1
angles=( 4 7  10 12 13 13.5 14 14.5 15 15.5 16 )
system=S2
my=31
mx=31
layers=10
wavelength=12000
output_steps=1000
output_id=1
output_intervals=50000
output_id_step=${#angles[@]}
type=random
csvs_file=${type}_${system}_${my}gy_${mx}gx_${wavelength}WL_${run}_input.csv
csvs2_file=${type}_S12_${my}gy_${mx}gx_${wavelength}WL_${run}_input.csv

for angle in ${angles[*]}
do
	out_file=${type}_${system}_${angle}_${my}gy_${mx}gx_${wavelength}WL_${run}
	
	case "$1" in
		run)
			if pgrep -f $out_file > /dev/null 
			then
				echo -e $out_file' is already running'
			else
				echo -e 'starting '$out_file
				sed -e 's/%run%/'${run}'/g' -e "s/%layers%/$layers/" -e "s/%angle%/$angle/" -e "s/%system%/$system/" -e "s/%my%/$my/" -e "s/%mx%/$mx/" -e "s/%wavelength%/$wavelength/" configs/${type}_template.xml > configs/${out_file}.xml
				nohup ./GranFrixrm configs/$out_file.xml > consoles/$out_file.out &
				
			fi
			;;
		kill)
			echo -e 'killing '$out_file
			pkill -9 -f $out_file
			;;
		skip)
			echo -e 'skipping '$out_file
			echo ' ' > $out_file/skip
			;;
		output)
			((output_id++))
			echo -e 'outputting '$out_file
			echo -e $output_steps'\n'$output_intervals'\n'$output_id'\n'$output_id_step > $out_file/output
			;;
		csvs)
			echo -e 'copying csvs '$out_file
			cp -f $out_file/vol_flux*.csv csvs
			cp -f $out_file/veloci*.csv csvs
			echo ${angle}','${out_file}','500100000 ',' 2000100000 >> csvs/$csvs_file
			
			;;
		 csvs2)
                        echo -e 'setting double '$out_file
                        echo ${angle}','${type}_S1_${angle}_${my}gy_${mx}gx_${wavelength}WL_${run}','500100000 ',' 2000100000 >> csvs/$csvs2_file
			echo ${angle}','${type}_S2_${angle}_${my}gy_${mx}gx_${wavelength}WL_${run}','500100000 ',' 2000100000 >> csvs/$csvs2_file

                        ;;

		clear)
			echo -e 'clearing '$out_file
			rm -rf $out_file
			;;
	esac
done
rm -f csvs/*${type}_${system}_13_${my}gy_${mx}gx_${wavelength}WL_${run}.csv
