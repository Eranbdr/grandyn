#!/bin/bash

run=1
angles=(13)
system=S2
my=11
mx=11
layers=10
wavelength=8000
output_steps=100
output_id=0
output_intervals=10000
output_id_step=${#angles[@]}
type=random
csvs_file=${type}_${system}_${my}gy_${mx}gx_${wavelength}WL_${run}_input.csv

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
				sed -e's/%layers%/'${layers}'/g' -e 's/%run%/'${run}'/g' -e "s/%angle%/$angle/" -e "s/%system%/$system/" -e "s/%my%/$my/" -e "s/%mx%/$mx/" -e "s/%wavelength%/$wavelength/" configs/${type}_template.xml > configs/${out_file}.xml
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
			cp -f $out_file/*.csv csvs
			echo ${angle}','${out_file}','500100000 ',' 2000100000 >> csvs/$csvs_file	
			;;
	esac
done

