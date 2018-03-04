#!/bin/bash

run=4
angles=(  4 7 8 9  10 11 12  13 14 15 16 )
#angles=( 7 9  12 13 14)
system=S2
amp1=
layers=10
wavelength=1000000
output_steps=1001
output_id=0
output_intervals=50000
output_id_step=${#angles[@]}
type=rd_sin_top_decay
#type=random
csvs_file=${type}_${system}_${amp}g_${wavelength}WL_${run}_input.csv

for angle in ${angles[*]}
do
	out_file=${type}_${system}_${angle}_${amp}g_${wavelength}WL_${run}
	act_dir=${type}_${system}_${angle}_${amp}gy_${amp}gx_${wavelength}WLy_${wavelength}WLx_${run}
	case "$1" in
		run)
			if pgrep -f $out_file > /dev/null 
			then
				echo -e $out_file' is already running'
			else
				echo -e 'starting '$out_file
				sed -e 's/%run%/'${run}'/g' -e "s/%layers%/$layers/" -e "s/%angle%/$angle/" -e "s/%system%/$system/" -e "s/%my%/$amp/" -e "s/%mx%/$amp/" -e "s/%wavelength_y%/$wavelength/g" -e "s/%wavelength_x%/$wavelength/g" -e "s/%wavelength%/$wavelength/g" configs/${type}_template.xml > configs/${out_file}.xml
				nohup ./GranFrixrm configs/$out_file.xml > consoles/$out_file.out &
				
			fi
			;;
		 res_run)
			
                        if pgrep -f $out_file_r > /dev/null
                        then
                                echo -e $out_file_r' is already running'
                        else
                                echo -e 'starting '$out_file_r
                                sed -e 's/%run%/'${run}'/g' -e "s/%layers%/$layers/" -e "s/%angle%/$angle/" -e "s/%system%/$system/" -e "s/%my%/$my/" -e "s/%mx%/$mx/" -e "s/%wavelength_y%/$wavelength_y/" -e "s/%wavelength_x%/$wavelength_x/" configs/${type}_restart_template.xml > configs/${out_file_r}.xml
                                nohup ./GranFrixrm configs/$out_file_r.xml > consoles/$out_file_r.out &

                        fi
                        ;;

		kill)
			echo -e 'killing '$out_file
			pkill -9 -f $out_file
			;;
		skip)
			echo -e 'skipping '$act_file
			echo ' ' > $act_dir/skip
			;;
		output)
			echo -e 'outputting '$act_file
			echo -e $output_steps'\n'$output_intervals'\n'$output_id'\n'$output_id_step > $act_dir/output
			((output_id++))
			;;
		csvs)
			echo -e 'copying csvs '$act_dir
			cp -f $out_file/vol_flux*.csv csvs
			cp -f $out_file/veloci*.csv csvs
			echo ${angle}','${out_file}','500000000 ',' 1300000000 >> csvs/$csvs_file
			
			;;
		clear)
			echo -e 'clearing '$act_dir
			rm -rf $act_dir
			;;
	esac
done

