#!/bin/bash

run=4
angles=(  4 7  10  12  13 13.5 14 14.5 15 15.5 16 )
#angles=(7) 

if [ "${4}" = 'first' ] ; then
	angles=( 7 12 13.5 14.5 15.5  )
elif [ "${4}" = 'second' ] ; then 
	angles=( 4 10 13 14 15 16)
elif [ -n ${4} ] ; then
        angles=(${4})
fi

system=${5}
amp=${1}
layers=10
wavelength=${2}
output_steps=1001
output_id=0
output_intervals=50000
output_id_step=${#angles[@]}
#type=rd_xy_sin_top_decay
#type=rd_xy_sin_full_decay
type=random
csvs_file=${type}_${system}_${amp}g_${wavelength}WL_${run}_input.csv

for angle in ${angles[*]}
do
	out_file=${type}_${system}_${angle}_${amp}g_${wavelength}WL_${run}
	act_dir=${type}_${system}_${angle}_${amp}gy_${amp}gx_${wavelength}WLy_${wavelength}WLx_${run}
	case "$3" in
		run)
			if pgrep -f $out_file > /dev/null 
			then
				echo -e $out_file' is already running'
			else
				# replace all the place holdersin confg template
				echo -e 'starting '$out_file
				sed  -e 's/%run%/'${run}'/g' -e "s/%layers%/$layers/" -e "s/%angle%/$angle/" -e "s/%system%/$system/" -e "s/%my%/$amp/" -e "s/%mx%/$amp/" -e "s/%wavelength_y%/$wavelength/g" -e "s/%wavelength_x%/$wavelength/g" -e "s/%wavelength%/$wavelength/g" configs/${type}_template.xml > configs/${out_file}.xml
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
			echo ' ' > $out_file/skip
			;;
		output)
			echo -e 'outputting '$act_file
			echo -e $output_steps'\n'$output_intervals'\n'$output_id'\n'$output_id_step > $act_dir/output
			((output_id++))
			;;
		csvs)
			echo -e 'copying csvs '$act_dir
			cp -f $out_file/vol_flux*.csv csvs
			cp -f $out_file/flux*.csv csvs
			cp -f $out_file/veloci*.csv csvs
			echo ${angle}','${out_file}','500000000 ',' 2000000000 >> csvs/$csvs_file
			
			;;
		clear)
			echo -e 'clearing '$act_dir
			rm -rf $act_dir
			;;
	esac
done

