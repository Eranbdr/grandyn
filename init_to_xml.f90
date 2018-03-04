program init_to_xml
   use xml_data_config_template
   integer :: f
   
   
   call write_xml_file_config_template( 'config/out_config.xml', 2000 )
end program