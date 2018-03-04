   if ( present(errout) ) errout = error
end subroutine

   write(info%lun,'(a)') '</config_template>'
   call xml_close(info)
end subroutine


end subroutine

end module
