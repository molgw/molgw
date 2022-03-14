#!/bin/bash

name=`git rev-parse HEAD`

echo "!!****m* DoNOF/m_gitver" > m_gitver.F90
echo "!! NAME" >> m_gitver.F90
echo "!!  m_gitver" >> m_gitver.F90
echo "!!" >> m_gitver.F90
echo "!! FUNCTION" >> m_gitver.F90
echo "!!  Module to print git version information ">> m_gitver.F90
echo "!! ">> m_gitver.F90
echo "!! COPYRIGHT ">> m_gitver.F90
echo "!! This file is distributed under the terms of the ">> m_gitver.F90
echo "!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .">> m_gitver.F90
echo "!! ">> m_gitver.F90
echo "!! ">> m_gitver.F90
echo "!! PARENTS ">> m_gitver.F90
echo "!! ">> m_gitver.F90
echo "!! CHILDREN ">> m_gitver.F90
echo "!! " >> m_gitver.F90
echo "!! " >> m_gitver.F90
echo "!! SOURCE " >> m_gitver.F90
echo "  " >> m_gitver.F90
echo "module m_gitver " >> m_gitver.F90 
echo " " >> m_gitver.F90
echo " implicit none" >> m_gitver.F90  
echo " " >> m_gitver.F90  
echo "! private ::" >> m_gitver.F90  
echo "!!***" >> m_gitver.F90  
echo " " >> m_gitver.F90  
echo " public :: gitversion" >> m_gitver.F90  
echo "!!***" >> m_gitver.F90  
echo " " >> m_gitver.F90  
echo "contains " >> m_gitver.F90  
echo "!!*** " >> m_gitver.F90  
echo " " >> m_gitver.F90  
echo "!!****f* DoNOF/gitversion" >> m_gitver.F90  
echo "!! NAME " >> m_gitver.F90  
echo "!! gitversion" >> m_gitver.F90  
echo "!!" >> m_gitver.F90  
echo "!! FUNCTION" >> m_gitver.F90  
echo "!!  Get the SHA string from Github version " >> m_gitver.F90 
echo "!! " >> m_gitver.F90
echo "!! INPUTS " >> m_gitver.F90
echo "!!   sha = string that on exit contains the SHA Github Key" >> m_gitver.F90
echo "!! " >> m_gitver.F90
echo "!! OUTPUT " >> m_gitver.F90
echo "!! " >> m_gitver.F90
echo "!! PARENTS " >> m_gitver.F90
echo "!!   " >> m_gitver.F90
echo "!! CHILDREN " >> m_gitver.F90
echo "!! " >> m_gitver.F90
echo "!! SOURCE " >> m_gitver.F90
echo " " >> m_gitver.F90
echo "subroutine gitversion(sha)" >> m_gitver.F90
echo "!Arguments ------------------------------------ " >> m_gitver.F90
echo "!scalars " >> m_gitver.F90
echo "!arrays " >> m_gitver.F90
echo "character(100),intent(inout)::sha" >> m_gitver.F90  >> m_gitver.F90
echo "!Local variables ------------------------------ " >> m_gitver.F90
echo "!scalars " >> m_gitver.F90
echo "!arrays " >> m_gitver.F90
echo "!************************************************************************ " >> m_gitver.F90
echo " " >> m_gitver.F90
echo "  write(sha,'(a)')'$name'" >> m_gitver.F90
echo " " >> m_gitver.F90
echo "end subroutine gitversion" >> m_gitver.F90
echo "!!***" >> m_gitver.F90
echo " " >> m_gitver.F90
echo "end module m_gitver"  >> m_gitver.F90
echo "!!***"  >> m_gitver.F90

