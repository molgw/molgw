#!/bin/bash

name=`git rev-parse HEAD`

echo "!!****m* DoNOF/gitver" > gitver.F90
echo "!! NAME" >> gitver.F90
echo "!!  gitver" >> gitver.F90
echo "!!" >> gitver.F90
echo "!! FUNCTION" >> gitver.F90
echo "!!  Subroutine to print git version information ">> gitver.F90
echo "!! ">> gitver.F90
echo "!! COPYRIGHT ">> gitver.F90
echo "!! This file is distributed under the terms of the ">> gitver.F90
echo "!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .">> gitver.F90
echo "!! ">> gitver.F90
echo "!! ">> gitver.F90
echo "!! PARENTS ">> gitver.F90
echo "!! ">> gitver.F90
echo "!! CHILDREN ">> gitver.F90
echo "!! " >> gitver.F90
echo "!! " >> gitver.F90
echo "!! SOURCE " >> gitver.F90
echo "  " >> gitver.F90
echo "!!****f* DoNOF/gitversion" >> gitver.F90  
echo "!! NAME " >> gitver.F90  
echo "!! gitversion" >> gitver.F90  
echo "!!" >> gitver.F90  
echo "!! FUNCTION" >> gitver.F90  
echo "!!  Get the SHA string from Github version " >> gitver.F90 
echo "!! " >> gitver.F90
echo "!! INPUTS " >> gitver.F90
echo "!!   sha = string that on exit contains the SHA Github Key" >> gitver.F90
echo "!! " >> gitver.F90
echo "!! OUTPUT " >> gitver.F90
echo "!! " >> gitver.F90
echo "!! PARENTS " >> gitver.F90
echo "!!   " >> gitver.F90
echo "!! CHILDREN " >> gitver.F90
echo "!! " >> gitver.F90
echo "!! SOURCE " >> gitver.F90
echo " " >> gitver.F90
echo "subroutine gitversion(sha)" >> gitver.F90
echo "!Arguments ------------------------------------ " >> gitver.F90
echo "!scalars " >> gitver.F90
echo "!arrays " >> gitver.F90
echo "character(100),intent(inout)::sha" >> gitver.F90  >> gitver.F90
echo "!Local variables ------------------------------ " >> gitver.F90
echo "!scalars " >> gitver.F90
echo "!arrays " >> gitver.F90
echo "!************************************************************************ " >> gitver.F90
echo " " >> gitver.F90
echo "  write(sha,'(a)')'$name'" >> gitver.F90
echo " " >> gitver.F90
echo "end subroutine gitversion" >> gitver.F90
echo "!!***" >> gitver.F90
echo " " >> gitver.F90
echo "!!***"  >> gitver.F90

