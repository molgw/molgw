!!****m* DoNOF/m_noft_driver
!! NAME
!!  m_noft_driver
!!
!! FUNCTION
!!  Module prepared to perform all procedures required for occ. and orbital optmization  
!!
!!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  Nfrozen |            Npairs_p_sing     |              Nvirt                               = NBF
!!  Nfrozen |         Npairs + Nsingleocc  |     Ncoupled*Npairs                   + Nempty   = NBF
!!                           | Nsingleocc  |   NBF_occ - Npairs_p_sing - Nfrozen   | Nempty   = NBF
!!                           Nbeta         Nalpha                                  NBF_occ
!!- - - - - - - - - - - - - - -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
!!
!! PARENTS
!!
!! CHILDREN
!!   m_optocc
!!   m_optorb
!!
!! SOURCE
module m_noft_driver

 use m_nofoutput
 use m_rdmd
 use m_integd
 use m_elag
 use m_optocc
 use m_optorb

 implicit none

 private :: read_restart,echo_input,occtogamma
!!***

 public :: run_noft
!!***

contains

!!***
!!****f* DoNOF/run_noft
!! NAME
!! run_noft
!!
!! FUNCTION
!!  Run NOFT procedures 
!!
!! INPUTS
!! INOF_in=PNOFi functional to use
!! Ista_in=Use PNOF7 (Ista_in=0) or PNOF7s (Ista_in=1)
!! NBF_occ_in=Number of orbitals that are occupied
!! Nfrozen_in=Number of frozen orbitals that remain with occ=2.0 
!! Npairs_in=Number of electron pairs
!! Ncoupled_in=Number of coupled orbitals per electron pair MINUS ONE
!! Nbeta_elect_in=Number of beta electrons (N/2 for spin compensated systems)
!! Nalpha_elect_in=Number of beta electrons (N/2 for spin compensated systems)
!! iERItyp_in=Index organization used for ERIs ({ij|lk}, <ij|kl>, and (ik|jl))
!! imethocc=Method used for OCC opt. CG (1) or L-BFGS (2)
!! imethorb=Method used to opt. orbs. currently only F_diag (1)
!! itermax=Max. number of global iters
!! iprintdmn=Print opt. 1,2-DMNs 
!! iprintints=Print ERIs in MO basis
!! itolLambda=Tol for Lambda_pq - Lambda_qp* is 10**-itolLambda
!! ndiis=Numb. of iter. used in orb. opt. to call DIIS
!! tolE_in=Tolerance on energy convergence
!! Vnn=Fixed nuclear-nuclear interaction energy
!! NO_COEF=Guess NO coefs (probably HF ones)
!! AOverlap_in= Overlap of atomic orbs. (matrix)
!! mo_ints=External subroutine that for given NO_COEF updates the hCORE and ERImol matrices
!! ofile_name=Name of the output file
!! lowmemERI=Logical parameter to decided whether to store only (NBF_tot,NBF_occ,NBF_occ,NBF_occ) part of the nat. orb. ERIs
!! restart=Logical parameter to decided whether we do restart
!! ireadGAMMAS, ireadOCC, ireadCOEF, ireadFdiag =Integer restart parameters that control the read of checkpoint files (true=1)
!! 
!! OUTPUT
!! Occ=Array containing the optimized occ numbers
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine run_noft(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,iERItyp_in,imethocc,imethorb,itermax,iprintdmn,iprintints,&
&  itolLambda,ndiis,Enof,tolE_in,Vnn,NO_COEF,AOverlap_in,Occ_inout,mo_ints,ofile_name,&
&  lowmemERI,restart,ireadGAMMAS,ireadOCC,ireadCOEF,ireadFdiag)   ! Optional
!Arguments ------------------------------------
!scalars
 integer,optional,intent(in)::ireadGAMMAS,ireadOCC,ireadCOEF,ireadFdiag
 logical,optional,intent(in)::restart,lowmemERI
 integer,intent(in)::INOF_in,Ista_in,imethocc,imethorb,itermax,iprintdmn,iprintints,itolLambda,ndiis
 integer,intent(in)::NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in
 integer,intent(in)::Nbeta_elect_in,Nalpha_elect_in,iERItyp_in
 real(dp),intent(in)::Vnn,tolE_in
 real(dp),intent(inout)::Enof
 external::mo_ints
!arrays
 character(len=100),intent(in)::ofile_name
 real(dp),dimension(NBF_tot_in),intent(inout)::Occ_inout
 real(dp),dimension(NBF_tot_in,NBF_tot_in),intent(in)::AOverlap_in
 real(dp),dimension(NBF_tot_in,NBF_tot_in),intent(inout)::NO_COEF
!Local variables ------------------------------
!scalars
 logical::ekt,diagLpL,restart_param
 integer::iorb,iter
 real(dp)::Energy,Energy_old,Vee,hONEbody
 type(rdm_t),target::RDMd
 type(integ_t),target::INTEGd
 type(elag_t),target::ELAGd
!arrays
 character(len=10)::coef_file
 character(len=200)::msg
 character(8)::date
 character(10)::time
 character(5)::zone
 integer,dimension(8)::tvalues
!************************************************************************

 diagLpL=.true.; restart_param=.false.;

 ! Initialize output
 call init_output(ofile_name)

 ! Write Header
 call write_header()

 ! Print user defined parameters used in this run
 if(present(restart)) then
  if(present(ireadGAMMAS).and.present(ireadCOEF).and.present(ireadFdiag).and.present(ireadOCC)) then
   restart_param=.true.
   call echo_input(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,imethorb,itermax,iprintdmn,&
&  iprintints,itolLambda,ndiis,tolE_in,restart=restart,ireadGAMMAS=ireadGAMMAS,&
&  ireadOCC=ireadOCC,ireadCOEF=ireadCOEF,ireadFdiag=ireadFdiag)
  else
   write(msg,'(a)') 'Warning! Asking for restart but the restart parameters are unspecified (not restarting).' 
   call write_output(msg)
   restart_param=.false.
   call echo_input(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,imethorb,itermax,iprintdmn,&
&  iprintints,itolLambda,ndiis,tolE_in)
  endif
 else 
  restart_param=.false.
  call echo_input(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,imethorb,itermax,iprintdmn,&
&  iprintints,itolLambda,ndiis,tolE_in)
 endif

 ! Initialize RDMd, INTEGd, and ELAGd objects.
 call rdm_init(RDMd,INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in,&
& Nbeta_elect_in,Nalpha_elect_in)
 if(present(lowmemERI)) then
  call integ_init(INTEGd,RDMd%NBF_tot,RDMd%NBF_occ,iERItyp_in,AOverlap_in,lowmemERI=lowmemERI)
 else
  call integ_init(INTEGd,RDMd%NBF_tot,RDMd%NBF_occ,iERItyp_in,AOverlap_in)
 endif
 call elag_init(ELAGd,RDMd%NBF_tot,diagLpL,itolLambda,ndiis,imethorb,tolE_in)

 ! Check for the presence of restart files. Then, if they are available read them (only if required)
 if(restart_param) then
  write(msg,'(a)') ' '
  call write_output(msg)
  call read_restart(RDMd,ELAGd,NO_COEF,ireadGAMMAS,ireadOCC,ireadCOEF,ireadFdiag)
  write(msg,'(a)') ' '
  call write_output(msg)
 endif

 ! Occ optimization using guess orbs. (HF, CORE, etc).
 write(msg,'(a)') ' '
 call write_output(msg)
 iter=-1;
 call mo_ints(RDMd%NBF_tot,RDMd%NBF_occ,INTEGd%NBF_jkl,NO_COEF,INTEGd%hCORE,INTEGd%ERImol)
 call INTEGd%eritoeriJK(RDMd%NBF_occ)
 call opt_occ(iter,imethocc,RDMd,Vnn,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K) ! Also iter=iter+1
 Energy_old=Energy

 ! Orb. and occ. optimization
 do
  ! Orb. optimization
  call ELAGd%clean_diis()
  call opt_orb(iter,imethorb,ELAGd,RDMd,INTEGd,Vnn,Energy,NO_COEF,mo_ints)
  if(imethorb==1) then ! For F diag method, print F_pp elements after each global iteration
   call ELAGd%print_Fdiag(RDMd%NBF_tot)
  endif
  call RDMd%print_orbs_bin(NO_COEF)

  ! Occ. optimization
  call opt_occ(iter,imethocc,RDMd,Vnn,Energy,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K) ! Also iter=iter+1
  call RDMd%print_gammas()

  ! Check convergence
  if(abs(Energy-Energy_old)<ELAGd%tolE) then
   Energy_old=Energy
   exit
  endif
  Energy_old=Energy
  
  ! Check maximum number of iterations
  if(iter>itermax) exit

 enddo

 ! Print optimized 1,2-RDMs
 if(iprintdmn==1) call RDMd%print_dmn(RDMd%DM2_J,RDMd%DM2_K) 

 ! Print hCORE and ERImol integrals in the last (opt) NO_COEF basis (if lowmemERI=.false. NBF_tot, otherwise only NBF_occ)
 if(iprintints==1) then
  if(present(lowmemERI)) then
   if(lowmemERI) then
    write(msg,'(a,f10.5,a)') 'Comment: Printing only the occ. ERI in the opt. nat. orb. basis'
    call write_output(msg)
   endif
  else
  endif
  call INTEGd%print_int()
 endif

 ! Print final diagonalized INTEGd%Lambdas values
 call ELAGd%diag_lag(RDMd,INTEGd,NO_COEF)

 ! Print final Extended Koopmans' Theorem (EKT) values
 if(RDMd%Nsingleocc==0) call ELAGd%diag_lag(RDMd,INTEGd,NO_COEF,ekt=ekt)

 ! Print optimal occ. numbers and save them in Occ_inout array
 write(msg,'(a)') ' '
 call write_output(msg)
 RDMd%occ(:)=2.0d0*RDMd%occ(:)
 write(msg,'(a,f10.5,a)') 'Total occ ',sum(RDMd%occ(:)),'. Optimized occ. numbers '
 call write_output(msg)
 do iorb=1,(RDMd%NBF_occ/10)*10,10
  write(msg,'(f12.6,9f11.6)') RDMd%occ(iorb:iorb+9)
  call write_output(msg)
 enddo
 iorb=(RDMd%NBF_occ/10)*10+1 
 write(msg,'(f12.6,*(f11.6))') RDMd%occ(iorb:) 
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)
 Occ_inout=0.0d0
 Occ_inout(1:RDMd%NBF_occ)=RDMd%occ(1:RDMd%NBF_occ)

 ! Print optimized nat. orb. coef.
 coef_file='NO_COEF'
 call RDMd%print_orbs(NO_COEF,coef_file)
 call RDMd%print_orbs_bin(NO_COEF)
 
 ! Print final Energy and its components (occs are already [0:2])
 hONEbody=0.0d0
 do iorb=1,RDMd%NBF_occ
  hONEbody=hONEbody+RDMd%occ(iorb)*INTEGd%hCORE(iorb,iorb)
 enddo
 Vee=Energy-hONEbody
 Enof=Energy+Vnn
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a,f15.6,a,i6,a)') 'Final NOF energy = ',Enof,' a.u. after ',iter,' global iter.'
 call write_output(msg)
 write(msg,'(a,f15.6,a)') 'hCORE            = ',hONEbody,' a.u.'
 call write_output(msg)
 write(msg,'(a,f15.6,a)') 'Vee              = ',Vee,' a.u.'
 call write_output(msg)
 write(msg,'(a,f15.6,a)') 'Vnn              = ',Vnn,' a.u.'
 call write_output(msg)
 write(msg,'(a)') ' '
 call write_output(msg)

 ! Free all allocated RDMd, INTEGd, and ELAGd arrays
 call ELAGd%free() 
 call INTEGd%free()
 call RDMd%free() 

 ! Write Footer
 call write_footer()
 
 ! Close unit 313 used for output file
 call close_output()

end subroutine run_noft
!!***

!!***
!!****f* DoNOF/echo_input
!! NAME
!! echo_input
!!
!! FUNCTION
!!  Echo all parameters employed in this run of this module
!!
!! INPUTS
!! INOF_in=PNOFi functional to use
!! Ista_in=Use PNOF7 (Ista_in=0) or PNOF7s (Ista_in=1)
!! NBF_occ_in=Number of orbitals that are occupied
!! Nfrozen_in=Number of frozen orbitals that remain with occ=2.0 
!! Npairs_in=Number of electron pairs
!! Ncoupled_in=Number of coupled orbitals per electron pair MINUS ONE
!! Nbeta_elect_in=Number of beta electrons (N/2 for spin compensated systems)
!! Nalpha_elect_in=Number of beta electrons (N/2 for spin compensated systems)
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine echo_input(INOF_in,Ista_in,NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,&
&  Ncoupled_in,Nbeta_elect_in,Nalpha_elect_in,imethocc,imethorb,itermax,iprintdmn,&
&  iprintints,itolLambda,ndiis,tolE_in,restart,ireadGAMMAS,ireadOCC,ireadCOEF,ireadFdiag)
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::restart
 integer,optional,intent(in)::ireadGAMMAS,ireadOCC,ireadCOEF,ireadFdiag
 integer,intent(in)::INOF_in,Ista_in,imethocc,imethorb,itermax,iprintdmn,iprintints,itolLambda,ndiis
 integer,intent(in)::NBF_tot_in,NBF_occ_in,Nfrozen_in,Npairs_in,Ncoupled_in
 integer,intent(in)::Nbeta_elect_in,Nalpha_elect_in
 real(dp),intent(in)::tolE_in
!arrays
!Local variables ------------------------------
!scalars
!arrays
 character(len=200)::msg
!************************************************************************
 write(msg,'(a)') ' '
 call write_output(msg)
 write(msg,'(a,i12)') ' NOF approximation in use          ',INOF_in
 call write_output(msg)
 if(INOF_in==7) then
  if(Ista_in==1) then
   write(msg,'(a,i12)') ' PNOF7s version selected Istat     ',Ista_in
  call write_output(msg)
  else
   write(msg,'(a,i12)') ' PNOF7  version selected Istat     ',Ista_in
  call write_output(msg)
  endif
 endif
 if(INOF_in==0) then
  write(msg,'(a)') ' Using Hartree-Fock approximation'
  call write_output(msg)
  write(msg,'(a)') ' L. Cohen and C. Frisberg, J. Chem. Phys, 65, 4234 (1976)'
  call write_output(msg)
 elseif(INOF_in==-1) then
  write(msg,'(a)') ' Using Muller-Baerends-Buijse approximation'
  call write_output(msg)
  write(msg,'(a)') ' A.M.K. Muller, Phys. Lett., 105A, 446 (1984)'
  call write_output(msg)
  write(msg,'(a)') ' M.A. Buijse and E.J. Baerends, Mol. Phys., 100, 401 (2002)'
  call write_output(msg)
 elseif(INOF_in==5) then
  write(msg,'(a)') ' Using PNOF5e approximation'
  call write_output(msg)
  write(msg,'(a)') ' M. Piris, X. Lopez, F. Ruiperez, J.M. Matxain, and J. M. Ugalde, J. Chem. Phys., 134, 164102 (2011)'
  call write_output(msg)
 elseif(INOF_in==7.and.Ista_in==0) then
  write(msg,'(a)') ' Using PNOF7(-) approximation'
  call write_output(msg)
  write(msg,'(a)') ' M. Piris, Phys. Rev. Lett., 119, 063002 (2017)'
  call write_output(msg)
  write(msg,'(a)') ' I. Mitxelena, M. Rodriguez-Mayorga, and M. Piris, Eur. Phys. J. B, 91, 109 (2018)'
  call write_output(msg)
 elseif(INOF_in==7.and.Ista_in==1) then
  write(msg,'(a)') ' Using PNOF7s approximation'
  call write_output(msg)
  write(msg,'(a)') ' M. Piris, Phys. Rev. A, 98, 022504 (2018)'
  call write_output(msg)
  write(msg,'(a)') ' M. Piris, Phys. Rev. A, 100, 032508 (2019)'
  call write_output(msg)
 else
  ! Nth
 endif
 write(msg,'(a,i12)') ' Numb. of basis functions          ',NBF_tot_in
 call write_output(msg)
 write(msg,'(a,i12)') ' Numb. of occ orbitals             ',NBF_occ_in
 call write_output(msg)
 write(msg,'(a,i12)') ' Numb. of frozen orbs (occ=2)      ',Nfrozen_in
 call write_output(msg)
 write(msg,'(a,i12)') ' Numb. of active e- pairs          ',Npairs_in
 call write_output(msg)
 write(msg,'(a,i12)') ' Numb. of "virtual" coupled orbs   ',Ncoupled_in
 call write_output(msg)
 write(msg,'(a,i12)') ' Numb. of singly occupied orbs     ',Nalpha_elect_in-Nbeta_elect_in
 call write_output(msg)
 if(imethocc==1) then
  write(msg,'(a,i12)') ' CG method used in occ opt.        ',imethocc
  call write_output(msg)
 else
  write(msg,'(a,i12)') ' L-BFGS method used in occ opt.    ',imethocc
  call write_output(msg)
 endif
 if(imethorb==1) then
  write(msg,'(a,i12)') ' F_diag method used in orb opt.    ',imethorb
  call write_output(msg)
  write(msg,'(a,e10.3)') ' Tolerance Lambda convergence        ',1.0d1**(-itolLambda)
  call write_output(msg)
  write(msg,'(a,i12)') ' Numb. of iter used in DIIS        ',ndiis
  call write_output(msg)
 else
  write(msg,'(a,i12)') ' Newton method used in orb opt.    ',imethorb
 call write_output(msg)
 endif
 write(msg,'(a,i11)') ' Max. number of global iterations   ',itermax
 call write_output(msg)
 write(msg,'(a,e10.3)') ' Tolerance Energy convergence        ',tolE_in
 call write_output(msg)
 write(msg,'(a,i12)') ' Print optimal 1,2-RDMs (true=1)   ',iprintdmn
 call write_output(msg)
 write(msg,'(a,i12)') ' Print last hCORE and ERImol ints  ',iprintints
 call write_output(msg)
 ! Check for the presence of restart files. If they are available, read them if required (default=not to read)
 if(present(restart)) then
  write(msg,'(a,i12)') ' Restart reading GAMMAs (true=1)   ',ireadGAMMAS
  call write_output(msg)
  write(msg,'(a,i12)') ' Restart reading DM1 file (true=1) ',ireadOCC
  call write_output(msg)
  write(msg,'(a,i12)') ' Restart reading COEFs  (true=1)   ',ireadCOEF
  call write_output(msg)
  if(imethorb==1) then
   write(msg,'(a,i12)') ' Restart reading F_pp   (true=1)   ',ireadFdiag
   call write_output(msg)
  endif
 endif
 write(msg,'(a)') ' '
 call write_output(msg)
 
end subroutine echo_input
!!***

!!***
!!****f* DoNOF/read_restart
!! NAME
!!  read_restart
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine read_restart(RDMd,ELAGd,NO_COEF,ireadGAMMAS,ireadOCC,ireadCOEF,ireadFdiag)
!Arguments ------------------------------------
!scalars
 integer,intent(in)::ireadGAMMAS,ireadOCC,ireadCOEF,ireadFdiag
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
!arrays
 real(dp),dimension(RDMd%NBF_tot,RDMd%NBF_tot),intent(inout)::NO_COEF
!Local variables ------------------------------
!scalars
 integer::iunit,istat,intvar,intvar1,icount
 real(dp)::doubvar
 real(dp),allocatable,dimension(:)::GAMMAS_in
 real(dp),allocatable,dimension(:,:)::NO_COEF_in
!arrays
 character(len=200)::msg
!************************************************************************

 ! Read NO_COEF for guess
 if(ireadCOEF==1) then
  allocate(NO_COEF_in(RDMd%NBF_tot,RDMd%NBF_tot))
  open(unit=iunit,form='unformatted',file='NO_COEF_BIN',iostat=istat,status='old')
  icount=0
  if(istat==0) then
   do
    read(iunit,iostat=istat) intvar,intvar1,doubvar
    if(istat/=0) then
     exit
    endif
    if(((intvar/=0).and.(intvar1/=0)).and.intvar*intvar1<=RDMd%NBF_tot*RDMd%NBF_tot) then
     NO_COEF_in(intvar,intvar1)=doubvar
     icount=icount+1
    else
     exit
    endif
   enddo
  endif
  if(icount==RDMd%NBF_tot*RDMd%NBF_tot) then
   NO_COEF(:,:)=NO_COEF_in(:,:)
   write(msg,'(a)') 'NO coefs. read from checkpoint file'
   call write_output(msg)
  endif
  close(iunit)
  deallocate(NO_COEF_in)
 endif

 ! Read GAMMAs (indep. parameters used to optimize occs.) from GAMMAS file
 if(ireadGAMMAS==1) then
  allocate(GAMMAS_in(RDMd%Ngammas))
  open(unit=iunit,form='unformatted',file='GAMMAS',iostat=istat,status='old')
  icount=0
  if(istat==0) then
   do 
    read(iunit,iostat=istat) intvar,doubvar
    if(istat/=0) then
     exit
    endif
    if((intvar/=0).and.intvar<=RDMd%Ngammas) then
     GAMMAs_in(intvar)=doubvar
     icount=icount+1
    else
     exit
    endif
   enddo
  endif
  if(icount==RDMd%Ngammas) then
   RDMd%GAMMAs_old(:)=GAMMAS_in(:)
   RDMd%GAMMAs_nread=.false.
   write(msg,'(a)') 'GAMMAs (indep. variables) read from checkpoint file'
   call write_output(msg)
  endif
  close(iunit)
  deallocate(GAMMAS_in)
 endif
 if(ireadOCC==1) then
  open(unit=iunit,form='unformatted',file='DM1',iostat=istat,status='old')
  icount=0
  if(istat==0) then
   do
    read(iunit,iostat=istat) intvar,intvar,doubvar
    if(istat/=0) then
     exit
    endif
    if((intvar/=0).and.intvar<=RDMd%NBF_occ) then
     RDMd%occ(intvar)=doubvar
     icount=icount+1
    else
     exit
    endif
   enddo
  endif
  if(icount==RDMd%NBF_occ) then
   if(ireadGAMMAS==1) then
    write(msg,'(a)') 'Comment: computing GAMMAs using occ. read from DM1 file'
    call write_output(msg)
   endif
   RDMd%occ(:)=0.5d0*RDMd%occ(:)
   call occtogamma(RDMd) 
   RDMd%occ=0.0d0   
   RDMd%GAMMAs_nread=.false.
   write(msg,'(a)') 'GAMMAs (indep. variables) calculated using DM1 file'
   call write_output(msg)
  endif
  close(iunit)
 endif

 ! Read diag. part of the F matrix for Lambda_pq - Lambda_qp* method
 if(ELAGd%imethod==1.and.ireadFdiag==1) then
  open(unit=iunit,form='unformatted',file='F_DIAG',iostat=istat,status='old')
  icount=0
  if(istat==0) then
   do 
    read(iunit,iostat=istat) intvar,doubvar
    if(istat/=0) then
     exit
    endif
    if((intvar/=0).and.intvar<=RDMd%NBF_tot) then
     ELAGd%F_diag(intvar)=doubvar
     icount=icount+1
    else
     exit
    endif
   enddo
  endif
  if(icount==RDMd%NBF_tot) then
   ELAGd%diagLpL=.false.
   write(msg,'(a)') 'F_pp elements read from checkpoint file'
   call write_output(msg)
  endif
  close(iunit)
 endif

end subroutine read_restart
!!***

!!***
!!****f* DoNOF/occtogamma
!! NAME
!!  occtogamma
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine occtogamma(RDMd)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!arrays
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5,iorb6,iorb7,iorb8
 real(dp)::argum,hole_iorb1 
!arrays
 real(dp),allocatable,dimension(:)::Holes
!************************************************************************

 allocate(Holes(RDMd%Npairs*(RDMd%Ncoupled-1)))
 do iorb=1,RDMd%Npairs
  iorb1 = RDMd%Nfrozen+iorb                                           ! iorb1=RDMd%Nfrozen+1,RDMd%Nbeta_elect
  RDMd%GAMMAs_old(iorb) = dacos(dsqrt(2.0d0*RDMd%occ(iorb1)-1.0d0))
  if(RDMd%Ncoupled/=1) then
   iorb2 = (RDMd%Ncoupled-1)*(iorb-1)+1
   iorb3 = (RDMd%Ncoupled-1)*iorb
   hole_iorb1 = 1.0d0 - RDMd%occ(iorb1)
   Holes(iorb2:iorb3) = hole_iorb1
   do iorb6=1,RDMd%Ncoupled-1
    iorb4 = (RDMd%Ncoupled-1)*(iorb-1)+iorb6                          ! iorb4=1,RDMd%pairs*(RDMd%Ncoupled-1)
    iorb5 = RDMd%Npairs+iorb4                                         ! iorb5=RDMd%pairs+1,RDMd%pairs*RDMd%Ncoupled
    iorb1 = RDMd%Nalpha_elect+RDMd%Ncoupled*(RDMd%Npairs-iorb)+iorb6  ! iorb1=RDMd%Nalpha_elect+1,RDMd%Nalpha_elect+RDMd%Ncoupled*RDMd%pairs-1         
    if(Holes(iorb4)>0.0d0) then
     argum=dsqrt(RDMd%occ(iorb1)/Holes(iorb4))
     if(argum>1.0d0)argum=1.0d0
     RDMd%GAMMAs_old(iorb5)=dasin(argum)
    else
     RDMd%GAMMAs_old(iorb5) = 0.0d0
    endif
    if(iorb6<RDMd%Ncoupled-1) then
     do iorb7=1,RDMd%Ncoupled-1-iorb6
      iorb8 = iorb4+iorb7                                             ! iorb4 < iorb8 < iorb*(RDMd%Ncoupled-1)
      Holes(iorb8) = Holes(iorb8) - RDMd%occ(iorb1)
     enddo
    endif
   enddo
  endif
 enddo
 deallocate(Holes)

end subroutine occtogamma
!!***

end module m_noft_driver
!!***
