use diver8

implicit none

   !**********  VARIABLES **********************
   !lectura ADN
   character*200   file_adn
   integer*8       file_size
   integer*8       nbases,enes
   character*2     tipo
   integer*8      ,allocatable::f(:)
   !Lista enlazada
   type corte
     integer*8    posic
     integer*8   ,allocatable::frec(:)
     integer*1    intact
     type(corte) ,pointer :: next,prev
   end type 
   type(corte)   ,pointer:: first,last
   !Lista segmentos
   type tseg
     integer*8    ini,fin
     type(corte) ,pointer::p
   end type
   type(tseg)    ,allocatable::S(:)
   !Enes
   integer*8      haygaps,ncontigs
   !segmentacion
   integer*8      nexper
   real*8         slini,slfin,sl     
   integer*8      new_segs,eliminados
   !salida resultados
   character*200  file_pro,file_segment_sizes   
   !varias
   integer*8       i,j  
   !test
   integer*8      nsegs,nsegs0,nsegs_good
   integer*8     ,allocatable::fcorte(:)
   
   
   integer*8      pcorte
   real*8         dmax,totdiv
   
   !***** GLOBALES **************   
   integer*8       ALF,m,ALF1   !tamano alfabeto y word size (m=0 para binario)
   integer*1      ,allocatable::ADN(:)
   logical         VERBOSE,WRITE_FRECS
   integer*8       INIADN,FINADN
   logical         HAY_RANDOM !Esta variable nos dice si se ha sustituido algun
                              !segmento de enes por nucleótidos aleatorios 
   integer*8      ,parameter::MIN_N_SIZE=0
   integer*8      NCORTES
   integer*4      UN          !unit for output
   !********************************************
   
   !*************** PARAMETROS *****************
   
   call get_parameters(file_adn,slini,slfin,nexper,ALF,tipo,file_pro,&
                       file_segment_sizes,INIADN,FINADN)
   ALF1=ALF-1
   if (ALF==4) then
     m=1
   else if (ALF==2) then
     m=0
   endif    
      
   !*** Opens file for the SCC profile ***************   
   if (file_pro/='') then
     UN = 7
     open(UN,file=trim(file_pro))
   else
     UN = 0
   endif  
   !***************************************************
   
   !**** Opens file for the list of segment sizes *****
   open(18,file=file_segment_sizes)
   !***************************************************
   
   
   !write(*,*) trim(file_adn)
   !write(*,*) slini,slfin
   !write(*,*) nexper
   !write(*,*) ALF,tipo
  
   !goto 999

   !************* LECTURA ADN ******************
   INQUIRE(FILE=file_adn, SIZE=file_size)
   if (FINADN==0) then
     allocate(ADN(file_size))
   else
     allocate(ADN(FINADN-INIADN+1))  
   endif  
   allocate(f(0:ALF1))
   
   
   
   call lee_adn(ADN,nbases,enes,file_adn,INIADN,FINADN,tipo,f)
   !**********************************************

   !***** INICIO DE LA LISTA  / ENES *************
   NCORTES=2
   call init_list(first,last,nbases,f)
   call INIT_RANDOM_SEED()
   call insert_enes(first,nbases,haygaps,enes)
  
   !***********************************************
   if (VERBOSE) then
     call write_info(0)
   endif 
   if (UN/=0) then
     call write_info(UN)
   endif
   !************************************************
  
   !call omp_set_num_threads(8)
  
     
   allocate(fcorte(0:ALF1))
   do j=1,nexper
     if (nexper/=1) then
       sl=slini+(j-1)*(slfin-slini)/(nexper-1) 
     else
       sl=slini
     endif    
     new_segs=1
     do while(new_segs/=0)
       new_segs=0
       eliminados=0
       call segment_list(first,S,nsegs) 
       !$OMP parallel do private(pcorte,fcorte,dmax)
       do i=1,nsegs
          call divmax(s(i),pcorte,fcorte,dmax)
          !calcula la significación y la compara con sl          
          if (sl_m_d(ALF,dmax,s(i)%fin-s(i)%ini+1)>sl) then
             call insert(s(i)%p,pcorte,fcorte)
             !$OMP critical
                NCORTES=NCORTES+1
             !$OMP end critical   
             new_segs=new_segs+1
          else
             s(i)%p%intact=1   
             eliminados=eliminados+1
          endif
       enddo
       !$OMP end parallel do
     enddo
     call segment_list(first,S,nsegs,.true.)
     totdiv=total_divergence(S,nsegs)
     if (VERBOSE) write(0,'(2f12.8,4i12)') sl,totdiv*sl,nsegs !,nsegs,NCORTES-nsegs
     
     !**** Print list of segment sizes *********
     !**** Only REAL segments, not n's *********
     nsegs_good=0
     do i=1,nsegs
       if (s(i)%p%intact/=-1) nsegs_good=nsegs_good+1
     enddo
     
     write(18,'(2f12.8,4i12)') sl,totdiv*sl,nsegs_good
     do i=1,nsegs
       if (s(i)%p%intact/=-1) then
         if (WRITE_FRECS) then
           write(18,'(5i14)') s(i)%fin-s(i)%ini+1,s(i)%p%next%frec(0:ALF1)-s(i)%p%frec(0:ALF1)
         else
           write(18,'(i14)') s(i)%fin-s(i)%ini+1
         endif  
       endif  
     enddo
     !******************************************
     !******************************************
     
     if (UN==0) then
        write(*,'(2f12.8,i12)') sl,totdiv*sl,nsegs
     else
        write(UN,'(2f12.8,i12)') sl,totdiv*sl,nsegs
     endif   
     call restore_list(first)
   enddo  
   
   !**********************************************
   
   close(18)
   
   
   
   
   

contains


!************************************************
!**** información de la secuencia ***************

subroutine write_info(UNIT)
   integer*4  unit
   
     write(unit,'(a)') '#_______________________________________'
     write(unit,'(a)') '#_______________________________________'
     write(unit,'(a)') '#'
     write(unit,'("# DNA file:          ",20(A))') trim(file_adn)
     write(unit,'("# Sequence length    ",i12)')   nbases
     if (INIADN/=1 .or. FINADN/=0) then
          write(unit,'("# Partial  from/to   ",2i12)') INIADN,FINADN
     endif
     write(unit,'("# Undefined          ",i12)')   enes
     if (ALF==2) then 
        write(unit,'("# Type               ",a2)')   tipo
     endif
     write(unit,'("# Frecs:             ",4i12)')  last%frec(0:ALF1)   
     write(unit,'("# Total Entropy      ",f12.8)')  entro(last%frec(0:ALF1))   
     write(unit,'(a)') '#'
     write(unit,'(a)') '#_______________________________________'
     write(unit,'(a)') '#_______________________________________'
     write(unit,'(a)') '#'
  
   
end subroutine write_info


!**************************************************
!*********** CALCULO DE DIVERGENCIA ***************
!**************************************************
FUNCTION total_divergence(S,nsegs)
  type(tseg) S(:)
  integer*8   i,nsegs
  integer*8   tamtot,tam
  integer*8  ,allocatable::ftot(:),f(:)
  real*8      diver,total_divergence
  
  allocate(ftot(0:ALF1))
  allocate(f(0:ALF1))
  diver=0.d0
  tamtot=0
  ftot=0
  do i=1,nsegs
    if (s(i)%p%intact/=-1) then
      tam=s(i)%fin-s(i)%ini+1
      tamtot=tamtot+tam
      f=s(i)%p%next%frec-s(i)%p%frec
      ftot=ftot+f
      diver=diver+entro(f)*tam
    endif
  enddo
  diver=entro(ftot)-diver/(1.d0*tamtot)
  total_divergence=diver
end function total_divergence  

  
  
  
!***************************************************
!***************************************************


!**************************************************
!**************** BUSCA CORTES ********************
!**************************************************
SUBROUTINE DiVMAX(s,pcorte,fcorte,dmax)
  
  type(tseg)  s
  integer*8   pcorte,fcorte(:)
  real*8      dmax
  
  integer*8  ,allocatable::f1(:),f2(:)
  integer*8   i,fin3
  real*8      d


  allocate(f1(0:ALF1))
  allocate(f2(0:ALF1))
  fin3=s%fin-3
  f1=0
  f2=s%p%next%frec-s%p%frec
  dmax=0.d0
  do i=s%ini,s%fin-1
    !if (adn(i).ne.-1) then   !**** no lo controlo porque no debe ocurrir
     f1(ADN(i))=f1(ADN(i))+1
     f2(ADN(i))=f2(ADN(i))-1	 
     d=jsdiver(f1,f2)
         !avoids segments shorter than 4 words
	 if ((d>=dmax) .and. (i-s%ini>3) .and. (i<fin3)) then
	   dmax=d
	   pcorte=i
	   fcorte=f1
	 endif
  enddo 
  fcorte=fcorte+s%p%frec
END SUBROUTINE DIVMAX
!*********************************************

!******************************************************************
!*************** OBTIENE LA LISTA DE SEGMENTOS ACTIVOS ************
!******************************************************************
subroutine segment_list(first,S,nsegs,full0)
   type(tseg)   ,allocatable::S(:)
   type(corte)  ,pointer::p,first
   integer*8     nsegs,i
   logical      ,optional::full0
   logical       full  
   ! full=.true. => gets all segments
   ! full=.false. => only those segmentable
   if (present(full0)) then
     full=full0
   else
     full=.false.
   endif    
   
   if (allocated(S)) deallocate(S)
   
   p=>first
   nsegs=0
   allocate(S(NCORTES))
   do while (associated(p%next))
     if (p%intact==0 .or. full) then   !1 => no se ha segmentado; -1 => son enes
       nsegs=nsegs+1
       S(nsegs)%ini=p%posic+1
       S(nsegs)%fin=p%next%posic
       S(nsegs)%p=>p
       !write(*,'(8i12)') p%posic+1,p%next%posic,p%next%posic-p%posic,p%next%frec-p%frec
     endif  
     p=>p%next
   enddo
end subroutine segment_list   



!***************************************************
!********** INSERTAR ENES **************************
!***************************************************
SUBROUTINE insert_enes(p0,size,haygaps,enes)
integer*8    size
type(corte),pointer::p,p0
integer*8    ini,ini_n,fin_n,hayenes,haygaps,i,j
integer*8    ftot(0:3),frec(0:3),enes

enes=0
HAY_RANDOM=.false.
ini=1
p=>p0
ftot=0
haygaps=0
do
   hayenes=findnext_n(size,ini,ini_n,fin_n)
   if (hayenes/=1) exit
   enes=enes+fin_n-ini_n
   frec=0
   do i=ini,ini_n
     if (ADN(i)/=-1) then
	   j=mod(adn(i),4)
       frec(j)=frec(j)+1
     endif
   enddo
   ftot=ftot+frec
   if (ini_n/=0) then
      call insert(p,ini_n,ftot) 
      NCORTES=NCORTES+1
      p=>p%next
   endif
   haygaps=1
   p%intact=-1           !contains n's and shouldn't be cut
   if (fin_n/=size) then
     call insert(p,fin_n,ftot)
     NCORTES=NCORTES+1
     p=>p%next
   endif
   ini=fin_n+1
enddo
if (HAY_RANDOM) then
  call ajusta_frecuencias(p0,enes)
endif
end SUBROUTINE insert_enes


!***************************************************

!***************************************************
!ini_seq = punto a partir de donde buscar
!ini     = inicio  de las enes - 1 (tal como se insertan los cortes)
!fin     = fin de las enes
FUNCTION findnext_n(size,ini_seq,ini,fin)
  integer*1  findnext_n,salida
  integer*8  LMIN
  integer*8  i0,i,ini,fin,size,suma,ini_seq
  real*8     r
  
  LMIN=MIN_N_SIZE
  i0=ini_seq
  777 continue
  ini=-1
  do i=i0,size
    if (ADN(i)==-1) then
      ini=i-1
      exit
    endif
  enddo

  if (ini/=-1) then
    salida=1
    do i=ini+1,size
      if (ADN(i)/=-1) then
        fin=i-1
        exit
      endif
      if (i==size) fin=size
    enddo
  else
    fin=0
    salida=0
  endif
  !si el trozo de n's tiene LMIN o menos se sustituyen por nucleótidos
  !aleatorios y se vuelve al principio a seguir buscando
  if (salida/=0.and.fin-ini<=LMIN.and.fin/=size) then
    HAY_RANDOM=.true.
    do i=ini+1,fin
      call random_number(r) 
      ADN(i)=r*4
    enddo
    i0=fin+1
    salida=0
    goto 777
  endif
  findnext_n=salida
end FUNCTION findnext_n
!********************************************
!como se pueden insertar trozos de hasta LMIN nucleotidos aleatorios
!hay que reajustar los vectores de frecuencias de los cortes
SUBROUTINE ajusta_frecuencias(first,enes)
  type(corte), pointer :: first,p
  integer*8  ini,fin,i,enes,f(0:3),j
 

  enes=0
  f=0
  p=>first
  do while (associated(p%next))
    ini=p%posic+1
    fin=p%next%posic
	if (p%intact==-1) then
	  enes=enes+(fin-ini+1)
    else 
	  do i=ini,fin
	    j=mod(ADN(i),4)
		f(j)=f(j)+1
	  enddo
    endif
    p=>p%next
	p%frec=f
  end do
END SUBROUTINE ajusta_frecuencias

!***************************************************



!***************************************************
!********** MANEJO DE LISTAS ***********************
!***************************************************
subroutine init_list(first,last,nwords,f_words)
  type(corte) ,pointer::first,last
  integer*8    nwords
  integer*8   ,dimension(:):: f_words

  !NCORTES=2  !variable global con el número de cortes
             !se incrementa cada vez que se inserta un corte  
             !no cuenta el último 

  allocate(first)
  allocate(last)
  allocate(first%frec(0:ALF1))
  allocate(last%frec(0:ALF1))

  first%posic=0
  first%frec=0
  first%next=>last
  first%intact=0
  nullify(first%prev)

  last%posic=nwords
  last%frec=f_words
  nullify(last%next)
  last%prev=>first
end subroutine init_list
!***********************************************
subroutine insert(p,posic,frec)
  type(corte), pointer :: p,newp
  integer*8               posic,frec(:)
  
  !NCORTES=NCORTES+1
  allocate(newp)
  allocate(newp%frec(0:ALF1))
  newp%posic=posic
  newp%frec=frec
  newp%next=>p%next
  newp%prev=>p
  newp%intact=0
  p%next%prev=>newp
  p%next=>newp
END SUBROUTINE insert


!**************************************************

!***************************************************

SUBROUTINE restore_list(first)
  type(corte),pointer :: first,p
  
  p=>first
  do while (associated(p%next))
    if (p%intact/=-1) p%intact=0  !Los segmentos de enes no se restauran
    p=>p%next
  enddo
end subroutine restore_list
  
  
!**************************************************



!***********************************************************
!********* LECTURA ADN BINARIO Y 4 BASES *******************
!***********************************************************
SUBROUTINE lee_adn(adn,size,enes,fen,ini,fin0,tipo0,f_words)
!************************************************************
character(*)   fen
integer*1      adn(:)
integer*8      size,letra,enes,ini,fin0,cont,n,i,f_words(0:)
integer*8      fin
character*83   linea
character*32   t_bases
character(2)   tipo,tipo0
logical        existe
integer*8      xletra


t_bases='ATCGNRYSWKMVHDBXatcgnryswkmvhdbx'
inquire(file=fen,exist=existe)
if (.not.existe) then
  size=-1
  goto 888
endif
open (101,file=fen)
if (fin0==0) then
   fin = 50400000000
else
   fin=fin0
endif
linea='kkk'
do while (linea(1:2) .ne. 'SQ' .and. linea(1:6) .ne. 'ORIGIN' .and. linea(1:1) .ne. '>')
   read (101,'(A83)') linea
enddo
size=0
enes=0
cont=0
f_words=0
tipo=tipo0
if (tipo(1:1)==tipo(2:2)) then
 xletra=index(t_bases,tipo(1:1))
 if (xletra>16) xletra=xletra-16
 tipo='XX'
endif
1  read(101,'(A83)',end=10) linea
   n=len_trim(linea)
   do i=1,n
      letra=index(t_bases,linea(i:i))
      if (letra>0) then
        cont=cont+1
        if (cont>fin) goto 10
        if (cont>=ini) then
          size=size+1
          if (letra>16) letra=letra-16
          if (letra>=5) then
             letra=5
             enes=enes+1
          endif
          if (letra>4) then
             letra=-1
          else
             if (tipo=='SW') then
        if ((letra==1).or.(letra==2)) then    ! AT =0, GC=1
                  letra=0
                else
                  letra=1
                endif
              else if (tipo=='RY') then
                if ((letra==1).or.(letra==4)) then
                  letra=0
                else
                  letra=1
                endif
              else if (tipo=='XX') then
                if (letra/=xletra) then
				   letra=0
                else 
				   letra=1
                endif
              else
                letra=letra-1
              endif
              f_words(letra)=f_words(letra)+1
          endif
          adn(size)=letra
        endif
	  endif
   enddo
   go to 1
10 continue
rewind(101)
close(101)

888 continue

end SUBROUTINE lee_adn

!****************************************************************

!***************** LEE PARAMETROS **************
 
 subroutine get_parameters(file_adn,slini,slfin,nexper,alf,tipo,file_pro,file_segment_sizes,INIADN,FINADN)
   integer*4      narg,i
   character*200  ss
 
   integer*8      INIADN,FINADN
   integer*8      alf 
   character*2    tipo
   integer*8      nexper
   character(*)   file_adn,file_pro,file_segment_sizes
   real*8         slini,slfin
   
   narg=iargc()
   if (narg<4) then
      write(0,*) 'ERROR: Required parameters missing'
      write(0,*) 
      call printhelp
      stop
   endif
   call getarg(1,file_adn)
   call getarg(2,ss)
   read(ss,*) slini
   call getarg(3,ss)
   read(ss,*) slfin
   call getarg(4,ss)
   read(ss,*) nexper
   
   !defaults
   alf=4
   tipo='4L'
   file_pro=''
   VERBOSE=.true.
   WRITE_FRECS=.false.
   file_segment_sizes='segment_sizes.dat'
   INIADN=1
   FINADN=0
   
   if (narg==4) return
   i=5  
   do 
     call getarg(i,ss)
     select case (trim(ss))
     case ('-a','-A')  
       i=i+1
       call getarg(i,ss)
       read(ss,*) alf
       if (alf==2) then
         i=i+1
         call getarg(i,tipo)
       else
         tipo='4L'  
       endif
     case ('-part','-PART')
       i=i+1
       call getarg(i,ss)
       read(ss,*) INIADN
       i=i+1
       call getarg(i,ss)
       read(ss,*) FINADN
     case ('-h','-H')
       call printhelp
       stop  
     case ('-o','-O')
       i=i+1
       call getarg(i,file_pro)
     case ('-s','-S')
       i=i+1
       call getarg(i,file_segment_sizes) 
     case ('+f','+F')
        WRITE_FRECS=.true.   
     case ('-v','-V')
        VERBOSE=.false.  
     endselect   
     i=i+1
     if (i>narg) exit  
   enddo
end subroutine get_parameters

 !**************************************************

 !******************** help ******************
 subroutine printhelp
   character*100 s   
   call getarg(0,s)
   write(0,*)
   write(0,*) 'Calculates the SCC profile of a DNA sequence'
   write(0,*) '(80 column GeneBank or FASTA file)'
   write(0,*)
   write(0,*) 'USAGE:'
   write(0,*) '    '//trim(s)//' <DNA file> <init. sl> <end sl> <# exp.> [options]'
   write(0,*)       
   write(0,*) 'OPTIONS: '      
   write(0,*) '-a  ALPHABET [type]    ALPHABET=4 => 4 bases'
   write(0,*) '                       ALPHABET=2 => binary, must include type: RY,SW or KM'
   write(0,*) 
   write(0,*) '-part <start> <end>    Reads DNA sequence only from start to end'
   write(0,*) '                       end=0 => reads from "start" to the end of the sequence'
   write(0,*) 
   write(0,*) '-o <File>              Output file with the SCC profile (Default: Std. output)'
   write(0,*) 
   write(0,*) '-s <File>              Output file with segment sizes (Default: "segment_sizes.dat")'
   write(0,*)
   write(0,*) '+f/-f                  Include frequencies in output (Default: -f = NO)'
   write(0,*) 
   write(0,*) '-v                     No verbose. Do not send extra info to Std. error'
   !write(0,*)
   !write(0,*) '-h              Prints this help screen and exit'
   write(0,*)
 end subroutine printhelp
 !**********************************************



!***********************************************

      !!!!!!!!!!!!!! BUILT IN RANDOM GENERATOR  !!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE INIT_RANDOM_SEED()
      INTEGER::I,N,CLOCK
      INTEGER,DIMENSION(:),ALLOCATABLE::SEED
          
      CALL RANDOM_SEED(SIZE=N)
      ALLOCATE(SEED(N))
      CALL SYSTEM_CLOCK(COUNT=CLOCK)
      SEED=CLOCK+37*(/ (I-1, I=1, N) /)
      CALL RANDOM_SEED(PUT=SEED)
      DEALLOCATE(SEED)
      RETURN
      END SUBROUTINE INIT_RANDOM_SEED

      !!! CALL RANDOM_NUMBER(x)   X real*8

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***************************************************

end
