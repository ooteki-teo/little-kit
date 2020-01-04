program extractband
    implicit none
    integer :: npara
    character(2) :: darg
    character(10) :: dargerange
    integer :: i, j, k, l, m, n, iost, ierr
    real :: ko(3), kn(3), dtmp, dsc,ef,rko(3),rkn(3)   
    real :: b1(3),b2(3), b3(3)
    real, allocatable :: kdist(:)
    integer :: nk, nb, nline
    character(1) :: customrange
    real, allocatable :: d(:), reslt(:, :)
    character(1) :: kname(20)
    real :: emin, emax
    integer :: nkline
    integer :: nkinline
    integer :: innk
    real, allocatable :: hklist(:,:), klist(:)
    real :: dk(3), dist, kold(3), knew(3),rdk(3)

    kname="K"
    
    npara=iargc()
    if (npara.ge.1) then
        do i=1,npara
        call getarg(i,darg)
        if (darg.eq."-f") then
            customrange="f"
        endif    
        if (darg.eq."-y") then
            customrange="c"
            call getarg(i+1,dargerange)
            read(unit=dargerange, fmt=*)emin
            call getarg(i+2,dargerange)
            read(unit=dargerange, fmt=*)emax
        endif
        if (darg.eq."-k") then
            call getarg(i+1, darg)
            read(unit=darg, fmt=*) innk
            do j=2,innk+1
                call getarg(i+j, kname(j-1))
            enddo
        endif
        if (darg.eq."-h") then
            write(*,*) "band.o accepted the following arguments:"
            write(*,*) "-h"
            write(*,*) "    help messages"
            write(*,*) "-f"
            write(*,*) "    plot E-fermi +- 2.eV range"
            write(*,*) "-y a1 a2"
            write(*,*) "    plot Eenergy range by custom [a1, a2]"
            write(*,*) "-k nk k1 k2 k3 ..."
            write(*,*) "    specified k point lable, character_len=1"
            stop
        endif
        enddo
    endif    

    ko=0d0
    kn=0d0
    dtmp=0d0
    ef=0d0

!-----------------------
! kpt
    
    dk=0d0
    dist=0d0

    call system("awk '{if (NR>4&&NF>=3) {print $1,$2,$3} }' KPOINTS > kline.dat")
    call system("awk 'END{print NR}' kline.dat > nk.dat")
    call system("grep -A 3 'reciprocal lattice vectors' OUTCAR | tail -n 3 |awk '{print $4,$5,$6}' > reci.dat")

    !> reciprocal lattice vectors
    open(unit=121,file="reci.dat")
    read(121,*)b1(:) 
    read(121,*)b2(:) 
    read(121,*)b3(:) 
    close(121)
    
    !> number of k lines
    open(unit=122, file="nk.dat")
    read(122,*)nkline
    nkline=nkline/2+1
    close(122)
    call system("rm nk.dat")

    allocate(kdist(nkline+1))
    allocate(hklist(nkline+1, 3))
    kdist=0d0
    hklist=0d0

    !> read high symmetry k points list
    open(unit=123, file="kline.dat")
    read(123,*)hklist(1,:)
    read(123,*)hklist(2,:)
    do i=3,nkline
        read(123,*)
        read(123,*) hklist(i,:)
    enddo
    close(123)
    
    do i=1, nkline
        write(*,*) "# kpt:",i
        write(*,*) hklist(i,:)
    enddo
    !> generate kpoints distance list
    do i=2,nkline+1
        dk=hklist(i,:)-hklist(i-1,:)
        call c2d(dk,rdk,b1,b2,b3)
        write(*,*) "# kpt:",i
        write(*,*) dk
        write(*,*) rdk
        dist=sqrt(rdk(1)**2+rdk(2)**2+rdk(3)**2)
        kdist(i)=kdist(i-1)+dist
    enddo
    
    

!----------------------
    
    
    open(unit=109, file="EIGENVAL", status="old")
    open(unit=110, file="ek.dat")

    !> fermi energy
    call system("grep E-fermi OUTCAR| awk '{print $3}' > fermienergy.dat")
    open(unit=111, file='fermienergy.dat')
    read(111,*) ef
    close(111)
    call system("rm fermienergy.dat")
    
    do i=1,5
        read(109, *) 
    end do

    read(109, *) i, nk, nb
    read(109, *)
    
    allocate(reslt(nk, nb), d(nk))
    reslt=0d0
    d=0d0
   
    do i=1,nk
        read(109, *) kn(1), kn(2), kn(3), dsc
        if (i.eq.1) ko=kn
        do j=1,nb
            read(109, *) m, reslt(i, j), dsc
        enddo
        call c2d(kn, rkn, b1,b2,b3)
        call c2d(ko, rko, b1,b2,b3)
        d(i)=dtmp+sqrt((rkn(1)-rko(1))**2 + (rkn(2)-rko(2))**2 + (rkn(3)-rko(3))**2 )   
        dtmp=d(i)
        ko = kn
        if (i.eq.nk) cycle
        read(109, *)
        
    enddo

    do i=1,nb   
        do j=1, nk
            write(110, '(f15.7, f15.7)')  d(j), reslt(j, i)-ef
        end do
        write(110, *) "  "
    end do
    close(109)
    close(110)
204 format('set arrow from ',F10.5,',',F10.5,' to ',F10.5,',',F10.5, ' nohead')
203 format(A3,'" ',F10.5,')')
202 format('set xtics (',20('"',A3,'" ',F10.5,','))
    if (customrange=='f') then
        emin=-2.00
        emax=2.00
    else if(customrange=="c") then
        write(*,*) "using custom range"
    else
        emin=minval(reslt)-ef
        emax=maxval(reslt)-ef
    endif
    
    open(111, file="ek.gnu")
    write(111, '(a)') 'set encoding iso_8859_1'
    write(111, '(a)') 'set size square '
    write(111, '(a)') 'set terminal png truecolor enhanced size 1920, 1680 font ",60"'
    write(111, '(a)') 'set output "ek.png"'
    write(111, '(a)') 'unset ztics'
    write(111, '(a)')'unset key'
    write(111, *)'emin=', emin
    write(111, *)'emax=', emax
    write(111, '(a)')''
    write(111, '(a)') 'set ylabel "Energy (eV)"'
    write(111, *) 'set xrange[',minval(d),':',maxval(d),']'
    write(111,202, advance="no") (kname(i), kdist(i), i=1, nkline-1)
    write(111, 203)kname(nkline),kdist(nkline)
    write(111,*)" "
    write(111, 204)minval(d), 0 ,maxval(d),0
    write(111, *)'set yrange[ emin:emax]'
    do i=1, nkline
        write(111, 204) kdist(i), emin, kdist(i), emax
    enddo
    write(111, '(a)') 'set title "BandStructure"'
    write(111, *) 'plot "ek.dat"',' u 1:2 w l lt 2 lw 2 lc rgb "#ff0000"'
    close(111)
    
    write(*,*)"type: gnuplot ek.gnu"
    call system("gnuplot ek.gnu")
     
end program extractband

subroutine c2d(b,a,r1,r2,r3)
    implicit none
    real,intent(in)::b(3), r1(3), r2(3), r3(3)
    real,intent(out)::a(3)
    
    a(1)=b(1)*r1(1)+b(2)*r2(1)+b(3)*r3(1)
    a(2)=b(1)*r1(2)+b(2)*r2(2)+b(3)*r3(2)
    a(3)=b(1)*r1(3)+b(2)*r2(3)+b(3)*r3(3)
    
end subroutine





