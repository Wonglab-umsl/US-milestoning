real*8 function logsum(a,b)
implicit none
	real*8 a,b
        if (a>b) then
                logsum=a+log(1+exp(b-a))
        else
            	logsum=b+log(1+exp(a-b))
        endif
	return
end function logsum


module stuff
implicit none
integer nrep,ndatamax,histdim(1:3),maxiter

real*8 res, beta, qmin(1:3),qmax(1:3),whamtol,mininf,timestart,timeend
parameter(mininf=-1000.0d0)
!1st index -- data, 2nd index -- simulation
!kbias -- biasing force constant, q0 -- center of harmonic potential
!q -- com data, etc.
real*8,allocatable::kbias(:),q0(:,:),q(:,:,:),lnz(:),zz(:),denom(:,:,:)
!q1bin,q2bin,ubin,simulation
integer,allocatable::hist(:,:,:,:),histtot(:,:,:),ndata(:)


contains

subroutine readcontrolfile(fname)
implicit none
character*255 fname
integer i,j,k,iostat,ndatatot
character*255 fname2
real*8 time,qq(1:3)
	open(unit=1,file=fname,status="old")
	read(1,*) nrep,ndatamax,timestart,timeend,res,beta,maxiter,whamtol
	allocate(kbias(1:nrep),q0(1:3,1:nrep),q(1:3,1:ndatamax,1:nrep),ndata(1:nrep))
	ndatatot=0
	do i=1,nrep
		read(1,*) kbias(i),q0(1,i),q0(2,i),q0(3,i)
		read(1,'(A)') fname2
		open(unit=2,file=fname2,status="old")
		write(*,*) "reading file ",fname2
		j=1
		do while (.true.)
			if (j>ndatamax) then
				write(*,*) "too much data on simulation ",i,j,ndatamax
				stop
			endif
			read(2,*,iostat=iostat) time,(qq(k),k=1,3)
			if (iostat/=0) exit
			if ((time>=timestart) .and. (time<timeend)) then
				q(:,j,i)=qq(:)
				j=j+1
			endif
		end do
		ndata(i)=j-1
		write(*,*) ndata(i)," data points read."
		close(unit=2)
		ndatatot=ndatatot+ndata(i)
	end do
	close(unit=1)
	write(*,*) "total ",ndatatot," data points read."
end subroutine readcontrolfile

subroutine createhistograms
implicit none
integer i,j,k,qbin(1:3),nbintot,ix,iy,iz
	!qmin(:) = minval(minval(q(:,:,:),3),2)
	!qmax(:) = maxval(maxval(q(:,:,:),3),2)
	do k=1,3
		qmin(k)=1000
		qmax(k)=-1000
		do i=1,nrep
			do j=1,ndata(i)
				if (q(k,j,i)>qmax(k)) qmax(k)=q(k,j,i)
				if (q(k,j,i)<qmin(k)) qmin(k)=q(k,j,i)
			end do
		end do	
		qmin(k)=res*floor(qmin(k)/res)
		qmax(k)=res*ceiling(qmax(k)/res)
		histdim(k) = int((qmax(k)-qmin(k))/res)
	end do
	allocate(hist(1:histdim(1),1:histdim(2),1:histdim(3),1:nrep))
	allocate(histtot(1:histdim(1),1:histdim(2),1:histdim(3)))
	!allocate(hist(1:nubin,1:nrep),histtot(1:nubin),histtot2(1:nubin,1:nq1bin,1:nq2bin),histtot3(1:nq1bin,1:nq2bin))
	hist=0
	do i=1,nrep
		write(*,*) "creating histogram for replica ",i
		do j=1,ndata(i)
			do k=1,3
				qbin(k) = int((q(k,j,i)-qmin(k))/res)+1
			end do
			hist(qbin(1),qbin(2),qbin(3),i) = hist(qbin(1),qbin(2),qbin(3),i) + 1
		end do
	end do
	histtot(:,:,:) = sum(hist(:,:,:,:),4)
	nbintot=0
	do ix=1,histdim(1)
	do iy=1,histdim(2)
	do iz=1,histdim(3)
		if (histtot(ix,iy,iz)>0) nbintot=nbintot+1
	end do
	end do
	end do 	
	write(*,*) "total ",nbintot," histogram bins have at least one data point"
	write(*,*) "a volume of ",real(nbintot)*res*res*res," A^3"
end subroutine createhistograms

subroutine dowham
implicit none
integer ix,iy,iz,j,k,iter
real*8 conv,ebias,uu,qcenter(1:3),qdiff(1:3)
real*8 logsum
	allocate(lnz(1:nrep),zz(1:nrep),denom(1:histdim(1),1:histdim(2),1:histdim(3)))
	lnz=0.0
	iter=0
	write(*,*) "beginning wham"
	do while (.true.)
		!denominator, eq. 7.3.10
		denom=mininf		
		do ix=1,histdim(1)
		do iy=1,histdim(2)
		do iz=1,histdim(3)
			if (histtot(ix,iy,iz)>0) then
				qcenter(1)=(ix-0.5)*res+qmin(1)
				qcenter(2)=(iy-0.5)*res+qmin(2)
				qcenter(3)=(iz-0.5)*res+qmin(3)
				do j=1,nrep
					qdiff(:)=qcenter(:)-q0(:,j)
					uu=0.5*kbias(j)*dot_product(qdiff,qdiff)
					denom(ix,iy,iz)=logsum(denom(ix,iy,iz),-beta*uu+log(real(ndata(j)))-lnz(j))
				end do
			endif
		end do
		end do
		end do

		!wham equation 7.3.10
		zz=mininf
		do j=1,nrep
			do ix=1,histdim(1)
                	do iy=1,histdim(2)
                	do iz=1,histdim(3)
				if (histtot(ix,iy,iz)>0) then
                        		qcenter(1)=(ix-0.5)*res+qmin(1)
                        		qcenter(2)=(iy-0.5)*res+qmin(2)
                       			qcenter(3)=(iz-0.5)*res+qmin(3)
					qdiff(:)=qcenter(:)-q0(:,j)
                                	uu=0.5*kbias(j)*dot_product(qdiff,qdiff)
					zz(j)=logsum(zz(j),-beta*uu+log(real(histtot(ix,iy,iz)))-denom(ix,iy,iz))
				endif
			end do
			end do
			end do
		end do
		!reset z(1)=0
		do j=2,nrep
			zz(j)=zz(j)-zz(1)
		end do	
		zz(1)=0.0
		conv=0.0
		do j=1,nrep
			conv=max(conv,abs(zz(j)-lnz(j)))
		end do
		iter=iter+1
		write(*,*) iter,conv
		lnz=zz
		write(*,*) lnz
		if ((conv<whamtol).or.(iter>maxiter)) exit
	end do
end subroutine dowham

!subroutine writepmf(outfile)
!implicit none
!character*255 outfile
!integer ix,iy,iz,k
!real*8 q(1:3)
!	open(unit=2,file=outfile,status="unknown")
!	do ix=1,histdim(1)
!	do iy=1,histdim(2)
!	do iz=1,histdim(3)
!		if (histtot(ix,iy,iz)>0) then
!			q(1)=(ix-0.5)*res+qmin(1)
!			q(2)=(iy-0.5)*res+qmin(2)
!			q(3)=(iz-0.5)*res+qmin(3)
!			!equation 7.3.9
!			write(2,'(3F12.4,I5,F12.6)') q(1),q(2),q(3),histtot(ix,iy,iz),&
!				-(1/beta)*(log(real(histtot(ix,iy,iz)))-denom(ix,iy,iz))
!		endif
!	end do
!	end do
!	end do
!end subroutine writepmf

!write in the DX format supported by VMD
subroutine writepmf_dx(outfile,dxfile)
implicit none
character*255 outfile,dxfile
integer ix,iy,iz,k,counter,ntot
real*8 q(1:3),pmfval,minpmf,d
real*8, allocatable::pmf(:,:,:),flatpmf(:),dist(:)
	ntot=histdim(1)*histdim(2)*histdim(3)
	allocate(pmf(1:histdim(1),1:histdim(2),1:histdim(3)),flatpmf(1:ntot),dist(1:nrep))
        open(unit=2,file=dxfile,status="unknown")
	open(unit=3,file=outfile,status="unknown")
	write(2,'("object 1 class gridpositions counts ",3I4)') histdim(1),histdim(2),histdim(3)
	write(2,'("origin ",3(F8.4,1X))') qmin(1),qmin(2),qmin(3)
	write(2,'("delta ",3F8.4)') res,0.0,0.0
	write(2,'("delta ",3F8.4)') 0.0,res,0.0
	write(2,'("delta ",3F8.4)') 0.0,0.0,res
	write(2,'("object 2 class gridconnections counts ",3I4)') histdim(1),histdim(2),histdim(3)
	write(2,'("object 3 classs array type double rank 0 items",I9," data follows")') ntot
        do ix=1,histdim(1)
        do iy=1,histdim(2)
        do iz=1,histdim(3)
                if (histtot(ix,iy,iz)>0) then
                        !q(1)=(ix-0.5)*res+qmin(1)
                        !q(2)=(iy-0.5)*res+qmin(2)
                        !q(3)=(iz-0.5)*res+qmin(3)
                        !equation 7.3.9
			pmfval=-(1/beta)*(log(real(histtot(ix,iy,iz)))-denom(ix,iy,iz))
			!write(3,'(3F12.4,I8,F12.6)') q(1),q(2),q(3),histtot(ix,iy,iz),pmfval
		else
			pmfval=1000
		endif
		pmf(ix,iy,iz)=pmfval
        end do
        end do
        end do
	!ensure minimum of pmf is zero
	minpmf=minval(pmf)
	counter=1
	do ix=1,histdim(1)
        do iy=1,histdim(2)
        do iz=1,histdim(3)
                if (histtot(ix,iy,iz)>0) then
                        q(1)=(ix-0.5)*res+qmin(1)
                        q(2)=(iy-0.5)*res+qmin(2)
                        q(3)=(iz-0.5)*res+qmin(3)
			d=sqrt(dot_product(q(:)-q0(:,1),q(:)-q0(:,1)))
			write(3,'(5F12.4,I8,F12.6)') q(1),q(2),q(3),interp_position(q,dist),d,histtot(ix,iy,iz),pmf(ix,iy,iz)-minpmf
		endif
		flatpmf(counter)=pmf(ix,iy,iz)-minpmf
		counter=counter+1
        end do
        end do
        end do
	close(unit=3)
	write(2,'(3F20.10)') (flatpmf(k),k=1,ntot)
	!write(2,*) 'attribute "dep" string "positions"'
	!write(2,*) 'object "regular positions regular connections"'
	!write(2,*) 'component "positions" value 1'
	!write(2,*) 'component "connections" value 2'
	!write(2,*) 'component "data" value 3'
	write(2,*) 'object "data" class field'
	close(unit=2)
	deallocate(pmf)
end subroutine writepmf_dx

real*8 function interp_position(q,dist)
implicit none
real*8 q(1:3),dist(:)
integer irep,minrep,ireplo,irephi
real*8 qdiff(1:3),vec(1:3),vec2(1:3),frac
	do irep=1,nrep
		qdiff(:)=q(:)-q0(:,irep)
		dist(irep)=sqrt(dot_product(qdiff,qdiff))
	end do
	minrep=minloc(dist,1)
	if (minrep==1) then
		ireplo=minrep
		irephi=minrep+1
	else if (minrep==nrep) then
		ireplo=minrep-1
		irephi=minrep
	else if (dist(minrep-1)<dist(minrep+1)) then
		ireplo=minrep-1
		irephi=minrep
	else 
		ireplo=minrep
		irephi=minrep+1
	end if
	!project q onto the the line connecting q0(ireplo) and q0(irephi) and express as a fraction of the distance
	vec(:)=q(:)-q0(:,ireplo)
	vec2(:)=q0(:,irephi)-q0(:,ireplo)
	frac=dot_product(vec,vec2)/dot_product(vec2,vec2)
	interp_position=real(ireplo)+frac
end function interp_position

end module stuff

program wham3d
use stuff
implicit none
character*255 fname,dxfname
	call getarg(1,fname)
	call readcontrolfile(fname)
	call createhistograms
	call dowham
	call getarg(2,fname)
	call getarg(3,dxfname)
	call writepmf_dx(fname,dxfname)
	deallocate(kbias,q0,q,lnz,zz,hist,histtot,ndata)
end program wham3d
