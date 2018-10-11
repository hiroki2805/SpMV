subroutine heapsort(n,array)
  implicit none
  integer,intent(in) :: n
  integer,intent(inout) :: array(1:n)

  integer ::i,k,j,l
  integer :: t

  if(n.le.0)then
     write(6,*)"Error, at heapsort"; stop
  endif
  if(n.eq.1)return

  l=n/2+1
  k=n
  do while(k.ne.1)
     if(l.gt.1)then
        l=l-1
        t=array(l)
     else
        t=array(k)
        array(k)=array(1)
        k=k-1
        if(k.eq.1) then
           array(1)=t
           exit
        endif
     endif
     i=l
     j=l+l
     do while(j.le.k)
        if(j.lt.k)then
           if(array(j).lt.array(j+1))j=j+1
        endif
        if (t.lt.array(j))then
           array(i)=array(j)
           i=j
           j=j+j
        else
           j=k+1
        endif
     enddo
     array(i)=t
  enddo

  return
end subroutine heapsort
