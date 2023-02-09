      program sizeofint
      integer  p,i
      p=1024*1024
      i= p*p
      if (i>0) then
         call abort()
      else
         print *, 4
      endif
      end
