      subroutine velfactors(vxm,vym,vzm,vxa,vya,vza,
     1         vm,va,gammam,gammaa,uvxm,uvym,uvzm,uvxa,uvya,uvza)
c
      real*8 vxm,vym,vzm,vxa,vya,vza
      real*8 vm,va,gammam,gammaa
      real*8 uvxm,uvym,uvzm,uvxa,uvya,uvza
c
c speed and gamma (boost factor) are the same for m and mbar.
      vm=sqrt(vxm**2+vym**2+vzm**2)
      gammam=1./sqrt(1.-vm**2)
      va=sqrt(vxa**2+vya**2+vza**2)
      gammaa=1./sqrt(1.-va**2)
      print*, ' monopole velocity, boost ', vm, gammam
c
c unit velocity vectors:
        if(vm.le.0.00001) then
         uvxm=0.
         uvym=0.
         uvzm=0.
         uvxa=0.
         uvya=0.
         uvza=0.
        else
         uvxm=vxm/vm
         uvym=vym/vm
         uvzm=vzm/vm
         uvxa=vxa/va
         uvya=vya/va
         uvza=vza/va
        endif
c
      return
      end
