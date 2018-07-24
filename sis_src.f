c Yamila Miguel PhD thesis, Marzo 2011
c Referencias:
c Miguel, Y., Guilera, O. & Brunini, A., 2011b, MNRAS, 417, 314
c Miguel, Y., Guilera, O. & Brunini, A., 2011a, MNRAS, 412, 2113
c Miguel, Y. & Brunini, A., 2010, MNRAS, 403, 1935
c Miguel, Y. & Brunini, A., 2009, MNRAS, 392, 324
c Miguel, Y. & Brunini, A., 2008, MNRAS, 387, 463
c este programa contiene:
c acrecion de embriones en crecimiento oligarquico
c inicialmente espaciados en 10 radios de Hill.
c los planetesimales migran por el drag con el gas nebular
c acrecion gaseosa con la contraccion K-H
c embriones separados en < 3.5 radios de Hill se fusionan en uno solo
c pero pierden todo el gas que podrian haber acretado hasta el momento 
c aumento de la seccion eficaz de captura por friccin con la
c atmosfera planetaria, si la hay.
c los planetas pueden migrar  migracion con de tipo I y II.
c el gas se disipa exponencialmente.
c el perfil de disco considerado es el que obtuvimos de acuerdo con
c Andrews et al. 2009 y Isella et al. 2009
c----------------------------------------------------------------      
      implicit real*8(a-h,o-z)
      real*8 metal,qest,lsmin,lsmax,lsmmin,lsmmax,lsm,lsdt
      integer disipa,co,na,nn, ncol(100000)
      dimension emepla(100000),a(100000),erreh(100000),wcrit(100000)
      dimension sigmap(100000),indgas(100000),realo(100000)
      dimension emeiso(100000),emedot_crit(100000),w(100000) 
      dimension emegas(100000),rplanet(100000),prot(100000),ele(3)
      dimension in(100000),eleang(100000,3),ele1(3),ele2(3),obli(100000)
      dimension bb(100000),ele3(3)
      dimension lsac(100000)
      logical verb
      Character*30 ArchivoSalida 
      Character*30 SalidaCorta ! Agregada 
      Character*30 Gigantes 
      Character*30 ident
c----------------------------------------------------------------
      common/pasar/ro_ice,r_med,g,dt,pi,emesol,uacm,yearsec
     1     ,det,c1
      common/pla/emepla,b,erreh,s,amin,emestar
      common/others/apert,fpert,constmigI     
cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
c     CANTIDADES VARIABLES EN CADA SIMULACION
cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      alfa=1.d-3                ! alfa de Shakura-Sunyaev (parametro adimesional que caracteriza a la viscosidad)
      emecrit= 10.d0            ! masa critica para la acrecion de gas
c     disipa= 0                 ! si disipa= 0, cuando chocan 2 embriones, la atmosfera
c     ! se conserva. si disipa= 1 la atmosfera se pierde   
cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
cmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      r_min=1.d4                ! radio del planetesimal mas pequeno
      r_max=1.d7                ! radio del planetesimal mas grande
      wra=2.5                   ! valor absoluto del indice de distribucion de radios debe ser >2    
      ro_ice=1.5d0              ! densidad de los planetesimales mas alla de la snow-line
      ro_ro= 3.d0               ! idem pero para los que estan mas cerca de la snow-line
      b=10.d0                   ! separacion entre embriones oligarcas
      superpos= 3.5d0           ! superposicion de la zona de alimentacion
c----------------------------------------------------------------
      nsal=100                  ! espaciado nsal*dt de las salidas
      tini=1.                   ! tiempo de inicio de la simulacion
      tfin=$Tfin                ! tiempo de la simulacion   
      s=0.05d0                  ! escala de altura del disco
      vluz=3.d10                ! velocidad de la luz(cm/seg)
      cappa=4.d0                ! opacidad media de la atmosfera
      pesomol=2.8d0             ! peso molecular del gas de la atmosfera
      dt=100.d0                 ! paso de tiempo en anos
      efe=0.01d0                ! radio del planeta en terminos del radio de Hill
c----------------------------------------------------------------
      nplanetas=1               ! en realidad es el numero de sistemas planetarios!
      verb=$Verb
c     parametros para la perturbacion
      apert=$Apert
      fpert=$Fpert
c     lee las variables del sistema
      emed=$Emed
      rc=$Rcar
      emestar=$Emestar
      qest=$Qest
      emetal=$Emetal
      taugas=$Taugas
      ArchivoSalida=$ArchivoSalida
      SalidaCorta=$SalidaCorta
      Gigantes=$Gigantes
      ident="$Identifier"
c     **********************************
c     defino como es el perfil del disco, vamos a probar con 0.5, 1 y 1.5
      gama=$Gama
c     **********************************
      constmigI=$ConstMigI
c----------------------------------------------------------------
      if(verb)then
         open(333,file=ArchivoSalida) !archivo de salida. 
      endif
      open(444,file=SalidaCorta) !archivo de salida. 
      open(555,file=Gigantes) !archivo de salida.
c----------------------------------------------------------------
c     cantidades auxiliares y ctes:
      pi=4.d0*datan(1.d0)
      uacm=1.5d13
      yearsec=365.d0*86400.d0
      emesol=2.d33
      emet=5.97d27              !masa de la tierra en gramos, valor correcto
      radtie=6.37814d8
      g=6.672d-8                !en cgs
      ctpi=4.18879d0            ! =(4/3)*pi
c     calculo del radio medio de los planetesimales
      eme_min=4.d0*pi*ro_ice*r_min**3/3.d0
      eme_max=4.d0*pi*ro_ice*r_max**3/3.d0
      x1=4.d0-wra    !radio minimo de la distribucion de planetesimales
      x2=wra-1.d0    !radio maximo de la distribucion de planetesimales
      eme_med=(x2/x1)*((eme_max**x1)*(eme_min**x2))**0.333333 ! 
      r_med=(3.d0*eme_med/4.d0*pi*ro_ice)**0.3333333
      acre_gas=7.5d-4*pesomol**4*vluz
      enumber=dexp(1.d0)
c     parametros usados para estudiar las colisiones
      cont=0.0
      dseed=145301987.d0 !numero semilla 


      do ii=1, nplanetas
c----------------------------------------------------------------
c     inicializo los pasos del indice del tiempo en logspace
c     para el output, con lsnr = no. de snapshots
c     lsnr,lsi,lsj,lsk,lsnew,lsold son int, los demas son reales
c----------------------------------------------------------------
	  lsnr = 500
	  lsmin=tini
	  lsmax=tfin
	  lsmmin=log10(lsmin)
	  lsmmax=log10(lsmax)
	  lsm=lsmmin
	  lsdt=(lsmmax-lsmmin)/(lsnr-1.)
	  pasost=tfin/dt
	  if(lsnr.ge.int(pasost/10.))then
		  print *,'Atencion: lsnr debe ser menor que npasost/10'
		  stop
	  endif
	  lsidx=1
  	  do lsi=1,lsnr+1
  		lsnew=10.**lsm
  		if(lsi.ge.2)then
  			if(lsnew.ne.lsold)then
  				lsac(lsidx)=lsold
  				lsidx=lsidx+1				
  			endif
  		endif
  		lsm=lsm+lsdt
  		lsold=lsnew
  		enddo
	  lsk=1
c-----------------------------------------------------------------
c     generacion de los planetas iniciales separados b * Rh
c     y del disco inicial
c-----------------------------------------------------------------



c     inicializo algunos valores
         n=1
         t=tini
c     lee la escala de tiempo de disipacion de la nebulosa

c     lee el Rc que nos va a delimitar el disco
         if(emestar.ge.0.5d0)then 
            elestar= (emestar)**4.0 !luminosidad de la estrella en
c     funcion de la del Sol   
         else
            elestar= (emestar)**2.8d0
         endif

         amin=0.07d0*sqrt(elestar) !limite interno del disco
         a(1)=amin
         rice= 2.7d0*emestar**2 ! linea de condensacion de los hielos
         emestar=emesol*emestar !masa de la estrella en gramos
c     factores numericos varios usados mas adelante  
         faktor=3.9*b**(2./5.)*dsqrt(g)*emestar**(1./6.)
         faktor=faktor/ro_ice**(9./15.)
         faktor=faktor/eme_med**(2./15.)
c     ----------------------
c     densidad de gas (en realidad la calculo mas adelante, porque aca no la necesito,
c     solo calculo el sigmag_0 para poder hallar el sigma de solidos
         sigmag_0=emed*(2.0d0-gama)*emesol/(2.*pi*rc**2.*2.25d26)!esta en g/cm2
c     -------------------------------------------------------------------
         sigmad_0=0.0149*sigmag_0*10.**emetal !esta en gr/cm2
c     --------------------
c     calculo de amax, que es el a que contiene el 95% de la masa del disco
         amax=3.0d0**(1.0d0/(2.0d0-gama))*rc !en las mismas unidades que rc, en UAs
c     --------------------------------------


 100     if(a(n).gt.rice)then
            fice=1.d0
            ropla=ro_ice        !densidad del planeta
         else
            fice=0.25d0
            ropla=ro_ro
         endif

c     -------------------
c     densidad de solidos
c     aqui entra el bloque de perturbacion
c     -------------------

         sigmap(n)=sigmad_0*fice*(a(n)/rc)**(-gama)*
     &        exp(-(a(n)/rc)**(2.-gama))*
     &     (1+apert*cos(2.*pi*a(n)/(fpert*0.03361386903*a(n)**(5./4.)))) !densidad de solidos en g/cm2
c	  print *, sigmad_0,fice,a(n),rc,-gama, sigmap(n)

         a_cm=a(n)*uacm         !grilla en cm
         emepla(n)=eme_med*
     &        sigmap(n)*pi*a_cm**2*b/(1.5d0*emestar)**0.333333
         emepla(n)=emepla(n)**(3.d0/5.d0) !masa inicial del planeta, 
c     esta es la necesaria para que comience el OG (Thommes et al-Ida & Makino)
         rplanet(n)=(emepla(n)/(ctpi*ropla))**0.3333333 !radio del planeta en cm
         emegas(n)=0.d0         !inicializo la masa de gas inicial del planeta a 0
c     momentos angulares iniciales de los planetas (inicialmente no rotan)
         eleang(n,1)=0.d0!Lx componente x del planeta n
         eleang(n,2)=0.d0!Ly del planeta n
         eleang(n,3)=0.d0!Lz del planeta n
         obli(n)=90.d0
         realo(n)=0.d0
         ncol(n)=0       !este es el numero de colisiones sufridas, inicialmente es 0
         erreh(n)=a(n)*(emepla(n)/(3.d0*emestar))**0.33333333!radio de hill del planeta en UA
         da=a(n)*b*(0.667d0*emepla(n)/emestar)**0.33333333 !!la separacion entre los embriones es de b*rh, pero el rh de las dos masas por eso el 2/3
         indgas(n)=0 
         a(n+1)=a(n)+da
         n=n+1
         if(a(n).lt.amax)goto 100
         do i=1,n
            in(i)=i
         enddo
         close(88)
		 it=1
c*****************************************************************
c*****************************************************************
c     CALCULO DE LA ACRECION PLANETARIA
c*****************************************************************
c*****************************************************************
 200     continue               ! lazo por paso de tiempo
         do lt=1,nsal
            det=dexp(t/taugas)
            do k=1,n
               i=in(k)
c     ----------------------------------------------------------------
c     densidad volumetrica del gas en el plano medio del disco 
               ro_0=8.33d-13*(emestar/emesol)**0.5*sigmag_0/rc**1.25*
     &              (a(i)/rc)**(-gama-1.25)*exp(-(a(i)/rc)**(2-gama))*
     &              elestar**0.125 !esta es la INICIAL y esta en g/cm3
               rogas=ro_0/det   !esta es la que VA QUEDANDO por disipacion
c     calcula la probabilidad de eyeccion mas que de acrecion
               if(emegas(i).gt.0.d0)then
                  x=(efe*1.4d0*(emestar/emepla(i))**(2./3.))**2 !VER DE IL04
			   elseif(emepla(i).eq.0.d0)then
				   x=0.d0
               else
                  x=(rplanet(i)*emestar/(2.d0*emepla(i)*a(i)*uacm))**2
               endif
               if(x.gt.1.d0)x=1.d0
               if(x.lt.0.1d0)x=0.d0
               dmdt=faktor*sigmap(i)*rogas**(2./5.)*emepla(i)**(2./3.)
               dmdt=dmdt*yearsec/(a(i)*uacm)**(1./10.) ! masa acretada por ano
c     agrandamiento de la seccion eficaz si el planeta tiene al menos
c     media masa terrestre de gas ligado, pero con un limite maximo
c     de 100 veces (segun Chambers, Icarus 180 (2006) 496-513)
               rhc=1.d0
               if(emegas(i).gt.0.5d0.and.sigmap(i).ne.0.d0)then         
                  rcoc=acre_gas*a(i)**1.5*yearsec*(emepla(i)/emestar)**2
                  rcoc=(rcoc/(cappa*r_med*sigmap(i)))**0.5 !cociente de los radios**2
                  if(rcoc.gt.rhc)rhc=rcoc
                  if(rhc.gt.100.d0)rhc=100.d0 !si el cociente es mayor que 100 stop
               endif   
               dmdt=dmdt*rhc    !aumento en la tasa de acrecion por este efecto
               emedot_crit(i)=dmdt*1.d6/emet !coef para el calculo de la masa critica
               dm=dmdt*dt       !en gr
               emepla(i)=emepla(i)+dm*x
               erreh(i)=a(i)*(emepla(i)/(3.d0*emestar))**0.33333333 !nuevo r de hill
			   if(dm.ne.0..or.erreh(i).ne.0.)then
               sigmap(i)=sigmap(i)-dm/(2.d0*pi*b*a(i)*erreh(i)*uacm**2) !le resta a 
c     la densidad de solidos lo que ya comio el planeta
               rplanet(i)=rplanet(i)*(1.d0+dm*x/(3.*emepla(i))) !nuevo radio del planeta en cm, lo calculo haciendo el rf=ri+dr y sale. 
	   	       endif
c     ................................................................................
c     calculo el nuevo momento angular ganado debido a la acrecion de planetesimales en una masa equivalente a dm
               arriba=6.54d-7*ro_0**0.4*emestar**(1./6.)*
     &              emepla(i)**(2./3.)
               abajo=ropla**(11./15.)*(a(i)*1.5d13)**(0.1)
               eleang(i,3)=eleang(i,3)+(arriba/abajo)*dm 
c     ................................................................................
c     calcula la masa maxima de gas que podria llegar a acretar, 
c     calculando todo el gas que hay disponible a partir de este momento
c     (pues algo de gas ya se disipo).    
               if(indgas(i).ne.1)then!xq la emeiso se calcula una sola vez
                  rogas=rogas*taugas/(taugas-t)
                  rogas=rogas*(1.d0-det/enumber)
                  emeiso(i)=rogas*s*4*pi*b/
     @                 ((3.d0*emestar)**0.333333)
                  emeiso(i)=emeiso(i)*(a(i)**1.25d0)*uacm*(a(i)*uacm)**2
                  emeiso(i)=emeiso(i)**1.5
               endif
c     acrecion gaseosa en terminos de la contraccion K-H
c     standard segun Ida y Lin (usa esos mismos coeficientes)
               em=emepla(i)/emet!em es la masa en masas terrestres
               crit_gas=emecrit*emedot_crit(i)**0.25 !masa critica(IL)
               if(em.gt.crit_gas.and.emepla(i).lt.emeiso(i))then
                  indgas(i)=1         
c     -------------------------------------------------------------
c     tasa de acrecion de gas
c                  tkh=0.16449d10/em**1.911294 !coef viejos fit de andre
c     estos son los nuevos valores firteados de lo ultimo de Andre!
                  tkh=8.352913363d10/em**4.889865!fit 17/12/09
                  tasa=em/tkh
c     -------------------------------------------------------------
                  if(tasa.lt.1.d-6)then
                     em=em+tasa*dt
                     emegas(i)=emegas(i)+tasa*dt
                  else          !esto es porque no dejamos que acrete 
c     mas de una masa terrestre por ano
                     em=em+dt*em/1.d6
                     emegas(i)=emegas(i)+dt*em/1.d6
                  endif 
                  emepla(i)=em*emet!la masa vuelve a estar en gr
               endif
               erreh(i)=a(i)*(emepla(i)/(3.d0*emestar))**0.33333333
               if(sigmap(i).lt.1.d-20)sigmap(i)=0.d0
c     calculo el w, wcrit y prot porque ahora adquieren L por acrecion de planetesimales
			   if(emepla(i).ne.0..or.rplanet(i).ne.0.)then
               w(i)=dsqrt(eleang(i,1)**2+eleang(i,2)**2+eleang(i,3)**2)/
     &              (0.4d0*emepla(i)*rplanet(i)**2) 
               wcrit(i)=(g*emepla(i)/(rplanet(i)**3.))**0.5 
               prot(i)=6.28d0/(w(i)*3600.d0)
			   endif
            enddo
c     ********************************************************************
c     me fijo para cada planeta si no supera el wcrit, en caso de superarlo, 
c     el planeta desaparece. Primero me fijo para todos hasta n-1
            kk=1
 777        if(kk.ne.n)then
               i=in(kk)
               if(w(i).gt.wcrit(i))then
                  if(emepla(i)/5.97d27.le.10.0d0.and.a(i).gt.0.d0)then
                     if(prot(i).lt.1.0d15)then
                        elemodulo=sqrt(eleang(i,1)**2+eleang(i,2)**2+
     &                       eleang(i,3)**2)
                    endif
                  endif
                  do k=kk,n-1
                     j=in(k)
                     jm1=in(k+1)
                     emepla(j)=emepla(jm1)
                     a(j)=a(jm1)
                     erreh(j)=erreh(jm1)
                     sigmap(j)=sigmap(jm1)
                     emeiso(j)=emeiso(jm1)
                     emegas(j)=emegas(jm1)
                     rplanet(j)=rplanet(jm1)
                     eleang(j,1)=eleang(jm1,1)
                     eleang(j,2)=eleang(jm1,2)
                     eleang(j,3)=eleang(jm1,3)
                     prot(j)=prot(jm1)
                     obli(j)=obli(jm1)
                     realo(j)=realo(jm1)
                     wcrit(j)=wcrit(jm1)
                     w(j)=w(jm1)
                     ncol(j)=ncol(jm1)
                  enddo
                  n=n-1
                  goto 777
               else
                  kk=kk+1
                  goto 777
               endif
            else
               i=in(n)
               if(w(i).gt.wcrit(i))then
                  if(emepla(i)/5.97d27.le.10.0d0.and.a(i).gt.0.d0)then
                     if(prot(i).lt.1.0d15)then
                        elemodulo=sqrt(eleang(i,1)**2+eleang(i,2)**2+
     &                       eleang(i,3)**2)
                     endif
                  endif
                  n=n-1
               endif
            endif
c     *******************************************************************
c     analisis de posibles acreciones entre embriones     
c     cuando estan a menos de superpos*RH
            do l=1,n-1
               i=in(l)
               im1=in(l+1)
               da=superpos*(erreh(i)+erreh(im1))
               if(da.gt.a(im1)-a(i))then
c     ******************************************************************
c     ******************************************************************
c     Incluimos el calculo de la nueva oblicuidad despues de la colision
c     oblicuidades.f
                  a1=a(im1)*1.5d13
                  a2=a(i)*1.5d13
                  delta=a1-a2   !separacion entre las orbitas
                  eme2=emepla(i)-emegas(i)*emet
                  eme1=emepla(im1)-emegas(im1)*emet
                  rpla1=rplanet(im1)
                  rpla2=rplanet(i)
                  ele2(1)=eleang(i,1)
                  ele2(2)=eleang(i,2)
                  ele2(3)=eleang(i,3)
                  ele1(1)=eleang(im1,1)
                  ele1(2)=eleang(im1,2)
                  ele1(3)=eleang(im1,3)
                  call spin(dseed,eme1,eme2,rpla1,rpla2,a1,a2,delta,
     &                 ele1,ele2,ele,ar,cont,emestar)
                  eme=eme2+eme1
                  ncol(i)=ncol(i)+1
                  rplanet(i)=rpla2*(eme/eme2)**0.333333
                  emepla(i)=eme
                  a(i)=ar/1.5d13
                  erreh(i)=a(i)*(emepla(i)/(3.d0*emestar))**0.3333333
                  sigmap(i)=(sigmap(i)+sigmap(im1))/2.d0
c     una vez que ya se produjo el choque calculamos la oblicuidad con la
c     que queda el planeta resultante.
                  eleang(i,1)=ele(1)
                  eleang(i,2)=ele(2)
                  eleang(i,3)=ele(3)
c     escribe la oblicuidad en grados y el periodo de rotacion en horas
                  dinercia=0.4d0*emepla(i)*rplanet(i)**2 !momento de inercia de una esfera
                  w(i)=dsqrt(eleang(i,1)**2+eleang(i,2)**2+
     &                 eleang(i,3)**2)/dinercia !mod velocidad angular, 
                  wcrit(i)=(g*emepla(i)/(rplanet(i)**3.))**0.5 !en 1/seg 
c     la obtiene usando que L=Iw
                  prot(i)=6.28d0/(w(i)*3600.d0) !prot=2pi/(w.seg en una hora), me da el periodo
c     de rotacion
                  obli(i)=datan(eleang(i,3)/dsqrt(eleang(i,1)**2+
     &                 eleang(i,2)**2))*180.d0/3.1416 !en grados
                  realo(i)=90.0d0-obli(i) !oblicuidad en grados
c     ******************************************************************
c     *******************************************************************
c     cambia los indices de los planetas que estan mas afuera
c     pues, como resultado de la acrecion de dos, hay un planeta menos
                  do k=im1,n-1
                     j=in(k)
                     jm1=in(k+1)
                     emepla(j)=emepla(jm1)
                     a(j)=a(jm1)
                     erreh(j)=erreh(jm1)
                     sigmap(j)=sigmap(jm1)
                     emeiso(j)=emeiso(jm1)
                     emegas(j)=emegas(jm1)
                     rplanet(j)=rplanet(jm1)
                     eleang(j,1)=eleang(jm1,1)
                     eleang(j,2)=eleang(jm1,2)
                     eleang(j,3)=eleang(jm1,3)
                     prot(j)=prot(jm1)
                     obli(j)=obli(jm1)
                     realo(j)=realo(jm1)
                     wcrit(j)=wcrit(jm1)
                     w(j)=w(jm1)
                     ncol(j)=ncol(jm1)
                  enddo
                  n=n-1
c     ahora, si ademas el valor de w supera al wcrit, entonces el nuevo planeta  desaparece. Para hacerlo desaparecer cambiamos nuevamente los indices de los planetas que siguen.
                  if(w(i).gt.wcrit(i))then
                    if(emepla(i)/5.97d27.le.10.0d0.and.a(i).gt.0.d0)then
                        if(prot(i).lt.1.0d15)then
                           elemodulo=sqrt(eleang(i,1)**2+eleang(i,2)**2+
     &                          eleang(i,3)**2)
                        endif
                     endif
                     do k=i,n-1
                        j=in(k)
                        jm1=in(k+1)
                        emepla(j)=emepla(jm1)
                        a(j)=a(jm1)
                        erreh(j)=erreh(jm1)
                        sigmap(j)=sigmap(jm1)
                        emeiso(j)=emeiso(jm1)
                        emegas(j)=emegas(jm1)
                        rplanet(j)=rplanet(jm1)
                        eleang(j,1)=eleang(jm1,1)
                        eleang(j,2)=eleang(jm1,2)
                        eleang(j,3)=eleang(jm1,3)
                        prot(j)=prot(jm1)
                        obli(j)=obli(jm1)
                        realo(j)=realo(jm1)
                        wcrit(j)=wcrit(jm1)
                        w(j)=w(jm1)
                        ncol(j)=ncol(jm1)
                     enddo
                     n=n-1   
                  endif
               endif
            enddo
c     termino de chequear las posibles interacciones entre embriones
c     migracion de planetesimales debido al drag con el gas nebular
c            call migra_planetes(n,a,sigmap,in)
c     migracion de planetas Tipo II debido la interaccion Planeta-Disco
            rm=10.d0*dexp(0.4d0*t/taugas) !esta en UA
            call migrar(n,a,rm,alfa,in,gama,rc,sigmag_0)
c     migracion Tipo I de planetas debido a efecto de marea con el disco de gas
            call migrarI(n,a,in,alfa,gama,rc,elestar,sigmag_0)
c     como los planetas pueden migrar y pasar mas adentro o mas
c     afuera que otros, se  puede  armar  lio en la busqueda
c     de acrecion entre embriones. Para evitarlo se 
c     los indexa y se ordenan los indices in(i) 
c     por distancia a la estrella central
            if(n.eq.1)goto 999
            if(n.lt.1)goto 666
            do k=1,n
               bb(k)=in(k)
            enddo
            call indexx(n,a,in)
 999        continue
 		 lscond=0
		 cond=0
c	 print*, t
         do k=1,n
            i=in(k)
c            elemodulo=sqrt(eleang(i,1)**2+eleang(i,2)**2+
c     &           eleang(i,3)**2)
c    solo adopto como planetas objetos con masas menores a 10000 masas terrestres y que son mayores 
c    que la masa de Mercurio
            if(emepla(i)/emet.le.26000.0d0.and.emepla(i)/emet.gt.
     &           0.05d0)then
               if(a(i).gt.0.0d0)then
		 cond=1
c      	         if(t-tini.eq.1d5.or.t-tini.eq.1d6.or.t-tini.eq.1d7.
c     &           or.t-tini.eq.2d7)then
                if(t-tini.eq.2d7)then
                   write(444,*)ident,',',it,',',t,',',a(i),',',
     &             emegas(i),',',emepla(i)/emet,',',rplanet(i)/radtie,
     &             ',',emestar,',',rc,',',qest,',',sigmag_0,',',
     &             emed,',',gama,',',apert,',',fpert,',',constmigI,
     &             ',',emetal,',',taugas
		 endif
		 if(it.eq.lsac(lsk))then
                   lscond=1	
		   if(verb)then
                     write(333,*)it,t,a(i),emegas(i),emepla(i)/emet,
     &               rplanet(i)/radtie,emestar,rc,qest,sigmag_0,
     &               emed,gama,apert,fpert,constmigI,emetal,taugas
		   endif
		   if(emepla(i)/emet.gt.20.0d0)then
                     write(555,*)ident,t,a(i),emegas(i),emepla(i)/emet
		   endif
 		 endif	
               endif
            endif			
         enddo
		 if(lscond.eq.1)lsk=lsk+1
		 if(cond.eq.1)it=it+1
c     actualiza el tiempo   
            t=t+dt     
         enddo 
         if(t.lt.tfin)goto 200  
         call flush(20)
         print*,ident," run finished correctly - at tfin"
 666     continue
	close(333)
      enddo
 600  format(d15.4,1x,d12.4,1x,d15.4,1x,d12.4,1x,d12.4,1x,d12.4,
     &     1x,i4,1x,i4,1x,d15.4,1x,d12.4,1x,d12.4,
     &     1x,d15.4)
      stop
      end
c==========================================================================
c==========================================================================
c     subrutina que calcula el cambio de momento angular l por la colision 
c     todas las unidades seran c.g.s.
c     entra: dseed: numero semilla
c     m1,m2: masas planetarias
c     rpla1,rpla2: radios planetarios
c     a1,a2: semiejes de los planetas
c     delta: separacion entre los planetas
c     l1,l2: vector momento angular de spin de los planetas
c     sale   l: vector momento angular del planeta resultante
      subroutine spin(dseed,eme1,eme2,rpla1,rpla2,a1,a2,delta,
     &     ele1,ele2,ele,a,cont,emestar)
      implicit real*8(a-h,o-z)
      dimension ele1(3),ele2(3),ele(3)
      dimension v(3),pa(3),torque(3)
      dseed=dseed+cont
      cont=cont+1.0
      a=(a1+a2)/2.d0!la nueva posicion, en el medio de las anteriores,[cm]
c     calcula la velocidad de escape
      call vesc(eme1,eme2,rpla1,rpla2,ve) 
c     calcula la velocidad relativa de colision
      call velo(a,delta,ve,vrel,emestar)
c     se asume que el mas grande es el 2
      rpla=rpla2
      empla=eme1
c     el momento angular de spin se conserva
      do i=1,3
         ele(i)=ele1(i)+ele2(i)
      enddo
c     pero si el mas grande es el 1...
      if(rpla1.gt.rpla2)then
         rpla=rpla1
         empla=eme2
      endif
c     calcula las componentes de la vel. relativa y el punto de aplicacion
      call vr(dseed,rpla,vrel,v,pa) 
c     calcula el torque en el planeta mayor debido al impacto
      call torquesub(empla,v,pa,torque)
c     calcula la variacion de momento angular
      call dl(ele,torque)
      return
      end
c--------------------------------------------------------------
c     subrutina que calcula las tres componentes de la velocidad
c     y el punto de aplicacion del impulso de la colision
c     se supone que el menos masivo es el proyectil
c     (el que le pega al mas masivo)
      subroutine vr(dseed,rpla,vrel,v,pa)
c     entra:
c     dseed: numero semilla
c     rpla: radio del planeta impactado
c     vrel: velocidad relativa
c     sale:
c     v: componentes de la velocidad de impacto
c     pa: punto de aplicacion del impulso
      implicit real*8(a-h,o-z)
      real*8 lambda,mu,l
      dimension v(3),pa(3)
      pi=4.d0*datan(1.d0)
c     calcula beta y lambda
      lambda=2.d0*pi*rnd(dseed)!es un angulo al azar
      clambda=dcos(lambda)
      slambda=dsin(lambda)
c     calcula beta
      xxx=1.d0-2.d0*rnd(dseed)!es un numero al azar   
      beta=dasin(xxx)!es un angulo al azar
      cb=dcos(beta)
      sb=dsin(beta)
c     calcula el punto de aplicacion (en coordenadas esfericas con 
c     centro en el centro del planeta mas masivo, este punto es al azar)
      pa(1)=rpla*cb*clambda
      pa(2)=rpla*cb*slambda
      pa(3)=rpla*sb
c     calcula los angulos eta y mu segun la receta de Henon
      eta=2.d0*pi*rnd(dseed)
      se=dsin(eta)
      ce=dcos(eta)
      mu=dasin(dsqrt(rnd(dseed)))
      sm=dsin(mu)
      cm=dcos(mu)
c     calcula l y theta
      l=dacos(sb*cm-cb*sm*ce)
      sl=dsin(l)
      cl=dcos(l)
      st=sm*se/sl
      ct=(cm*sb+sm*cb*ce)/sl
      theta=datan2(st,ct)
c     calcula s
      s=lambda-theta
      v(1)=-vrel*sl*dcos(s)!componente x de la velocidad de impacto
      v(2)=-vrel*sl*dsin(s)!vy de impacto
      v(3)=-vrel*cl!vz de impacto
      return
      end
c-----------------------------------------------------
c     subrutina que calcula la velocidad relativa entre
c     dos planetas separados una distancia delta, ambos en orbitas
c     coplanares, suponiendo que delta es << que los semiejes de ambos 
c     planetas.
      subroutine velo(a,delta,ve,v,emestar)
c     a: semieje promedio de ambos semiejes
c     delta: separacion entre los dos planetas
c     ve: velocidad de escape mutua (se debe sumar para calcular
c     la velocidad de colision)
c     v: velocidad tipica de colision (aproximadamente)
      implicit real*8(a-h,o-z)
      gm=6.672d-8*emestar!G+Mstar
      omega=dsqrt(gm/a**3)
      v=omega*delta/2.d0
      v=v+ve
      return
      end
c-------------------------------------------------------
c     subrutina que dadas la velocidad, la masa del impactor
c     y el punto de aplicacion, calcula el torque
      subroutine torquesub(eme,v,pa,torque)
c     aca m es la masa del mas chico, ya que ese es el embrion que 
c     aplica la fuerza, el que genera el torque
      implicit real*8(a-h,o-z)
      dimension v(3),pa(3),torque(3)
      torque(1)=eme*(v(3)*pa(2)-v(2)*pa(3))
      torque(2)=eme*(v(1)*pa(3)-v(3)*pa(1))
      torque(3)=eme*(v(2)*pa(1)-v(1)*pa(2))
c     este es el torque por delta t. Consideramos que la vinicial 
c     es cero frente a la de impacto. Por eso usamos directamente v para el
c     calculo
      return
      end
c------------------------------------------------------
c     subrutina que cambia el vector momento angular de spin 
c     de un planeta, dado un cierto torque
      subroutine dl(ele,torque)
      implicit real*8(a-h,o-z)
      dimension ele(3),torque(3)
      do i=1,3
         ele(i)=ele(i)+torque(i)!recordemos que tau es torque por delta t
      enddo
      return
      end
c--------------------------------------------------------
c     subrutina que calcula la velocidad de escape de ambos planetas
      subroutine vesc(eme1,eme2,r1,r2,ve)
      implicit real*8(a-h,o-z)
      g=6.668d-8
      ve=dsqrt(2.d0*g*(eme1+eme2)/(r1+r2))
      return
      end
C---------------------------------------------------------------------------
C GENERACION DE UN NUMERO AL AZAR:
      DOUBLE PRECISION FUNCTION RND(DSEED)
      DOUBLE PRECISION DSEED,D2A32
      D2A32=2.D0**32
      DSEED=DMOD(3125.*DSEED,D2A32)
      RND=DSEED/D2A32
      RETURN
      END
c==================================================================
c==================================================================
c     subrutina que migra planetesimales debido al efecto del drag
c     (Ver Thommes et al. Icarus 161, 2003)
      subroutine migra_planetes(n,a,sigmap,in)
      implicit real*8(a-h,o-z)
      dimension a(100000),sigmap(100000),emepla(100000)
      dimension dadt(100000),erreh(100000),in(100000)
      common/pasar/ro_ice,r,g,dt,pi,emesol,uacm,yearsec,
     1     det,c1
      common/pla/emepla,b,erreh,s,amin,emestar
      do j=1,n-1
         i=in(j)
         eta1=1.6d-3*a(i)**1.5
         vk1=dsqrt(g*emestar/(a(i)*uacm))
         rogas1=(fac_g*1.4d-9*a(i)**(-2.75))/det
         tau1=2.6667d0*ro_ice*r/(rogas1*vk1)
         rh1=(emepla(i)/(3.d0*emestar))**0.33333
         eq1=2.7d0*rh1*(r*ro_ice/(b*a(i)*uacm*rogas1))**0.2
         corch1=dsqrt(eta1**2+0.75d0*eq1**2)*(eta1+eq1**2)
         dadt(i)=-2.d0*a(i)*uacm*corch1/tau1
      enddo
      do j=1,n-2
         i=in(j)  
         a1=a(i)*uacm
         a2=a(i+1)*uacm 
         da=a2-a1
         ds2=a2*sigmap(i+1)*dadt(i+1)
         ds1=a1*sigmap(i)*dadt(i)
         ds=-(ds2-ds1)/(a(i)*uacm)
         sigmap(i)=sigmap(i)+ds*dt/(a(i)*uacm)
         if(sigmap(i).lt.1.d-20)sigmap(i)=0.d0
      enddo
      return
      end
c     subrutina que migra planetas gigantes con migracion tipo II
c     (Ver Ida y Lin, The Astrophysical Journal, 604:388, 2004)
      subroutine migrar(n,a,rm,alfa,in,gama,rc,sigmag_0)
      implicit real*8(a-h,o-z)
      dimension emepla(100000),a(100000),erreh(100000)
      dimension in(100000)
      common/pasar/ro_ice,r,g,dt,pi,emesol,uacm,yearsec,
     .     det,c1
      common/pla/emepla,b,erreh,s,amin,emestar
      common/others/apert,fpert,constmigI
      sigmag=sigmag_0*(rm/rc)**(-gama)*exp(-(rm/rc)**(2-gama))*
     &     (1+apert*cos(2.*pi*rm/(fpert*0.03361386903*rm**(5./4.))))
      sigmag=sigmag/det 
      om=dsqrt(g*emestar/(rm*uacm)**3)
      ememe=emestar/emesol
      alfa3=alfa/1.d-3
      coef=ememe*alfa3*5.97d27*30.d0
      coef2=1.2d2*ememe*5.97d27
      do j=1,n
         i=in(j)
         if(a(i).gt.amin)then   ! no migra si llego al borde interno
            eme_gap=coef2*a(i)**0.75 !masa para que abra un gap IL04 en gr
            if(emepla(i).gt.eme_gap)then ! migra si abrio un gap (Tipo II)     
               sg=(a(i)-rm)/dabs(a(i)-rm)
               tau2=4.217d14*a(i)**0.5*(emepla(i)/5.97d27)*
     &              ememe**(-0.5)*rm**(-0.25)/(alfa*sigmag)!en anos
               dadt=sg*a(i)/tau2
               da=dadt*dt
               if(dabs(da).lt.0.1d0)a(i)=a(i)+dadt*dt !porque migra unicamente 
c     si pasa esto y no cuando el cambio en a es mayor??resultado numerico!
            endif
         endif
      enddo
      return
      end
      subroutine migrarI(n,a,in,alfa,gama,rc,elestar,sigmag_0)
      implicit real*8(a-h,o-z)
      dimension emepla(100000),a(100000),erreh(100000)
      dimension in(100000)
      common/pasar/ro_ice,r,g,dt,pi,emesol,uacm,yearsec,
     .     det,c1
      common/pla/emepla,b,erreh,s,amin,emestar
      common/others/apert,fpert,constmigI
      ememe=emestar/emesol
      alfa3=alfa/1.d-3 
      coef=ememe*alfa3*5.97d27*30.d0
      coef2=1.2d2*ememe*5.97d27
      do j=1,n
         i=in(j)
         if(a(i).gt.amin)then   ! no migra si llego al borde interno
            eme_gap=coef2*a(i)**0.75 !masa para que abra un gap IL04
            if(emepla(i).lt.eme_gap)then !no migra por migI si migra por migII
               sigmag=sigmag_0*(a(i)/rc)**(-gama)*
     &              exp(-(a(i)/rc)**(2.-gama))*
     &     (1+apert*cos(2.*pi*a(i)/(fpert*0.03361386903*a(i)**(5./4.))))
               sigmag=sigmag/det
               cteaux=2.7+1.1*
     &              (gama+(2-gama)*0.43429448*(a(i)/rc)**(2-gama))
               dadt=9.46d-10*cteaux*(emepla(i)/5.97d27)*ememe**0.5*
     &              sigmag*elestar**(-0.25)*a(i)
               dadt=-dadt*constmigI
               da=dadt*dt
               if(dabs(da).lt.0.1d0)a(i)=a(i)+dadt*dt
            endif
         endif
      enddo
      return
      end
C***************************************************************************
C ==> ESTA RUTINA GENERA UN VECTOR DE INDICES DE ACUERDO AL ORDEN
C ==> CRECIENTE DEL VECTOR ARRIN (NUMERICAL RECIPES)
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
      INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
      IF(L.GT.1)THEN
      L=L-1
      INDXT=INDX(L)
      Q=ARRIN(INDXT)
      ELSE
      INDXT=INDX(IR)
      Q=ARRIN(INDXT)
      INDX(IR)=INDX(1)
      IR=IR-1
      IF(IR.EQ.1)THEN
      INDX(1)=INDXT
      RETURN
      ENDIF
      ENDIF
      I=L
      J=L+L
20    IF(J.LE.IR)THEN
      IF(J.LT.IR)THEN
      IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
      ENDIF
      IF(Q.LT.ARRIN(INDX(J)))THEN
      INDX(I)=INDX(J)
      I=J
      J=J+J
      ELSE
      J=IR+1
      ENDIF
      GO TO 20
      ENDIF
      INDX(I)=INDXT
      GO TO 10
      END
