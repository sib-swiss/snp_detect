        Parameter        (NERR=0)

        Integer*02        ISEQ(1000000)
        Character*01      CSNP(1000000)
        Integer*02        IMAT(10000000)
        Character         CMAT(10000000)

        Integer           IALL(1000)
        Character         CALL(1000)
        Integer           KALL(1000)       

        Integer           IBEG(1000)       
        Integer           IEND(1000)       

* Sequence header

        Character*14      FSEQ(10000)
        Character*10      FCHR(10000)
        Character*01      CSTR(10000)
        Character*01      CTYP(10000)
        Integer           IPOS(10000)

        Logical           LREP

        Character*1024    RCIN

* SNP search parameters

        KTUP=7
        KDEL=10
        KNIB=2

        NMIN=2
        RMIN=0.1
	NCW=2
	NRW=2

        KLEN=20

        Call Repar
     *    (LREP,KTUP,KDEL,KNIB,NMIN,RMIN,NCW,NRW,KLEN,IRC)
        If(IRC.NE.0) then
           Write(NERR,'(/,
     *     ''Usage: snp_detect [options]  '',//,
     *     ''  where options are:'',/,
     *     ''    -r        generate report format '',/,
     *     ''    -k <int>  k-tuple length [7]'',/,
     *     ''    -d <int>  maximal SNP gap length [10]'',/,
     *     ''    -b <int>  nibbling range [2]'',/,
     *     ''    -n <int>  min. # of sequences per allele [2]'',/,
     *     ''    -c <int>  weight multiplier for chr sequences [2]'',/,
     *     ''    -u <int>  weight multiplier for rna sequences [2]'',/,
     *     ''    -f <int>  min. fraction of sequences per allele'', 
     *                  '' [0.1]'',/,
     *     ''    -s <int>  min. spacer length between SNPs [20]'',//, 
     *     ''    (reads from stdin, writes to stdout)'',// 
     *       )')
        Stop
        End if

* sequence input 

        NLEN=0
        NSEQ=0
        Read(5,'(A)',End= 50) RCIN

    1   If(RCIN(1:1).NE.'>') then 
           Read(5,'(A)',End= 50) RCIN
           Go to   1
        End if

        Call ParseHeader(   1,FSEQ,FCHR,CTYP,CSTR,IPOS,RCIN) 

        J1=0
    2   Read(5,'(A)',Iostat=IOS) RCIN
        If(RCIN(1:1).EQ.'>'.OR.IOS.NE.0) then
C          Print *,J1
           If(NLEN.EQ.0) then
              NLEN=J1
           Else if(J1.NE.NLEN) then
              Write(6,'(''Error unequal sequence length detected.'')')
              Stop
           End if
           If(IOS.EQ.0) 
     *        Call ParseHeader(NSEQ+2,FSEQ,FCHR,CTYP,CSTR,IPOS,RCIN) 

* k-tuple encoding / write one row of IMAT
* 
* rules for k-tuple encoding (e.g. k=3): 
*
* 1.  #(AAA)=0, #(CAA)=1, #(TTT)=255
* 2.  #(AAA-)=0, 
* 3.  #(AAN)=0, #(ANA)=0, #(NAA)=0 


           K1=(NSEQ*NLEN)+1
           I=0 
           K=4**(KTUP-2) 
           M=0
           N=0
           IMAT(K1)=-1
              I1=1
           If(ISEQ(1).GE.0) then
              N=ISEQ(1)
              I=1
           End if
           Do I1=2,NLEN
              K1=K1+1
              If     (ISEQ(I1).EQ.-2) then 
                 M=0
                 IMAT(K1)=-1
                 I=0
                 N=0
              Else if(ISEQ(I1).EQ.-1) then
                 M=M+1
                 If(M.EQ.KDEL) then 
                    Do I2=K1-KDEL+1,K1
                       IMAT(I2)=-1
                    End do
                    I=0
                    N=0
                 Else
                    IMAT(K1)=IMAT(K1-1)
                 End if
              Else
                 M=0
                 If     (I.EQ.0) then
                    N=ISEQ(I1) 
                    I=1
                    IMAT(K1)=-1
                 Else if(I.GT.K) then
                    N=N/4+I*ISEQ(I1)
                    IMAT(K1)=N
                 Else if(I.EQ.K) then
                    I=I*4
                    N=N+I*ISEQ(I1) 
                    IMAT(K1)=N
                 Else      
                    I=I*4
                    N=N+I*ISEQ(I1) 
                    IMAT(K1)=-1
                 End if
              End if
           End do

* border nibbling 

           K1=NSEQ*NLEN+1 
           K2=(NSEQ+1)*NLEN
           M=0
   21      Continue 
           If     (IMAT(K1).LE.-1) then
              M=M+1
              K1=K1+1
              If(K1.GT.K2) go to  22
           Else if(M.GE.KDEL) then
              M=0
              Do I2=1,KNIB
                 IMAT(K1)=-1
                 K1=K1+1
                 If(K1.GT.K2) go to  22
              End do
           Else
              M=0
              K1=K1+1
              If(K1.GT.K2) go to  22
           End if
           Go to  21
   22      Continue
                         
           K1=(NSEQ+1)*NLEN
           K2=NSEQ*NLEN+1 
           M=0
   23      Continue 
           If     (IMAT(K1).LE.-1) then
              M=M+1
              K1=K1-1
              If(K1.LT.K2) go to  24
           Else if(M.GE.KDEL) then
              M=0
              Do I2=1,KNIB
                 IMAT(K1)=-1
                 K1=K1-1
                 If(K1.LT.K2) go to  24
              End do
           Else
              M=0
              K1=K1-1
              If(K1.LT.K2) go to  24
           End if
           Go to  23
   24      Continue

* done

C          K1=NSEQ*NLEN
C          Do I1=1,NLEN
C             K1=K1+1
C             Print *,I1,K1,ISEQ(I1),IMAT(K1)
C          End do

           NSEQ=NSEQ+1
           If(IOS.NE.0) Go to  50
           J1=0

        Else

* continue to read sequence

              K1=J1
              If(NSEQ.GE.1) K1=K1+(NSEQ*NLEN)
           Do I1=1,256
              If(RCIN(I1:I1).NE.' ') then
              J1=J1+1
              K1=K1+1
                 If     (RCIN(I1:I1).EQ.'A'.OR.RCIN(I1:I1).EQ.'a') then
                    ISEQ(J1)=0
                    CMAT(K1)='A'
                 Else if(RCIN(I1:I1).EQ.'C'.OR.RCIN(I1:I1).EQ.'c') then
                    ISEQ(J1)=1
                    CMAT(K1)='C'
                 Else if(RCIN(I1:I1).EQ.'G'.OR.RCIN(I1:I1).EQ.'g') then
                    ISEQ(J1)=2
                    CMAT(K1)='G'
                 Else if(RCIN(I1:I1).EQ.'T'.OR.RCIN(I1:I1).EQ.'t') then
                    ISEQ(J1)=3
                    CMAT(K1)='T'
                 Else if(RCIN(I1:I1).EQ.'-'.OR.RCIN(I1:I1).EQ.'.') then
                    ISEQ(J1)=-1
                    CMAT(K1)='-'
                 Else
                    ISEQ(J1)=-2
                    CMAT(K1)='N'
                 End if
              End if
           End do      

        End if
        Go to   2 

  50    Continue 

* Edit sequences

           J2=NLEN
        Do I1=2,NSEQ
        Do I2=1,NLEN
           J2=J2+1
           If(CMAT(J2).NE.'-'.AND.CMAT(J2).EQ.CMAT(I2)) CMAT(J2)='.'
        End do
        End do

* SNP test 

        Do I1=1,NLEN

           CSNP(I1)='-'

* - generate ktup allele repertoire

           NALL=0
           KSEQ=0
	      J2=0
           Do I2=I1,NSEQ*NLEN,NLEN
	      J2=J2+1
              If(IMAT(I2).GE.0) then
C                If(I1.EQ.10) Print *,CTYP(J2)
		 If     (CTYP(J2).EQ.'C') then 
		    K3=NCW
		 Else if(CTYP(J2).EQ.'R') then 
		    K3=NRW
                 Else
		    K3=1
                 End if
		 KSEQ=KSEQ+K3
                 Do I3=1,NALL
                    If(IMAT(I2).EQ.IALL(I3)) then
                       KALL(I3)=KALL(I3)+K3
                       Go to  60
                    End if
                 End do
                 NALL=NALL+1
                 IALL(NALL)=IMAT(I2)
                 KALL(NALL)=K3
   60            Continue
              End if
           End do

* - determine frequence cut-off

           JCUT=RMIN*(KSEQ)
              
* - check allele frequencies

           JSNP=0
C          Print *,I1,KSEQ,NMIN,JCUT,NALL,(KALL(ii1),ii1=1,NALL)
           Do I2=1,NALL
              If(KALL(I2).GE.NMIN.AND.KALL(I2).GT.JCUT) JSNP=JSNP+1
           End do
           If(JSNP.GE.2) CSNP(I1)='N'
        End do

* expand snps 

           K=0
	Do I1=1,NLEN
	   If(CSNP(I1).EQ.'-') then
	      K=K+1
           Else 
	      Do I2=I1-MIN(K,KTUP-1),I1-1
		 CSNP(I2)='N'
              End do 
           End if
        End do

* generate single base allele repertoire 

        Do I1=1,NLEN
	If(CSNP(I1).EQ.'N') then

           NALL=0
           KSEQ=0
	      J2=0
           Do I2=I1,NSEQ*NLEN,NLEN
	      J2=J2+1
              If(IMAT(I2).GE.0) then
		 If     (CTYP(J2).EQ.'C') then 
		    K3=NCW
		 Else if(CTYP(J2).EQ.'R') then 
		    K3=NRW
                 Else
		    K3=1
                 End if
		 KSEQ=KSEQ+K3
                 Do I3=1,NALL
                    If(CMAT(I2).EQ.CALL(I3)) then
                       KALL(I3)=KALL(I3)+K3
                       Go to  70
                    Else if(CMAT(I2).EQ.'.') then
                       KALL( 1)=KALL( 1)+K3
                       Go to  70
                    End if
                 End do
                 NALL=NALL+1
                 CALL(NALL)=CMAT(I2)
                 KALL(NALL)=K3
   70            Continue
              End if
           End do
              
* - check allele frequencies

           JSNP=0
           JCUT=RMIN*KSEQ
C          Print *,I1,KSEQ,NMIN,JCUT,NALL,(KALL(ii1),ii1=1,NALL)
           Do I2=1,NALL
              If(KALL(I2).GE.NMIN.AND.KALL(I2).GT.JCUT) JSNP=JSNP+1
           End do
           If(JSNP.LT.2) CSNP(I1)='-'

        End if
        End do
* clean snps 
           
C          K=1
C       Do I1=NLEN,1,-1
C          If(CSNP(I1).EQ.'N') then 
C             If(K.LT.KTUP) CSNP(I1)='-'
C             If(CMAT(I1).NE.'-') K=K+1
C          Else
C             If(K.NE.1) CSNP(I1+1)='N'
C             K=1    
C          End if
C       End do
           
* count snps

        ISNP=0
        NL=KLEN
        Do I1=1,NLEN
           If(CSNP(I1).EQ.'N') then
              If(NL.GE.KLEN) then
                 ISNP=ISNP+1
                 IBEG(ISNP)=I1   
                 IEND(ISNP)=I1   
                 NL=0
              End if 
              IEND(ISNP)=I1   
           Else 
              NL=NL+1
           End if
        End do
              
* print SNP report

        If(LREP) then

        Do I1=1,ISNP

* - SNP regions

        Write(6,'(''SNP region'',I4,'': from,'',I6,'' to '',I6,''.'')')
     *     I1,IBEG(I1),IEND(I1)

           If(IEND(I1)-IBEG(I1).GT.60) then
              JBEG=IBEG(I1)          
              JEND=IEND(I1)          
           Else
              JBEG=(IBEG(I1)+IEND(I1))/2-30
              JBEG=MAX(1,JBEG+1)
              JEND=MIN(JBEG+59,NLEN)
           End if      

        Write(6,'('''')')
        Write(6,'(''                                      '',132A)')
     *    (CSNP(ii1),ii1=JBEG,JEND)

* - compute starting positions of sequences

           J2=0
        Do I2=1,NSEQ

              K=0
              J3=J2
           Do I3=1,JBEG-1
              J3=J3+1
              If(CMAT(J3).NE.'-') K=K+1
           End do
           If(CSTR(I2).EQ.'+') then
              NPOS=IPOS(I2)+K
           Else
              NPOS=IPOS(I2)-K
           End if

* sequence relevant ?

              K=0
           Do I3=J2+JBEG,J2+IEND(I1)
              If(CMAT(I3).NE.'-') K=K+1
           End do

              L=0
           Do I3=J2+IBEG(I1),J2+JEND
              If(CMAT(I3).NE.'-') L=L+1
           End do

           If(K.GT.0.AND.L.GT.0) then
           Write(6,'(A1,'' '',A14,'' '',A10,'' '',A1,I8,'' '',132A)')
     *        CTYP(I2),FSEQ(I2),FCHR(I2),CSTR(I2),NPOS,
     *       (CMAT(ii1),ii1=J2+JBEG,J2+JEND)
           End if

           J2=J2+NLEN  
        End do
        Write(6,'('''')')
        End do

* print SNP sequence

        Else
        Write(6,'(''>SNP_MAP SnpCount: '',I4,'' ..'')') ISNP
        Write(6,'((60A))')(CSNP(ii1),ii1=1,NLEN)
        End if

        Stop
        End
*----------------------------------------------------------------------*
        Subroutine  Repar
     *    (LREP,KTUP,KDEL,KNIB,NMIN,RMIN,NCW,NRW,KLEN,IRC)

        Character*62      CARG

        Logical           LREP

        IRC=0
        LREP=.FALSE.

        N1=Iargc()
        I1=1

    1   If(I1.GT.N1) go to  100
        Call GetArg(I1,CARG)
        If     (CARG(1:1).NE.'-') then
           Go to 900 
        Else if(CARG(1:2).EQ.'-r') then
           LREP=.TRUE.
        Else if(CARG(1:2).EQ.'-k') then
           I1=I1+1
           If(I1.GT.N1) go to 900
           Call GetArg(I1,CARG)
           Read(CARG,*,Err=900) KTUP
        Else if(CARG(1:2).EQ.'-d') then
           I1=I1+1
           If(I1.GT.N1) go to 900
           Call GetArg(I1,CARG)
           Read(CARG,*,Err=900) KDEL
        Else if(CARG(1:2).EQ.'-b') then
           I1=I1+1
           If(I1.GT.N1) go to 900
           Call GetArg(I1,CARG)
           Read(CARG,*,Err=900) KNIB
        Else if(CARG(1:2).EQ.'-n') then
           I1=I1+1
           If(I1.GT.N1) go to 900
           Call GetArg(I1,CARG)
           Read(CARG,*,Err=900) NMIN
        Else if(CARG(1:2).EQ.'-f') then
           I1=I1+1
           If(I1.GT.N1) go to 900
           Call GetArg(I1,CARG)
           Read(CARG,*,Err=900) RMIN
        Else if(CARG(1:2).EQ.'-c') then
           I1=I1+1
           If(I1.GT.N1) go to 900
           Call GetArg(I1,CARG)
           Read(CARG,*,Err=900) NCW
        Else if(CARG(1:2).EQ.'-u') then
           I1=I1+1
           If(I1.GT.N1) go to 900
           Call GetArg(I1,CARG)
           Read(CARG,*,Err=900) NRW
        Else if(CARG(1:2).EQ.'-s') then
           I1=I1+1
           If(I1.GT.N1) go to 900
           Call GetArg(I1,CARG)
           Read(CARG,*,Err=900) KLEN
        Else 
           Go to 900         
        End if
        I1=I1+1
        Go to   1

  100  Return

  900  IRC=1
       Go to 100
       End
*----------------------------------------------------------------------*
        Subroutine ParseHeader(NSEQ,FSEQ,FCHR,CTYP,CSTR,IPOS,RCIN) 

        Character*14      FSEQ(10000)
        Character*10      FCHR(10000)
        Character*01      CSTR(10000)
        Character*01      CTYP(10000)
        Integer           IPOS(10000)
        Character*1024    RCIN
        Integer           Rindex
         
* sequence type

        CTYP(NSEQ)='U'
        If     (RCIN(2:4).EQ.'chr') then 
           CTYP(NSEQ)='C'
        Else if(RCIN(2:4).EQ.'rna') then 
           CTYP(NSEQ)='R'
        Else if(RCIN(2:4).EQ.'est') then 
           CTYP(NSEQ)='E'
        Else if(RCIN(2:4).EQ.'ore') then 
           CTYP(NSEQ)='O'
        End if

* sequence ID

        RCIN(5:5)=' '
        J=Index(RCIN,'|')-1
        FSEQ(NSEQ)=RCIN(5:J)

* sequence chromatogram ID

        FCHR(NSEQ)=' '
           J=Index(RCIN,']')
           If(Index('abcdefghijklmnopqrstuvwxyz',RCIN(J+ 1:J+ 1)).NE.0 
     *   .AND.Index('abcdefghijklmnopqrstuvwxyz',RCIN(J+ 2:J+ 2)).NE.0 
     *   .AND.Index('0123456789'                ,RCIN(J+ 3:J+ 3)).NE.0 
     *   .AND.Index('0123456789'                ,RCIN(J+ 4:J+ 4)).NE.0 
     *   .AND.Index('abcdefghijklmnopqrstuvwxyz',RCIN(J+ 5:J+ 5)).NE.0 
     *   .AND.Index('0123456789'                ,RCIN(J+ 6:J+ 6)).NE.0 
     *   .AND.Index('0123456789'                ,RCIN(J+ 7:J+ 7)).NE.0 
     *   .AND.Index('.'                         ,RCIN(J+ 8:J+ 8)).NE.0
     *   .AND.Index('rsx'                       ,RCIN(J+ 9:J+ 9)).NE.0
     *   .AND.Index('0123456789'                ,RCIN(J+10:J+10)).NE.0 
     *       ) FCHR(NSEQ)=RCIN(J+1:J+10) 

* strand

        If(Index(RCIN,' minus strand').NE.0) then
           CSTR(NSEQ)='-'
        Else
           CSTR(NSEQ)='+'
        End if

* positions

        J=Rindex(RCIN,',')
        RCIN(J:J)=' '
        J=Rindex(RCIN,',')
        RCIN(J:J)=' '
        J=Rindex(RCIN,',')
        RCIN(J:J)=' '
        K=Rindex(RCIN(1:J-1),' ')

        Read(RCIN(K:),*) N1,N2,N3

        If(CSTR(NSEQ).EQ.'+') then
           IPOS(NSEQ)=N1
        Else
           IPOS(NSEQ)=N3-N1+1
        End if

C       Print *,FSEQ(NSEQ),' ',FCHR(NSEQ),' ',CSTR(NSEQ),' ',IPOS(NSEQ)

        Return
        End     
*----------------------------------------------------------------------*
        Function          Lblnk(STRING) 
        Character*(*)     STRING

        L=Len(STRING)
        Lblnk=0

        Do I1=L,1,-1
           If(STRING(I1:I1).NE.' ') then
              Lblnk=I1
              Exit
           End if
        End do 

        Return
        End
*----------------------------------------------------------------------*
        Integer Function  RINDEX(STR1,STR2) 

        Character*(*)     STR1
        Character*(*)     STR2

        L1=Len(STR1)
        L2=Len(STR2)
    
        Rindex=0
        If(L1.LT.L2) go to 100

        Do   9 I1=L1-L2+1,1,-1
              J2=1
           Do   8 I2=I1,I1+L2-1
              If(STR1(I2:I2).NE.STR2(J2:J2)) go to   9
    8      Continue
           Rindex=I1
           Go to 100 
    9   Continue

  100   Return

        End
*----------------------------------------------------------------------*
