#legge i file di mrbayes e capisce che modello usa. Non puo capire quante categorei ci sono in un eventuale gamma 
import sys
from random import randint
import re
import os
from math import log

def Table(x,cases=[]):
    table={}
    for case in cases:
        table[case]=0
    if len(cases)!=0:
        table['NA']=0
    for obs in x:
        if table.has_key(obs):
            table[obs]+=1
        else:
            if not cases:
                table[obs]=1
            else:
                table['NA']+=1
    return table

def findThinning(nexusfilename, burnin, Exsampling):
    """finding the """
    nexus=open(nexusfilename, 'r')
    nexusL=nexus.readlines()
    line=[l for l in nexusL if l.strip()[:5]=='mcmcp']
    if len(line)==0:
        line=[l for l in nexusL if l.strip()[:4]=='mcmc']
    line=line[0]
    line=[l.split('=') for l in line.split()[1:]]
    D={}
    D.update(line[:-1])
    ngen=int(D['ngen'])
    samplefreq=int(D['samplefreq']) 
    print ngen
    ngen=ngen/samplefreq
    print ngen
    ngen=ngen-int(burnin)
    thinning=ngen/int(Exsampling)
    print samplefreq,ngen, Exsampling,
    return  thinning
class readmcmc:
    def __init__(self, infilename, lseqlist='', nexusfilename='', pathmb='./mb', headtreefile=False):
        #generali
        #localpath="/".join(infilename.split("/")[:-1])
        #if localpath=="": localpath="./"
        ##I added one so that slash was not included in baseinfilename
        #baseinfilename=infilename[infilename.find(localpath)+1+len(localpath):]
        #localdir=[x for x in os.listdir(localpath) if x.find(baseinfilename)> -1]
        self.path=pathmb
        self.nexus=nexusfilename
        self.lseq=lseqlist
        if lseqlist=='':
            self.partitionsize()
        #prepara file p
        self.pfile=open(infilename+'.run1.p','r')
        while 1:
            temp=self.pfile.readline()
            #print temp
            if temp.split()[0]=='Gen':
                break
        self.testata=temp
            # conta numero di partizioni perche se incontro un paremetro con la parola "all" devo sapere quante partizioni ho subito
        flag=False
        num=''
        part=[]
        for i in temp:
            if i=='}':
                if num!='all':
                    part+=map(int,num.split(','))
                flag=False
                num=''
            if flag:
                num=num+i
            if i=='{':
                flag=True
        if len(part)==0:
            part.append(1)
        part.sort()
        self.partition=int(part[-1])
            #definisci il modello dell'analisi
        self.struttura={'tree':[None], 'molt':[None],'Rvalori':'listone','gamma':[None],'inv':[None],'basi':'listone'}
        sinonimi={'TL':'tree','kappa':'Rvalori','r':'Rvalori','m':'molt','alpha':'gamma','pinvar':'inv','pi':'basi'}
        self.modellopos={}
        count=0
        for i in self.testata.split()[2:]:
            if len(part)>1:
                I=i.split('{')
                parlist=I[1][0:-1].split(',')
                if parlist[0]=='all':
                    parlist=map(str,range(1,self.partition+1))
            else:
                I=[i]
                parlist=['1'] 
            for par in parlist:
                if not self.modellopos.has_key(par):
                    self.modellopos[par]={}
                    for pezzo in self.struttura.keys():
                        self.modellopos[par][pezzo]=self.struttura[pezzo] 
                avatar=I[0].split('(')[0]
                radice=sinonimi[avatar]
                if (avatar!='kappa') & (self.struttura[radice]!=[None]):
                    speclist=I[0].split('(')
                    if len(speclist)>1:
                        spec=speclist[1][:-1]
                    if self.modellopos[par][radice]=='listone':
                        self.modellopos[par][radice]={}
                    self.modellopos[par][radice][spec]=count
                else: self.modellopos[par][radice]=[count]
            count+=1 
        #prepara file t   
        ntree=len(re.findall('TL',self.testata))
        self.tlistfile=[]
        for part in range(1,ntree+1):
            if ntree==1: inter=''
            else: inter='.tree'+str(part)
            self.tlistfile.append(open(infilename+inter+'.run1.t','r'))
        for treefile in self.tlistfile:
            flag=False
            count=0
            self.head=''
            while 1:
                temp=treefile.readline()
                if headtreefile:
                    self.head+=temp
                if flag: count+=1
                if re.match(" +translate",temp)!=None:
                    flag=True
                if flag and (re.search(';',temp)!=None):
                    break
            self.taxa=str(count)
    def skipline(self):
        for treefile in self.tlistfile:
            T=treefile.readline().split()
            if T==[]: return 'EOF'
            elif T[0]!='tree':
                return 'EOF'
        self.pfile.readline()
    def readline(self):
        self.tree=[]
        for treefile in self.tlistfile:
            T=treefile.readline().split()
            if T==[]:
                return 'EOF'
            elif T[0]!='tree':
                return 'EOF'
            self.tree.append(T[-1])
        self.param=self.pfile.readline().split()[2:]
        return 'continua...'
    def partitionsize(self):
        #interroga mb sul numero e dimensioni delle partizioni 
        #creo un nuovo file nex che non ha sicuramente un commando di esecuzione mcmc o di quit
        filenex=open(self.nexus,'r')
        nexcom=filenex.readlines()
        newnexcom=[x for x in nexcom if x[0:5] not in ['mcmc;', 'quit;']]
        filenewnex=open('nexcom.nex','w')
        filenewnex.writelines(newnexcom)
        filenewnex.close() 
        #creo un command file per leggere stat di partizioni
        comfile=open('CC','w')
        comfile.write('log start filename=pippo.log replace\n')
        comfile.write('execute nexcom.nex;\n')
        comfile.write('charstat;\n')
        comfile.write('log stop;\n')
        comfile.write('quit;\n')
        comfile.close()
        os.popen(self.path+' < CC')
        #leggi la risposta
        partfile=open('pippo.log','r')
        part=partfile.readlines()
        partfile.close()
        #cancello il pipe modificat 
        #os.popen('rm nexcom.nex')
        #estraine il significato
        #inizializazione
        partname={}
        flag=False
        partitions={}
        partitioncount=[]
        #leggi le righe
        #la struttura del commando charstat di mb non deve cambiare troppo
        #il nome delle diverse partizioni deve essere in un rigo che comincia per "Partition"
        #la definizione della partizione in uso deve essere nella riga "current"
        # nessuna riga dopo la tabella deve avere piu di 6 elementi (5 numerazione python)
        for line in part:
            L=line.split()
            if len(L)==0: continue
            if flag and len(L)>5:
                if L[4]=='Included':
                    if partitions.has_key(L[goodpos]):
                        partitions[L[goodpos]][0]+=1
                    else:
                        partitioncount.append([1])
                        partitions[L[goodpos]]=partitioncount[-1]
            if L[0]=='Partition':
                partname[ L[2][1:-1] ] =L[1]
            if L[0]=='Current':
                curpar=partname[L[-1]]
            if L[0]=='#':
                flag=True
                count=0
                for l in L:
                    #print l
                    if l==curpar:
                       # aggiungo 2 perche la testata della tabella ha due colonne in meno che la tabella medesima
                       goodpos=count+2
                       #print goodpos,'goodpos', curpar
                       break
                    count+=1
        self.partition=len(partitioncount)
        LL=[]
        for i in partitioncount:
            LL=LL+i
        self.lseq=LL
                

class simulatorMol:
    def __init__(self, reader,pathevolver='./evolver'):
        self.reader=reader
        self.basiStat={}
        self.statdist=[] 
        self.M=[]
        self.motivi={}
        self.counter=0
        #dati sul programma di simulazione
        self.path=pathevolver
        self.formato='2' #formato puap di output per evolver 0,1:seqs or patterns in paml format (mc.paml); 2:paup format (mc.nex)
        self.outfile='mc.nex'
        self.seme=str(randint(1,1000000))
        self.modelli=[None,'4']+4*[None]+ ['7']
        self.Rforma=['0']*self.reader.partition
        self.ordini={'Rvalori':{'0':['MettiUno'],'4':[0],'7':['C<->T','A<->T','G<->T', 'A<->C','C<->G']},'basi':['T','C','A','G']}
        comfile=open('RR','w')
        comfile.write('5\n0\n')
    def skipline(self):
        self.reader.skipline()
    def writeout(self):
        #un rigo in avanti
        segnale=self.reader.readline()
        if segnale=='EOF': return 'EOF'
        self.preparaArchivioMotivi()
        #ci sono invarianti?
        pinv=0
        #lista dei motivi con doppio accesso 

        #Scrivi un documento di commando per partizione
        #NB sull'albero, la posizione presa da modellopos sul file param coincide con la numerazione degli alberi nel self.tree.
        # Questo e' garantitto solo se TL e' sempre il primo o i primi parametri nel file param dopo Gen e Lnl
        for part in map(str,range(1,self.reader.partition+1)):
            if self.reader.modellopos[part]['Rvalori']!='listone':
                self.Rforma[int(part)-1] = self.modelli[len(self.reader.modellopos[part]['Rvalori'])]
            outfile=open('MCbase.dat','w')
            if self.reader.modellopos[part]['inv'][0]!=None:
                pinv=float(self.reader.param[self.reader.modellopos[part]['inv'][0]])
            else: pinv=0
            #print self.reader.lseq
            outfile.write("\n".join([self.formato,self.seme,'']))
            lse=int(round((1-pinv)*int(self.reader.lseq[int(part)-1])))
            pse=int(self.reader.lseq[int(part)-1])-lse
            outfile.write(' '.join([self.reader.taxa,str(lse),"1"])+'\n')
            if self.reader.modellopos[part]['molt']==[None]:
                outfile.write("-1 \n")
            else:
                outfile.write(self.reader.param[ self.reader.modellopos[part]['molt'][0] ]+"\n")
            outfile.write(self.reader.tree[self.reader.modellopos[part]['tree'][0]]+'\n')
            outfile.write(self.Rforma[int(part)-1] +'\n')
            for param in self.ordini['Rvalori'][self.Rforma[int(part)-1]]:
                #print param, part, self.Rforma[int(part)-1], self.ordini['Rvalori'][self.Rforma[int(part)-1]]
                if param=='MettiUno': outfile.write('1 ')
                else: outfile.write(self.reader.param[self.reader.modellopos[part]['Rvalori'][param]]+' ')
            outfile.write('\n')
            
            if self.reader.modellopos[part]['gamma']!=[None]:
                a=self.reader.param[self.reader.modellopos[part]['gamma'][0]]
            else: a=None
            if a==None:
                outfile.write('0 4')
            else:
                outfile.write(a+' 4')
            outfile.write('\n')
            Basi={}
            for param in self.ordini['basi']:
                basi=self.reader.param[self.reader.modellopos[part]['basi'][param]]
                Basi[param]=float(basi)
                outfile.write(basi+' ')
            outfile.write("\n  T   C   A   G\n")
            outfile.close()
            #finito di scrivere il documento di commando lo eseguo
            self.evolve()
            #calcolo la distribuzione dei motivi presenti tenendo conto di eventuali siti invaribili
            #print self.outfile, pse, Basi
            self.multistatistica(self.outfile,pse,Basi)
            #print self.M
        #dalla distribuzione dei motivi ricavo un valore scalare che archivio
        #self.puliscimotivi()
        self.statdist.append(self.monostatistica())
        #preparo lettore motivi per accogliere i nuovi dati
        #print "aggiungo un posto"
        #print self.M[0]
        
        #print self.M[0], sum([x[-1] for x in self.M])
        return 'continua,...'
            #default={'molt':-1,'Rforma':None,'Rvalori':None,'gamma':0,'inv':0,'basi':None}
    def evolve(self):
        #esegue il programma di simulazione
        os.popen(self.path+' < RR')
    def leggialineamentopaml(self, infilename):
        from Bio.Alphabet import IUPAC
        import Bio.Align.Generic
        infile=open(infilename,'r')
        alpha=IUPAC.IUPACAmbiguousDNA()
        AL=Bio.Align.Generic.Alignment(alpha)
        while 1:
            line=infile.readline()
            if line=='': break
            if line[0]!=' ':
                L=re.split(' ',line,1)
                if len(L)==2:
                    #print L
                    L[1]=re.sub(' ','',L[1])[:-1]
                    AL.add_sequence(L[0],L[1])
        return AL
    def leggialineamentoNexus(self, infilename):
        from Bio.Nexus.Nexus import Nexus
        from Bio.Alphabet import IUPAC
        import Bio.Align.Generic
        alpha=IUPAC.IUPACAmbiguousDNA()
        AL=Bio.Align.Generic.Alignment(alpha)
        N=Nexus()
        N.read(infilename)
        for seqname in N.matrix:
            AL.add_sequence(seqname,N.matrix[seqname].tostring())
        return AL
##    def leggialineamento(self, infilename):
##        from Bio.Alphabet import IUPAC
##        import Bio.Align.Generic
##        infile=open(infilename,'r')
##        alpha=IUPAC.IUPACAmbiguousDNA()
##        AL=Bio.Align.Generic.Alignment(alpha)
##        flag=False
##        filetext=infile.read()
##        #print filetext
##        S=re.search('matrix\n', filetext).end()
##        filetext=filetext[S:]
##        S=re.search('\n[\t ]*;', filetext).start()
##        filetext=filetext[:S]
##        L=filetext.split('\n')
##        for line in L:
##             print line
##             pos=re.match('[\t ]+',line).end()
##             L=re.split('[\t ]+',line[pos:],1)
##             
##             if L[0][0]=='[': continue
##             if len(L)==2:
##                if len(L[0])==0:
##                     L=re.split('[\t ]+',L[1],1)
##                if len(L[0])==0:continue
##                #print L
##                L[1]=re.sub(' ','',L[1])
##                AL.add_sequence(L[0],L[1])
##        return AL
##    def leggialineamento(self, infilename):
##        from Bio.Alphabet import IUPAC
##        import Bio.Align.Generic
##        infile=open(infilename,'r')
##        alpha=IUPAC.IUPACAmbiguousDNA()
##        AL=Bio.Align.Generic.Alignment(alpha)
##        flag=False
##        filetext=infile.read()
##        #print filetext
##        S=re.search('matrix\n', filetext).end()
##        filetext=filetext[S:]
##        S=re.search('\n[\t ]*;', filetext).start()
##        filetext=filetext[:S]
##        L=filetext.split('\n')
##        for line in L:
##            line=''.join(re.split('\[.+\]',line))
##            try:
##                name,seq=re.split('[\t ]+',line)
##            except ValueError:
##                continue
##            AL.add_sequence(name,seq)
##        return AL
    def leggialineamento(self, infilename):
        from Bio import AlignIO
        AL=AlignIO.read(infilename,'nexus')
        return AL
##    def multistatistica(self,infilename,pse=0,Basi=0):
##        #calcolo della distribuzione dei motivi di sequenza
##        #AL=self.leggialineamentoNexus(infilename)
##        AL=self.leggialineamento(infilename)
##        #print AL.get_alignment_length()
##        #print AL.get_alignment_length(),'AL'
##        #print "ccicio"
##        for i in range(AL.get_alignment_length()):
##            if self.motivi.has_key(AL.get_column(i)):
##                diff=self.counter+1-len(self.motivi[AL.get_column(i)])
##                #print "s", diff
##                if diff > 0:
##                    self.motivi[AL.get_column(i)]+=diff*[0]
##                #print self.M
##                self.motivi[AL.get_column(i)][self.counter]+=1
##            else:
##                self.M.append((self.counter*[0])+[1])
##                self.motivi[AL.get_column(i)]=self.M[-1]
##                #print self.motivi[AL.get_column(i)], self.counter
##        #print self.M
##        #se vi sono siti invariabili li aggiungo al conteggio
##        #Codici da controllorare
##        if pse>0:
##            B=0
##            for b in self.ordini['basi'][:-1]:
##                Basi[b]=int(round(Basi[b]*pse))
##                #print Basi
##                B+=Basi[b]
##            Basi[self.ordini['basi'][-1]]=int(pse-B)
##            #print pse, B
##            #print Basi
##            for b in self.ordini['basi']:
##                motivo=b*int(self.reader.taxa)
##                if self.motivi.has_key(motivo):
##                    """????che numero devo mettere in get column???"""
##                    diff=self.counter+1-len(self.motivi[motivo])
##                    if diff > 0:
##                        self.motivi[motivo]+=diff*[0] 
##                    #print self.motivi[motivo]
##                    #print 'B',Basi[b]
##                    self.motivi[motivo][self.counter]+=Basi[b]
##                else:
##                    self.M.append((self.counter*[0])+[Basi[b]])
##                    self.motivi[motivo]=self.M[-1]
    def multistatistica(self,infilename,pse=0,Basi=0):
        #calcolo della distribuzione dei motivi di sequenza
        
        AL=self.leggialineamento(infilename)
        
        for i in range(AL.get_alignment_length()):
            if self.motivi.has_key(AL[:,i]):
                self.motivi[AL[:,i]][-1]+=1
            else:
                self.M.append((self.counter)*[0])
                self.M[-1][-1]+=1
                self.motivi[AL[:,i]]=self.M[-1]
                #print self.motivi[AL.get_column(i)], self.counter
        #print self.M
        #se vi sono siti invariabili li aggiungo al conteggio
        #Codici da controllorare
        if pse>0:
            B=0
            for b in self.ordini['basi'][:-1]:
                Basi[b]=int(round(Basi[b]*pse))
                #print Basi
                B+=Basi[b]
            Basi[self.ordini['basi'][-1]]=int(pse-B)
            #print pse, B
            #print Basi
            for b in self.ordini['basi']:
                motivo=b*int(self.reader.taxa)
                if self.motivi.has_key(motivo):
                    """????che numero devo mettere in get column???"""
                    diff=self.counter-len(self.motivi[motivo])
                    if diff > 0:
                        self.motivi[motivo]+=diff*[0] 
                    #print self.motivi[motivo]
                    #print 'B',Basi[b]
                    self.motivi[motivo][-1]+=Basi[b]
                else:
                    self.M.append(((self.counter-1)*[0])+[Basi[b]])
                    self.motivi[motivo]=self.M[-1]
    def preparaArchivioMotivi(self):
        for i in self.motivi.keys():
            diff=self.counter+1-len(self.motivi[i])
            if diff > 0:
                self.motivi[i]+=diff*[0]
            #sys.stdout.write(str(len(self.motivi[i])) )
        #print self.counter, diff
        self.counter+=1
        #print "ciccio", self.counter, self.M[0], self.motivi[i]
    def monostatistica(self, rep=-1):
        N=0
        n=0
        n=sum([x[rep]*log(x[rep]) for x in self.M if x[rep]!=0])
        N=sum([x[rep] for x in self.M])
##        for m in self.M:
##           mm=m[rep]
##           if mm==0:continue
##           n+=mm*log(mm)
##           N+=mm
        #print '\n',N, n
        S=n-N*log(N)
        print S
        return S
    def conteggioBasi(self):
        for m in self.motivi:
            L=len(m)
            break
        self.countBasi=[]
        for i in xrange(L):
            self.countBasi.append(Table('',cases=['A','C','G','T']))
        
        for m in self.motivi:
            for bn in xrange(len(m)):
                try:
                    self.countBasi[bn][m[bn].upper()]+=self.motivi[m]
                except KeyError:
                    self.countBasi[bn]['NA']+=self.motivi[m]
        return self.countBasi
    def Basi3stat(self):
        GC=[]
        GCskew=[]
        ATskew=[]
        for taxa in self.countBasi:
            GC.append((taxa['G']+taxa['C'])/float(taxa['G']+taxa['C']+taxa['A']+taxa['T']))
            GCskew.append((taxa['G']-taxa['C'])/float(taxa['G']+taxa['C']))
            ATskew.append((taxa['A'])/float(taxa['T']+taxa['A']))
        return GC, GCskew, ATskew
