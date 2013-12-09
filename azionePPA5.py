from PPA4 import readmcmc, simulatorMol, findThinning
import sys
import math

com={'-b':'1', '-Motivi':False}
count=1
key=None
for i in sys.argv:
    if key!=None:
        com[key]=i
        key=None
    if i[0]=="-":
        key=i
spiegazione=""" 
 -i radice comune dei nomi dei file di uscita di mrbayes
 -ii file di input di mb
 -mb eseguibile per mb
 -e eseguibile per il simulatore
 -s  numeri di campioni della distribuzione a posterori
 -b burnin 
 -Motivi True/False :scrivi su file la distribuzione dei motivi 
"""
if len(com)<6:
    print spiegazione
#Produci le distribuzioni predittive
MB=readmcmc(infilename=com['-i'],nexusfilename=com['-ii'], pathmb=com['-mb'])
sim=simulatorMol(reader=MB,pathevolver=com['-e'])
t=findThinning(com['-ii'], com['-b'],com['-s'] )
thin=range(t)
thin.reverse()
burn=range(int(com['-b'])+1)
while 1:
    if 0==burn.pop(): break
    sim.skipline()
while 1:
    T=thin.pop(0)
    thin.append(T)
    if T==0:
        segnale=sim.writeout()
        if segnale=='EOF':break
        
    else: sim.skipline()

Mexpc=sim.statdist
print "Observed Value"
sim.preparaArchivioMotivi()
sim.multistatistica(infilename=com['-ii'],pse=0,Basi=0)
Mobs=sim.monostatistica()
#sim.conteggioBasi()
#sim.puliscimotivi()
motivi=sim.M
#Calcolo delle statistiche
##SA=0
##SAv=0
##SAm=0
##
##print len(motivi),'motivi'
##for mot in motivi:
##    m=0
##    for tom in mot[:-1]:
##        m+=tom
##       
##    print (len(mot))
##    m=m/(len(mot)-1)
##    m=math.pow(m-mot[-1],2)
##    SAm+=m
##    v=0
##    for tom in mot[:-1]:
##        v+=math.pow(tom-m,2)
##    v/(len(mot)-1)
##    SAv+=v
##    SA+=m+v
##SA=math.pow(SA,0.5)
P=0
meanE=0
devEO=0
for i in Mexpc:
    if i>=Mobs: P+=1
    meanE+=i
    devEO+=math.pow(i-Mobs,2)

meanE=meanE/len(Mexpc)
varE=0
for i in Mexpc:
    varE+=math.pow(i-meanE,2)

varE=varE/len(Mexpc)
devEO=devEO/len(Mexpc)
#statistica bolback
P=P/float(len(Mexpc))
#statistica L
Ls=math.pow(devEO+varE,0.5)
L =math.pow(math.pow(meanE-Mobs,2)+varE,0.5)
#uscita="\n".join(['SA =\t'+str(SA),'SAm=\t'+str(SAm),'SAv=\t'+str(SAv),'L =\t'+str(L),'Ls =\t'+str(Ls),'varSim=\t'+str(varE),'P =\t'+str(P),'meanE=\t'+str(meanE),'Mobs=\t'+str(Mobs)+'\n'])
output=open(com['-i']+'.checkmodel','w')
#output.write('SA =\t'+str(SA))
#output.write('SAm=\t'+str(SAm))
#output.write('SAv=\t'+str(SAv))
##output.write('L =\t'+str(L)+'\n')
##output.write('Ls =\t'+str(Ls)+'\n')
##output.write('Bolback stat')
##
##output.write('P =\t'+str(P)+'\n')
##output.write('meanE=\t'+str(meanE)+'\n')
##output.write('Mobs=\t'+str(Mobs)+'\n')
##output.write('varSim=\t'+str(varE)+'\n')

output.write('<res>\n')
output.write('\t<L description="Ibrahim statistics">'+str(L)+'</L>\n')
output.write('\t<Ls description="Ibrahim statistics">'+str(Ls)+'</Ls>\n')
output.write('\t<P description="Bolback statistics">'+str(P)+'</P>\n')
output.write('\t<meanE description="Mean Expected Value of MaxLikelihood">'+str(meanE)+'</meanE >\n')
output.write('\t<Mobs description="Observed MaxLikelihood">'+str(Mobs)+'</Mobs>\n')
output.write('\t<varSim description="Variance of Expected MaxLikelihood">'+str(varE)+'</varSim>\n')
output.write('\t<ExpDist description="Distribution of ExpectedValue">'+",".join(map(str,Mexpc))+'</ExpDist>\n')
output.write('</res>')

#output.write(uscita)
if com['-Motivi']:
    uscitafile=open('distM.txt','w')
    KK=sim.motivi.keys()
    KK.sort()
    for i in KK:
        uscitafile.write(str(i)+' ')
        for j in sim.motivi[i]:
            uscitafile.write(str(j)+' ')
        uscitafile.write('\n')
