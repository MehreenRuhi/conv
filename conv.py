#Finds convergent amino acids between a species (or group) and another species (or group of species)



import os
cwd=os.getcwd()
import re
import sys
fname=sys.argv[1]



def consensus(group,dictionary):
    cons=[]
    x=group[group.find("(")+1:group.find(")")]
    f=x.split(',')
    for i in f:
        cons.append(dictionary[i])
    return consensoos(cons)
   
def namechange(group):
    import re
    x=group[group.find("(")+1:group.find(")")]
    slashes= x.replace(',','/')
    b=re.sub("([\(\[]).*?([\)\]])", "\g<1>\g<2>", group)
    d= re.sub('[()]','',b)
    i=slashes+str(',')+d
    return i
   
def consensoos(cons):
    same=''
    for i in range(len(cons[0])):
        consensus=[]
        for j in range(len(cons)):
            consensus+=cons[j][i]
        if '-' in consensus:
            same+= '-'  
        elif len(set(consensus))==1:
            same+=str(consensus[1]) 
        else:
            same+= 'X'
    return same

def outgroupconsensus(group, dictionary):
    import re
    cons=[]
    b=re.sub("([\(\[]).*?([\)\]])", "\g<1>\g<2>", group)
    d= re.sub('[()]','',b)
    x=d.split(',')
    for i in x:
        if i !='':
            cons.append(dictionary[i])
    return consensoos(cons)

def diverged(f,r):
    s=''
    for i in range(len(f)):
        if r[i]=='-' or f[i]=='-':
            s+= '-'
        elif r[i]=='X' or f[i]=='X':
            s+= 'X'
        elif r[i]==f[i]:
            s+= 'U'
        else:
            s+= str(f[i])
    return s
    
     
def diverge(clade,species):
    divergedspecies=[]
    divergedsequence=[]
    for index in range(len(species)):       
        b=[x for i, x in enumerate(clade) if i!=index]
        consensus= consensoos(b)
        compto= clade[index]
        divtest=diverged(compto,consensus)
        divergedsequence.append(divtest)
    return divergedsequence

def same(x,y):
    splitx=x.split(',')
    splity=y.split(',')
    for i in splitx:
        if i in splity:
            return True
            
def makegroup(x):
    c=0
    species=[]
    retlist=[]
    for i in range(len(x)):
        species=[]
        species.append(x[c])
        species.append(x)
        c+=1
        speciesend=species[1:]
        specbeg=species[0]
        spec=[item for sublist in speciesend for item in sublist]
        d= [specbeg]+spec
        yes= [",".join(d)]
        retlist.append(yes)           
    return retlist

def countpossible(g,j):
    c=0
    for i in range(len(g)):
        if g[i] !='X' and j[i]!='X'and g[i]!='-' and j[i]!='-':
            c+=1
    return c

def converge(g,j):
    s=''
    for i in range(len(g)):
        if g[i]=='-' or j[i]=='-':
            s+= '-'
        elif g[i]==j[i]:
            s+= str(g[i])
        else:
            s+= 'X'
    return s


def writetofile(x,xx, iden,out):
	results2=open(cwd+'/whole_'+out,'a')   
	results2.write(str(xx))
	results2.write('=possible_converge')
	results2.write(',')
	results2.write(iden)
	results2.write('\n') 
	chars=set('ARNDCEQGHILKMFPSTWYV')
	if any((c in chars) for c in x):
		results=open(cwd+'/conv_'+out,'a')
		results.write(str(x))
		results.write(',')
		results.write(str(x.count('-')))
		results.write('=missing')
		results.write(',')
		results.write(str(len(x)))
		results.write('=length')
		results.write(',')
		results.write(str(x.count('X')))
		results.write('=not_diverged_or_not_converged')
		results.write(',')
		results.write(str(xx))
		results.write('=possible_converge')
		results.write('\n')
		results.write(iden)
		results.write('\n')

        
def convergence(results,species,iden,out,sconv):
	shortspec=species
    	shortres=results
	if sconv=='0010':
    		conv=converge(shortres[0][0],shortres[1][0])   
    		possible=countpossible(shortres[0][0],shortres[1][0])
   	elif sconv=='01':
                conv=converge(shortres[0],shortres[1])
                possible=countpossible(shortres[0],shortres[1])
	elif sconv=='010':
                conv=converge(shortres[0],shortres[1][0])
                possible=countpossible(shortres[0],shortres[1][0])
	elif sconv=='001':
                conv=converge(shortres[0][0],shortres[1])
                possible=countpossible(shortres[0][0],shortres[1])
    	writetofile(conv,possible,iden,out)

def finddivergence2(grouping,dictionary,iden,out):
    results=[]
    species=[]
    for i in grouping:
       clade=[]
       div=[]
       if '(' in i:
           f= consensus(i,dictionary)
           d=namechange(i)
           r=outgroupconsensus(i,dictionary)
           divgd=diverged(f,r)
           results.append(divgd)
           species.append(d)
       else:
           x=i.rstrip(',').split(',')
           for j in x:
               clade.append(dictionary[j.rstrip('*').lstrip(' ')])
           g=diverge(clade,x)
           results.append(g)
           mg=makegroup(x)
           species.append(mg)
    convergence(results,species,iden,out)

def finddivergence(grouping,dictionary,iden,out):
    	results=[]
    	species=[]	
    	sconv=''
	if '(' in grouping[0] and '(' in grouping[1]:
		sconv='01'
	elif '(' in grouping[0]:
		sconv='010'
	elif ')' in grouping[1]:
		sconv='001'
	else:
		sconv='0010'       
	for i in grouping:       
		clade=[]
       		div=[]
       		if '(' in i:
        		f= consensus(i,dictionary)
			d=namechange(i)
           		r=outgroupconsensus(i,dictionary)
           		divgd=diverged(f,r)
           		results.append(divgd)
           		species.append(d)
       		else:
           		x=i.rstrip(',').split(',')
           		for j in x:
               			clade.append(dictionary[j.rstrip('*').lstrip(' ')])
           		g=diverge(clade,x)
           		results.append(g)
           		mg=makegroup(x)
           		species.append(mg)
	convergence(results,species,iden,out,sconv)


codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

gene=[]
geneedited=[]
species=[]
dictionary={}
dictionary2={}
groups2=[sys.argv[2],sys.argv[3]]

groupspart1=groups2[0].split(',!')[0].replace('_',',').rstrip("'").lstrip("'")
groupspart2=groups2[1].split(',!')[0].replace('_',',').rstrip("'").lstrip("'")

f=open(cwd+'/'+fname,'r')   

groups=[groupspart1,groupspart2]
out=groups[0]+'_'+groups[1]
iden=fname.split('.')[0].split('/')[1]

for i in f:
	if i[0]!='>':
		prot=''
		cds=i.rstrip('\n')
		
		for n in range(0,len(cds),3):
			
			if cds[n:n+3] in codontable:
				prot+= codontable[cds[n:n+3]]
			elif '-' in cds[n:n+3]:
				prot+='-'
			elif 'N' in cds[n:n+3]:
				prot+='-' 
			else:
				prot+='X'
		gene.append(prot)
	else:
		if '_' in i:
			species.append(i.split('_')[1].rstrip('\n'))

		else:
			species.append(i.split('.')[1].rstrip('\n'))
dictionary = dict(zip(species, gene))


editedmegseq=[]
megseq=dictionary["Megaladapis2xMasked"]

for eachspecies in species:
	editedprot=''
	oldprotseq=dictionary[eachspecies]
	for eachamino in range(len(megseq)):
		if megseq[eachamino]=='-':
			editedprot+='-'
		else:
			editedprot+= oldprotseq[eachamino]
	geneedited.append(editedprot)
dictionary2=dict(zip(species,geneedited))
grouptotlist=[]

for group1or2 in groups:
	allspeciesalone=group1or2.replace(')',',').replace('(','').split(',')
	for eachone in allspeciesalone:
		grouptotlist.append(eachone)
if len(grouptotlist)==len(set(grouptotlist)):

	finddivergence(groups, dictionary2,iden,out)
else:
	pass

f.close()




