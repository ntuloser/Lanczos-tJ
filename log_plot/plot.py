import pylab

Energy=[]
f=open('log_opt_by_Fi','rw')
for line in f:
    line=line[:-1]
    Lst=line.split(',')
    a=Lst[0]
    Energy.append(float( Lst[0][3:] ))

x=range(0,len(Energy))
pylab.plot(x[:], Energy[:], '.-', lw=1.5,color='r')

Energy=[]
f=open('log','rw')
for line in f:
    line=line[:-1]
    Lst=line.split(',')
    a=Lst[0]
    Energy.append(float( Lst[0][3:] ))
x=range(0,len(Energy))
pylab.plot(x[:], Energy[:], '.-', lw=1.5,color='g')
    
pylab.show()
