##Opening the output data file produced from the Intermol Conversion Package
filename=input('What is the output data from the Intermol conversion?')
head=[]
data=open(filename,"r")
for k in data:
    head=head+[k]
##Reading the Contents of the file
for i in head:
    if len(i.split())<2:
        continue
    if i.split()[1]=='atoms':
        Natom=int(i.split()[0])
    if i.split()[1]=='bonds':
        Nbond=int(i.split()[0])
    if i.split()[1]=='angles':
        Nangle=int(i.split()[0])
    if i.split()[1]=='dihedrals':
        Ndihedral=int(i.split()[0])
    if i.split()[1]=='impropers':
        Nimproper=int(i.split()[0])
    if i.split()[1]=='atom':
        Natomtype=int(i.split()[0])
    else:
        continue
i=0
while i < len(head):
    if head[i] == 'Atoms\n':
        atoms=head[i:i+2+Natom]
    if head[i] == 'Bond Coeffs\n':
        bondcoeff=head[i:i+2+Nbond]
        i=i+1
    if head[i] == 'Bonds\n':
        bonds=head[i:i+2+Nbond]
        i=i+1
    if head[i] == 'Angle Coeffs\n':
        anglecoeff=head[i:i+2+Nangle]
        i=i+1
    if head[i] == 'Angles\n':
        angles=head[i:i+2+Nangle]
        i=i+1
    if head[i] == 'Dihedral Coeffs\n':
        dihedralcoeff=head[i:i+2+Ndihedral]
        i=i+1
    if head[i] == 'Dihedrals\n':
        dihedrals=head[i:i+2+Ndihedral]
        i=i+1
    if head[i] == 'Improper Coeffs\n':
        impropercoeff=head[i:i+2+Nimproper]
        i=i+1
    if head[i] == 'Impropers\n':
        impropers=head[i:i+2+Nimproper]
        i=i+1
    else:
        i=i+1
##Intermol seperates each bond,angle, dihedral and improper angle individually even if created from the same molecule file in GROMACS. In Lammps it is mor euseful to group these parameters by type so the data file is read and identical parameters are grouped.
diffbondcoeff=[bondcoeff[2]]
for i in bondcoeff[3:]:
    counter=0
    for j in diffbondcoeff:
        if i.split()[1:]!=j.split()[1:]:
            counter=counter+1
        else:
            continue
    if counter==len(diffbondcoeff):
        diffbondcoeff=diffbondcoeff+[i]
    else:
        continue
diffanglecoeff=[anglecoeff[2]]
for i in anglecoeff[3:]:
    counter=0
    for j in diffanglecoeff:
        if i.split()[1:]!=j.split()[1:]:
            counter=counter+1
        else:
            continue
    if counter==len(diffanglecoeff):
        diffanglecoeff=diffanglecoeff+[i]
    else:
        continue
diffdihedralcoeff=[dihedralcoeff[2]]
for i in dihedralcoeff[3:]:
    counter=0
    for j in diffdihedralcoeff:
        if i.split()[1:]!=j.split()[1:]:
            counter=counter+1
        else:
            continue
    if counter==len(diffdihedralcoeff):
        diffdihedralcoeff=diffdihedralcoeff+[i]
    else:
        continue
##The numbering for the interactions in the data file are replaced with their respective interaction type.
nbonds=[]
for i in bonds[2:]:
    temp=i.split()
    counter=1
    for j in diffbondcoeff:
        if bondcoeff[int(temp[1])+1].split()[1:]!=j.split()[1:]:
            counter=counter+1
        else:
            temp[1]=str(counter)
    nbonds=nbonds+[temp]
nangles=[]
for i in angles[2:]:
    temp=i.split()
    counter=1
    for j in diffanglecoeff:
        if anglecoeff[int(temp[1])+1].split()[1:]!=j.split()[1:]:
            counter=counter+1
        else:
            temp[1]=str(counter)
    nangles=nangles+[temp]
ndihedrals=[]
for i in dihedrals[2:]:
    temp=i.split()
    counter=1
    for j in diffdihedralcoeff:
        if dihedralcoeff[int(temp[1])+1].split()[1:]!=j.split()[1:]:
            counter=counter+1
        else:
            temp[1]=str(counter)
    ndihedrals=ndihedrals+[temp]
##Energy units are converted from KJ mol$^{-1}$ to eV used in LAMMPS metal units.
dbondcoeff=[]
counter=1
for i in diffbondcoeff:
    temp=i.split()
    temp[0]=str(counter)
    temp[2]=str(round((float(temp[2])*0.043375),4))
    dbondcoeff=dbondcoeff+[temp]
    counter=counter+1
danglecoeff=[]
counter=1
for i in diffanglecoeff:
    temp=i.split()
    temp[0]=counter
    temp[2]=str(round((float(temp[2])*0.043375),4))
    danglecoeff=danglecoeff+[temp]
    counter=counter+1
ddihedralcoeff=[]
counter=1
for i in diffdihedralcoeff:
    temp=i.split()
    temp[0]=counter
    temp[2]=str(round((float(temp[2])*0.043375),4))
    temp[3]=str(round((float(temp[3])*0.043375),4))
    temp[4]=str(round((float(temp[4])*0.043375),4))
    temp[5]=str(round((float(temp[5])*0.043375),4))
    temp[6]=str(round((float(temp[6])*0.043375),4))
    ddihedralcoeff=ddihedralcoeff+[temp]
    counter=counter+1
##Heading of data file is curated with new parameter types, then written to new data file.
Dimatype=head[12:29]
Intro=['Altered in Python\n']+head[1:8]+[str(Natomtype)+" atom types\n"]+[str(len(dbondcoeff))+" bond types\n"]+[str(len(danglecoeff))+" angle types\n"]+[str(len(ddihedralcoeff))+" dihedral types\n"]+Dimatype
filehandle=open('newdata.lmp', 'w')
for i in Intro:
    filehandle.write(i)
for i in dbondcoeff:
    for j in i:
        filehandle.write('  '+ str(j) )
    filehandle.write('\n')
filehandle.write('\n' + "Angle Coeffs\n" + '\n')
for i in danglecoeff:
    for j in i:
        filehandle.write('  '+ str(j) )
    filehandle.write('\n')
filehandle.write('\n' + "Dihedral Coeffs\n" + '\n')
for i in ddihedralcoeff:
    for j in i:
        filehandle.write('  '+ str(j) )
    filehandle.write('\n')
filehandle.write('\n')    
for i in atoms:
    filehandle.write(i)
filehandle.write('\n' + "Bonds\n" + '\n')
for i in nbonds:
    for j in i:
        filehandle.write('  '+ str(j) )
    filehandle.write('\n')
filehandle.write('\n' + "Angles\n" + '\n')
for i in nangles:
    for j in i:
        filehandle.write('  '+ str(j) )
    filehandle.write('\n')
filehandle.write('\n' + "Dihedrals\n" + '\n')
for i in ndihedrals:
    for j in i:
        filehandle.write('  '+ str(j) )
    filehandle.write('\n')
filehandle.close()
