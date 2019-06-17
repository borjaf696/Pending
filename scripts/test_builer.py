#Borja :)
import sys, os
from random import randint

def readReferences(path):
    files = []
    if os.path.isdir(path):
        for r,d,f in os.walk(path):
            for file in f:
                files.append(path+file)
        print 'Ficheros: ', files
    return files

def buildRandomSequences(numSeqs = 20):
    trans = dict({0:'A',1:'C',2:'G',3:'T'})
    seqs = []
    for i in range(numSeqs):
        seq = ''
        length = randint(100, 1000)
        fList = [randint(0,3) for t in range(length)]
        for j in fList:
            seq += trans[j]
        seqs.append(seq)
    return seqs

def readFiles(files):
    references = []
    for file in files:
        print file
        with open(file,'r') as f:
            fullRead = ''
            for i,line in enumerate(f.readlines()):
                if line[0] == '>' and i != 0:
                    references.append(fullRead)
                    fullRead = ''
                elif line[0] != '>':
                    fullRead += line.rstrip()
            references.append(fullRead)
    return references

def addChanges(references, seqs):
    for seq in seqs:
        ref = randint(0, len(references)-1)
        place = randint(0, len(references[ref]))
        references[ref] = references[ref][0:place]+seq+references[ref][place:len(references[ref])]

def writeNewSequences(references, output, seqs):
    with open(output+'output.fasta','w+') as f:
        for i,r in enumerate(references):
            f.write('>Seq'+str(i)+'\n')
            f.write(r+'\n')
    with open(output+'seqs.fasta','w+') as f:
        for i,r in enumerate(seqs):
            f.write('>Seq'+str(i)+'\n')
            f.write(r+'\n')

if __name__=='__main__':
    print 'Running: ', sys.argv[0]
    print 'References: ', sys.argv[1]
    print 'Output: ', sys.argv[2]
    files = readReferences(sys.argv[1])
    seqs = buildRandomSequences()
    references = readFiles(files)
    original = references
    addChanges(references, seqs)
    writeNewSequences(references, sys.argv[2], seqs)