# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 18:40:37 2015

@author: dominic
"""

def smalting(fastq1, fastq2, assembly, outfile):
    import os
    os.system('smalt index -s 1 tmpindex ' + assembly) # Produces Index File
    os.system('smalt map -x -f sam ' + fastq1 + ' ' + fastq2 + ' tmpindex -o ' + outfile) # Produce mapping sam file
    os.system('rm tmpindex')
    return 0

def samming(samfile, outfile):
    import os
    os.system('samtools view -b -S -o tmpbamfile ' + samfile) # Converts Sam to Bam
    os.system('samtools sort tmpbamfile tmpsrtdbam')
    os.system('rm tmpbamfile')
    os.system('samtools mpileup tmpsrtdbam > ' + outfile)
    os.system('rm tmpbamfile')
    return 0
    

def SetFastaHeader(infasta, outdir):
    from Utility import GetFasta
    from Bio import SeqIO
    import re
    import os
    """Changes all headers in fasta to include
    the file name"""
    counter = 1 #Contig Counter used for ID
    fastdat = GetFasta(infasta)
    for sequence in fastdat:
        # os.path.basename used to take filename from user input
        # Genome Name taken from file name using re
        # re returns 'genomename.' indexing then used to remove the last character (.)
        sequence.id = re.match(r'.+\.',  os.path.basename(infasta)).group()[:-1] + '_' + 'contig_'  + str(counter)
        sequence.description = sequence.id
        counter += 1
    outfile = open(outdir + '/' + os.path.basename(infasta), 'w')
    SeqIO.write(fastdat, outfile, 'fasta')
    outfile.close()
    
    


if __name__ == '__main__':
    import sys
    import os
    from numpy import array, nonzero
    
    sysarray = array(sys.argv)
    fastq1 = sys.argv[nonzero(sysarray == '-1')[0] + 1]
    fastq2 = sys.argv[nonzero(sysarray == '-2')[0] + 1]
    outfile = sys.argv[nonzero(sysarray == '-o')[0] + 1]
    assembly = sys.argv[nonzero(sysarray == '-a')[0] + 1]
    
# Make Output Folders #########################################################

    try:
        os.makedirs(outfile + '/Assemblies/')
    except OSError:
        pass
    try:
        os.mkdir(outfile + '/Mappings/')
    except OSError:
        pass
    

    SetFastaHeader(assembly, outfile + '/Assemblies')
    smalting(fastq1,fastq2,assembly,outfile + '/' + os.path.basename(assembly + 'SAM'))
    samming(outfile + '/Mappings' + os.path.basename(assembly) + 'SAM', outfile + '/' + os.path.basename(assembly + 'BAM'))