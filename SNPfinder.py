#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 15:59:28 2015

@author: dominic
"""

###############################################################################
# Blast Functions #############################################################
###############################################################################

def NucleotideBlast(query, subject, output):
    from os import system
    """Using a query from stdin blasts against a subject"""
    system("blastn -query " + query + " -subject " + subject + 
    " -outfmt 5 >> " + output) #Outfmt 5 specifices xml output
    return 

    
def BestHitblastn(query, subject, output): 
    from os import system
    """This function simply runs blast.
    I prefer the syntax of this over NCBI's python module"""
    system("blastn -query " + str(query) + " -subject " + str(subject) + \
    " -outfmt 6 -max_target_seqs 1 > " + str(output)) #Outfmt 6 specifices tsv output
    return 0
    
def AllGenomeBlast(query, GenomeDir, output):
    from os import listdir
    """Uses the previous blast function to blast query against all 
    subject genomes in a directory"""
    for genome in listdir(GenomeDir):
        NucleotideBlast(query, GenomeDir +'/'+ genome, output)
        print "Blast " +query+ " vs " +genome+ " successfull"
    return 0

###############################################################################
# Accessing output functions ##################################################
###############################################################################

def GetBlastRef(inblast):
    """Produces an array from a blast result tsv file"""
    from numpy import array
    from csv import reader
    ref = array([line for line in reader(open(inblast), delimiter='\t')])
    return ref

def GetBlast(blastfile):
    from Bio.Blast import NCBIXML as Nxml
    """In XML format sbjct details (name) are contained in alignment objects
    hit details (sequence, %ID, hitlength etc) are contained in HSP objects
    HSP objects are a sub object of alignment objects (alignment[0].hsps)
    If there are multiple hits for a single sbjct (contig) then multiple HSP
    objects will be stored in the same alignment object.
    It is therefore neccesarry to extract all HSP objects into a list which
    is indexed by another list containing the sbjct names.  This also deals with
    no hit objects recorded in blast XML format.
    """
    hitindex = []
    hitdata = []
    for ii in Nxml.parse(open(blastfile)): # Blast record iterator
        if ii.alignments != []: # Empty lists show no hits and are ignored - This is the sbjct detail level
            for jj in ii.alignments[0].hsps: # Goes a level deeper - This is the hit detail level
                hitindex.append(ii.alignments[0].hit_def) 
                hitdata.append(jj)
    return hitindex, hitdata
    
###############################################################################
# Blast hit extraction and SNP identification #################################
##############################################################################
#

def GetInContigPosition(sbjct_start, sbjct_end, contigcounter):
    if sbjct_start < sbjct_end:
        output = str(int(sbjct_start) + int(contigcounter))
    else:
        output = str(int(sbjct_start) - int(contigcounter))
    return output
            
def ParaOrOrtho(blastfile, singlehitdata, referenceGenome, tmpseqfile, tmpblastfile): 
    import os
    tmpseqfi = open(tmpseqfile, 'w')
    tmpseqfi.write('>tmpfasta\n' + singlehitdata.sbjct)
    tmpseqfi.close()
    BestHitblastn(tmpseqfile, referenceGenome, tmpblastfile)
    os.system('rm ' + tmpseqfile)
    bestblasthit = GetBlastRef(tmpblastfile)
    os.system('rm ' + tmpblastfile)
    return os.path.basename(blastfile) == bestblasthit[0][1]    
    
    
    


def snpsNseqs(singlehitdata, singlehitindex, alignment_end, blastfile, referenceGenome, tmpseqfile, tmpblastfile):
    """In order to produce good sequences for alignments the sequences must be gapped
    to include gaps at the beginning of the sequence"""
    if ParaOrOrtho(blastfile, singlehitdata, referenceGenome, tmpseqfile, tmpblastfile):
        cntgbpcounter = 1
        algnbpcounter = singlehitdata.query_start
        SNPList = []
        OutSeq = []
        for ii in '-'*(singlehitdata.query_start - 1): # How many gaps to add to the start of the sequence
            OutSeq.append(ii)
        for ii,jj in zip(singlehitdata.query, singlehitdata.sbjct):
            if ii == jj: # Nothing interesting
                OutSeq.append(jj)
                algnbpcounter += 1
                cntgbpcounter += 1
            elif ii == '-': # If gap in query (insertion)
                SNPList.append([singlehitindex, ii, jj, algnbpcounter, GetInContigPosition(singlehitdata.sbjct_start, singlehitdata.sbjct_end, cntgbpcounter) ,singlehitdata.query_start, singlehitdata.query_end, singlehitdata.sbjct_start, singlehitdata.sbjct_end]) #Position of insertion in alignment(doesn't appear in alignment) and contig(for coverage lookup)
                cntgbpcounter += 1 # Does not add to alignment counter as doesn't appear in alignment
            elif jj == '-': # if gap in sbjct (deletion)
                SNPList.append([singlehitindex, ii, jj, algnbpcounter, GetInContigPosition(singlehitdata.sbjct_start, singlehitdata.sbjct_end, cntgbpcounter), singlehitdata.query_start, singlehitdata.query_end, singlehitdata.sbjct_start, singlehitdata.sbjct_end])
                OutSeq.append(jj) # gap added to alignment
                algnbpcounter += 1 # Add position of gap in alignment, contig not added because gap is not present in contig
            
            else: # Last possibility, simple SNP
                SNPList.append([singlehitindex, ii, jj, algnbpcounter, GetInContigPosition(singlehitdata.sbjct_start, singlehitdata.sbjct_end, cntgbpcounter), singlehitdata.query_start, singlehitdata.query_end, singlehitdata.sbjct_start, singlehitdata.sbjct_end])
                OutSeq.append(jj)
                algnbpcounter += 1
                cntgbpcounter += 1
        
        if SNPList == []: # Adds entry to contig if no SNPs are found, this allows it to still be used to index with later
            SNPList.append([singlehitindex,0,0,0,0,singlehitdata.query_start, singlehitdata.query_end,singlehitdata.sbjct_start,singlehitdata.sbjct_end])
        for ii in '-' * (int(alignment_end) - singlehitdata.query_end):
            OutSeq.append(ii)
        return OutSeq, SNPList
    return 'Paralog', 'Paralog'
 
           
            
def snpsNseqs4wholeblast(blastfile, referenceGenome, tmpseqfile, tmpblastfile):
    """ Does snpsNseqs for every hit in a blast file"""
    from numpy import array
    import os
    hitcounter = 1
    hitindex, hitdata = GetBlast(blastfile)
    fullSeq = []
    fullSNP = []
    for ii,jj in zip(hitdata, hitindex):
        tmpSeq, tmpSNP = snpsNseqs(ii, jj, os.path.basename(blastfile).split('_')[-1], blastfile, referenceGenome, tmpseqfile, tmpblastfile)
        if tmpSeq != 'Paralog':
            print "hit " + str(hitcounter) + ": Ortholog"
            fullSeq.append(tmpSeq)
            fullSNP.append(tmpSNP)
            hitcounter += 1
        else:
            print "hit " + str(hitcounter) + ": Paralog"
            hitcounter += 1
    return array(fullSeq), fullSNP
    

def writefasta(seqarray, seqindex, outfile):
    handle = open(outfile, 'w')
    for ii, jj in zip(seqarray, seqindex):
        handle.write('>' + jj[0][0] + '\n' + ''.join(ii) + '\n')
    handle.close()
        
def SNPfilter(SNPlist, seqarray, outfile):
    from numpy import array, nonzero, unique
    from csv import writer
    GenomeIndex = array([ii[0][0].split('_')[0] for ii in SNPlist])
    goodSNP = []
    for SNP in SNPlist:
        for jj in SNP:
            if jj[1] != 'N' and jj[2] != 'N':
                lentesttmp = array([(SNPlist[i][0][0].split('_')[0], SNPlist[i][0][0]) for i in nonzero(seqarray[:,int(jj[3]) - 1] != '-')[0]])
                if len(unique(lentesttmp[:,0])) == len(unique(GenomeIndex)):
                    tmp = seqarray[nonzero(GenomeIndex == jj[0].split('_')[0])][:,int(jj[3]) - 1]
                    tmp = tmp[tmp!='-']
                    if len(unique(tmp)) == 1:
                        goodSNP.append(jj)
    f = open(outfile, 'w')
    outSNP = []
    for ii in goodSNP:
        if ii[1] != 0:
            writer(f).writerow(ii)
            outSNP.append(ii)
    f.close()
    return outSNP
                    

def ExtractHitMain(blastfile, reference, outdir, tmpseqfile, tmpblastfile):
    import os
    
        
    seqarray, seqindex = snpsNseqs4wholeblast(blastfile, outdir + '/AllGenes.fasta', tmpseqfile, tmpblastfile)
    print 'blast file ' +os.path.basename(blastfile)+  ' parsed'
    writefasta(seqarray, seqindex, outdir + '/Alignments/' + os.path.basename(blastfile) + '.fasta')
    print 'fasta ' +os.path.basename(blastfile)+ '.fasta written'
    goodSNPs = SNPfilter(seqindex, seqarray, outdir + '/SNPs/' + os.path.basename(blastfile) + '.csv')
    print str(len(goodSNPs)) + ' SNPs identified'
    return goodSNPs








###############################################################################
# Old Versions ################################################################
###############################################################################
##
#
#
#def GenMultifasta(blastfile, outdir):
#    import os
#    os.system('mkdir ' + outdir + '/Alignments')
#    os.system('mkdir ' + outdir + '/SNPs')
#    import os
#    index,blast = GetBlast(blastfile)
#    strt,end = os.path.basename(blastfile).split('_')[-2:]
#    if len(index) < 10:
#        for hitindex,contighits in zip(index,blast):
#            for hit in contighits:
#                aligncounter = 1
#                contigcounter = 1
#                outseq = ''
#                for querybase,sbjctbase in zip(hit.query,hit.sbjct):
#                                    
#                
#                    if querybase != sbjctbase and querybase != '-' and sbjctbase != '-': # If there is a difference which is not an indel
#                    #Write Out SNPs genomecontig, querybase, sbjctbase, in contig base position
#                        open(outdir + '/SNPs/' + os.path.basename(blastfile), 'a').write(hitindex[0].hit_def + ',' + querybase + ',' + sbjctbase +',' + str(aligncounter + hit.query_start-1)+','+ GetInContigPosition(hit.sbjct_start, hit.sbjct_end, contigcounter)+'\n')
#                        outseq += sbjctbase
#                        aligncounter += 1
#                        contigcounter += 1
#                    
#                    
#                    
#                    
#                    elif querybase != sbjctbase and sbjctbase == '-': # If there is a difference which is a deletion
#                        open(outdir + '/SNPs/' + os.path.basename(blastfile), 'a').write(hitindex[0].hit_def + ',' + querybase + ',' + sbjctbase +',' + str(aligncounter + hit.query_start-1)+','+GetInContigPosition(hit.sbjct_start, hit.sbjct_end, contigcounter)+'\n')
#                        outseq += sbjctbase
#                        aligncounter += 1
#                        
#                    
#                    
#                    elif querybase != sbjctbase and querybase == '-': # If there is a difference which is an insertion
#                        open(outdir + '/SNPs/' + os.path.basename(blastfile), 'a').write(hitindex[0].hit_def + ',' + querybase + ',' + sbjctbase +',' + str(aligncounter + hit.query_start-1)+','+GetInContigPosition(hit.sbjct_start, hit.sbjct_end, contigcounter)+'\n')
#                        contigcounter += 1
#                        
#                
#                    elif querybase == sbjctbase:
#                        outseq += sbjctbase
#                        aligncounter += 1#
#                        contigcounter += 1
#                    
#                    else:
#                        print 'HUH?'
#                    
#                open(outdir + '/Alignments/' + os.path.basename(blastfile), 'a').write('>' + hitindex[0].hit_def + '\n' + '-'*(hit.query_start-1) + outseq + ((int(end)-int(strt)+1)-hit.query_end)*'-' +'\n')#
#
#
###############################################################################
# Utility Functions ###########################################################
###############################################################################

def GetFasta(infasta):
    from Bio import SeqIO
    """Opens a fasta file and returns as object"""
    indat = SeqIO.parse(open(infasta), 'fasta')
    seqdat = [line for line in indat]
    if len(seqdat) == 1:
        seqdat = seqdat[0]
    return seqdat
    

    
    
def traverse(o, tree_types=(list, tuple)):
    if isinstance(o, tree_types):
        for value in o:
            for subvalue in traverse(value, tree_types):
                yield subvalue
    else:
        yield o


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
    




###############################################################################
# Split fasta #################################################################
###############################################################################

def splitgenome(Genome, window, outdir):
    from Bio import SeqIO
    outsplit = []
    for contig in GetFasta(Genome):
        counter = 0
        while len(contig) - counter > 1.5*window: # Prevents the last window being smaller than half a window
            tmp = contig[counter:counter+window]
            tmp.id = contig.id + '_' + str(counter+1) + '_' + str(counter + window)
            tmp.description = tmp.id
            outfile = open(outdir + tmp.id, 'w')
            SeqIO.write(tmp, outfile, 'fasta')
            outfile.close()
            counter += window
            outsplit.append(tmp)
            
        tmp = contig[counter:len(contig)]
        tmp.id = contig.id + '_' + str(counter+1) + '_' + str(len(contig))
        tmp.description = tmp.id
        outfile = open(outdir + tmp.id, 'w') # Change dir to user input
        SeqIO.write(tmp, outfile, 'fasta')
        outfile.close()
        outsplit.append(tmp)
    return outsplit

###############################################################################
# Coverage Functions ##########################################################
###############################################################################


def extractcov(Genome, Contig, Pos, tmpfile, extrapos=10):
    """Fast method of getting coverage areas of interest in a file, faster than iterating
    over whole coverage map in python"""
    from csv import reader
    import os
    for ii in range(-int(int(extrapos)/2), int(int(extrapos)/2 + 1)):

    #    os.system("awk 'lines>0 {print; --lines} /contig_" + Contig + "\t" + str(int(Pos)-extrapos/2) + "\t/ {lines=" + str(extrapos/2) + "}' " + Genome)  
    
    
        os.system('cat ' + Genome + ' | grep -P "contig_' + Contig + '\t' + str(int(Pos)+ii) +'\t" >> ' + tmpfile)    
    CntgCov = [line for line in reader(open(tmpfile), delimiter='\t')]
    print len(CntgCov)
    os.system('rm ' + tmpfile)
    return CntgCov


def GetCov(SNPdat, CovDir, outfile, tmpfile, extrapos=10, mincov=1):
    import re
    import os
    from csv import writer, reader
    CovList = []
    for SNP in [line for line in reader(open(SNPdat))]:
        print SNP

        genome = [x for x in os.listdir(CovDir) if re.search(SNP[0].split('_')[0], x)]
        fulltmp = extractcov(CovDir + '/' + genome[0], SNP[0].split('_')[-1], SNP[4], tmpfile, extrapos)
        extrapostmp = []
        for ii in fulltmp:
            if int(ii[1]) != int(SNP[4]):
                extrapostmp.extend([ii[4],ii[5]])
            else:
                if int(ii[3]) >= mincov:
                    actualpostmp = ii
                else:
                    break
        try:
            actualpostmp.extend(extrapostmp)
            CovList.append(actualpostmp)
        except NameError:
            print "SNP coverage too low or zero"
            pass
                
                        

    output = open(outfile, 'w')
    for line in CovList:
        
        if line != []:
            writer(output).writerow(line)
    output.close()
    return CovList












###############################################################################
# Main Function ###############################################################
###############################################################################

def SNPfinder(Gene, GenomeDir, FullReference, OutDir, tmpseqfile, tmpblastfile, CovDir, extrapos, mincov=1):
    import os
    import time
    
    print "Beginning SNPfinder for gene: " + Gene
    print "Blastn " + Gene + " vs subject genomes"
    strttime = time.time()
    AllGenomeBlast(Gene, GenomeDir, OutDir + '/' + os.path.basename(Gene))
    print "Blastn finished, time taken: " + str(time.time()-strttime)
    print "Starting hit extraction and SNP identification"
    midtime = time.time()
    ExtractHitMain(OutDir + '/' + os.path.basename(Gene), FullReference, OutDir, tmpseqfile, tmpblastfile)
    print "Hits extracted, time taken: " + str(time.time() - midtime)
    os.system('rm ' + OutDir + '/' + os.path.basename(Gene))
    print "Starting SNP coverage extraction"
    midtime = time.time()
    GetCov(OutDir + '/SNPs/' + os.path.basename(Gene) + '.csv', CovDir, OutDir + '/SNPCov/' + os.path.basename(Gene) + '.csv'  , OutDir + '/' + os.path.basename(Gene), extrapos, mincov)
    print "SNP coverage extracted, time taken " + str(time.time() - midtime)
    print "Total time taken: " + str(time.time() - strttime)
















if __name__ == '__main__':
    import sys
    import os
    from numpy import array, nonzero
    import multiprocessing as mp
    sysarray = array(sys.argv)
    
    
    # Determining positions of arguments
    
    
    OutDir = sys.argv[nonzero(sysarray == '-o')[0] + 1]
    
    CovDir = sys.argv[nonzero(sysarray == '-c')[0] + 1]
    
    
    # Optional input of -split
    if len(sysarray[sysarray == '-split']) > 0:
        try:
            window = int(sysarray[nonzero(sysarray == '-split')[0] + 1][0])
        except IndexError:
            "Incorrect split size provided, please input -split ####"
    else:
        window = 1000000
        
        
        
    if len(sysarray[sysarray == '-cpus']) > 0:
        try:
            cpus = int(sysarray[nonzero(sysarray == '-cpus')[0] +1][0])
        except IndexError:
            "Incorrect cpus provided, please input -cpus ####"
    else:
        cpus = 1
                
            
    
    # Optional input of -cpus
    
    if len(sysarray[sysarray == '-minc']) > 0:
        try:
            mincov = int(sysarray[nonzero(sysarray == '-minc')[0] + 1][0])
        except IndexError:
            "Incorrect split size provided, please input -split ####"
    else:
        mincov = 1
    
    
    
    
    try:
        RefGenome = sys.argv[nonzero(sysarray == '-r')[0] + 1]
    except IndexError:
        print "Please input complete reference genome fasta with: -r RefFasta"
        quit()
    except TypeError:
        print "Please input complete reference genome fasta with: -r RefFasta"
        quit()
        
    
    # Optional input of additional position to get coverage
    if len(sysarray[sysarray == '-xpos']) > 0:
        try:
            extrapos = int(sysarray[nonzero(sysarray == '-xpos')[0] +1][0])
        except IndexError:
            "Incorrect extra positions provide, plase input -xpos ####"
    else:
        extrapos = 10
    
    
    
    
    #print nonzero(sysarray == '-cpus')[0]
    #print type(sys.argv[nonzero(sysarray == '-cpus')[0] + 1])
     # If -cpus is supplied
    #if nonzero(sysarray == '-cpus')[0] != []: # If -cpus is supplied
    #    print nonzero(sysarray == '-cpus')[0]
    #print sys.argv[nonzero(sysarray == '-cpus')[0] + 1] == int
        
    #try:
    #    cpus = sys.argv[nonzero(sysarray == 'cpus')[0][0] + 1]
    #    print cpus
    #except IndexError:
    #        print "Incorrect cpus provided, please supply number: -cpus ###"
    #else:
    #        print "Incorrect cpus provided, please supply number: -cpus ###"
        
    
    
    GenomeDir = sys.argv[nonzero(sysarray == '-g')[0] + 1]
    
    
    
    
    
    
    
    
    
###############################################################################
# Making output directories ###################################################
###############################################################################
    
    try:
        os.makedirs(OutDir + '/SplitGenome')
    except OSError:
        pass
    try:
        os.mkdir(OutDir + '/Assemblies')
    except OSError:
        pass
    try:
        os.mkdir(OutDir + '/Alignments')
    except OSError:
        pass
    try:
        os.mkdir(OutDir + '/SNPs')
    except OSError:
        pass
    try:
        os.mkdir(OutDir + '/SNPCov')
    except OSError:
        pass
    
    
###############################################################################
#
###############################################################################

    
    
    
    # running split genome and mv fasta header
    try: 
        splitgenome(RefGenome, window, OutDir +'/SplitGenome/')
    except IOError:
        print 'Incorrect Reference Sequence'
        print "Please input complete reference genome fasta with: -r RefFasta"
        quit()
    
    os.system('cat ' + OutDir + '/SplitGenome/* >> ' + OutDir + '/AllGenes.fasta')
    
    
    # running chngHeader on assemblies, ensures inline with assemblies used for
    # coverage mapping.  If original assemblies are used will produce 
    # Identical assemblies from coverage mapping
    # If coverage mapped assemblies are used simply produces the same thing again
    for ii in os.listdir(GenomeDir):
        SetFastaHeader(GenomeDir + ii, OutDir + '/Assemblies/')    
    
    
    
    
    #threading SNPfinder proper
    
    geneins = os.listdir(OutDir + '/SplitGenome')

    while len(geneins) > 0:
        #print cpus
        #print range(cpus)
        processes = [mp.Process(target=SNPfinder, args=(OutDir + '/SplitGenome/' + str(geneins[x]), GenomeDir, RefGenome, OutDir, OutDir + '/tmpseqfile' + str(x), OutDir + '/tmpblastfile' + str(x), CovDir, extrapos, mincov)) for x in range(cpus)]
        #print processes
        #print len(processes)
        geneins = geneins[2:]
        for p in processes:
            p.start()
            
            
        for p in processes:
            p.join()
            
    
        
    

































































