#!/usr/bin/env python
#
# Simulate introduction of Rotorua fish into Waipa population

"""
*******************************************************************************
*                                                                             *
* Copyright 2012 Anna Friedlander, Elizabeth Heeg, and Peter Ritchie          *
* (anna.fr@gmail.com; elizabeth.heeg@vuw.ac.nz; peter.ritchie@vuw.ac.nz)      *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation, either version 3 of the License, or (at your option)   *
* any later version; see http://www.gnu.org/licenses/                         *
*                                                                             *
*******************************************************************************

Rotorua (codes 1-male; 2-female : 20 fish), and
Waipa (codes 3-male; 4-female : 50 fish) 
fish genotypes read in from FSTAT file

Waipa population cloned (sequentially) to get 1000 fish

1000 replications of the following simulation (each with a new random seed) are run:
  NO MIGRATION (*2):
  * Waipa fish population (1000 fish) undergo *mating scheme* for 100 gens/years
  MIGRATION:
  * Pre-simulation, 20 random Waipa fish removed and 20 Rotorua fish migrated 
    to Waipa pop (total post-migration pop = 1000)
  * Population undergoes *mating scheme* for 100 gens/years
NB: Random-mating = parents chosen randomly with replacement
Pop size kept constant in both conditions

Mating scheme:
  * All individuals <= MAX-AGE (5yrs; age of mortality) cloned to the next gen/year
  * Individuals > MAX-AGE "die" (are not copied to the next gen/year)
  * The remainder of the population is made up by random mating between 
    individuals between MATING-AGE (3yrs) and MAX-AGE (5yrs) (inclusive)

At the end of each replication, pairwise Fst calculated for all pre-evolution 
and post-evolution pops

Populations are output into GENEPOP files

The simulations have the following assumptions:
  Equal males and females make up the founder group (actual sex data used 10m-10f)
  Each individual makes an equal contribution to the next generation (random-mating)
  There is random-mating between the local population and the founder group (random-mating)
  The probability of offspring survival is equal (all offspring survive)
  The resident population in the Tongariro river has even sex ratios (actual sex data used 480m-520f)
   and they all have the same chances of reproductive success. (random-mating)
  The individuals in the founder group stay in the Tongariro river.
  Population size is kept constant
  Overlapping generations
  Maturity at 3 years, and Mortality at 5 years

"""

import simuOpt, os, sys, time, random
from simuPOP import *


def simuVanilla(replications=1000,gens=100,rsize=20,wsize=50,maxAge=5,matingAge=3):

    #read Fstat file in
    #NOTE THIS PROGRAM IS HARD-CODED TO TAKE A SPECIFIC FILE (WITH A PARTICULAR 
    #FORMAT); IT MUST BE ADAPTED FOR USE WITH ANOTHER FILE
    filein = open('ROTO_WPA', 'r')
    info = filein.readlines()
    filein.close()

    #get num subpops, num loci, pop size
    popinfo = info[0].split()
    subpops = popinfo[0]
    numloci = int(popinfo[1])
    popsize = len(info)-numloci-1
    
    #get loci names
    locinam = []
    for i in range(numloci):
        locinam.append(info[i+1].strip())

    #create Waipa and Rotorua populations
    Rotorua = Population(size=[rsize], loci=numloci, infoFields=['age'],
                         lociNames=locinam, subPopNames=['Rotorua'])
    Waipa   = Population(size=[wsize], loci=numloci, infoFields=['age'],
                         lociNames=locinam, subPopNames=['Waipa'])


    #initialise fish genomes and sex (using data from fstat file)
    idx=numloci+1
    sex = MALE
    for i in range(popsize):
        p0 = []; p1 = []
        nums = info[idx+i].split()
        #get genotype
        for j in range(1,numloci+1):
            p0.append(int(nums[j][0:2]))
            p1.append(int(nums[j][2:4]))
        geno = p0 + p1
        #get sex
        if(int(nums[0])%2 == 1): sex = MALE
        else:                    sex = FEMALE
        #set genotype and sex
        if(i<rsize):
           Rotorua.individual(i).setGenotype(geno)
           Rotorua.individual(i).setSex(sex)   
        else:
           Waipa.individual(i-rsize).setGenotype(geno) 
           Waipa.individual(i-rsize).setSex(sex)   


    #randomly initialise individual's age between 0 and maxAge (5)
    #NOTE THAT WAIPA INDIVIDUALS' AGES ARE SET *BEFORE* CLONING
    Rotorua.setIndInfo([random.randint(0, maxAge) for x in range(popsize)], 'age')
    Waipa.setIndInfo([random.randint(0, maxAge) for x in range(popsize)], 'age')

    #three "virtual" sub-pops: age < matingAge 
    #                          matingAge <= age <= maxAge
    #                          age > maxAge 
    Rotorua.setVirtualSplitter(InfoSplitter('age', cutoff=[matingAge, maxAge+0.1]))
    Waipa.setVirtualSplitter(InfoSplitter('age', cutoff=[matingAge, maxAge+0.1]))
 

    #clone the Waipa fish sequentially to grow each population to 1000 fish
    simu=Simulator(Waipa)
    simu.evolve(initOps = [InitSex()], 
                matingScheme = CloneMating(subPopSize=[1000]),
                gen = 1
                )

    #print out the cloned Waipa pop (only needs to be printed once!)
    gpout("WAIPA-CLONED.gen", simu.population(0))

    #set populations: waipa (1000 waipa fish), waipa2 (1000 waipa fish), 
    #"rotowai" (1000 waipa fish to receive 20 rotorua fish)
    pops = simu.population(0).clone()
    pops.addIndFrom(simu.population(0))
    pops.addIndFrom(simu.population(0))
    pops.addIndFrom(Rotorua)

    #open file and write headers
    fileout = open('simuOUT', 'w')
    fileout.write('preWaipa-preWaipa2\tpreWaipa-preRotowai\tpreWaipa2-preRotowai\t')
    fileout.write('postWaipa-postWaipa2\tpostWaipa-postRotowai\tpostWaipa2-postRotowai\t')
    fileout.write('preWaipa-postWaipa\tpreWaipa-postWaipa2\tpreWaipa-postRotowai\t')
    fileout.write('preWaipa2-postWaipa\tpreWaipa2-postWaipa2\tpreWaipa2-postRotowai\t')
    fileout.write('preRotowai-postWaipa\tpreRotowai-postWaipa2\tpreRotowai-postRotowai\n')

    ##### RUN SIMULATION ##########################################################
    #                                                                             #
    # 1000 replications (10 for test-runs) * 100 years (1 year = 1 generation)    #
    #                                                                             #
    # Population size kept constant                                               #
    # All individuals <= MAX-AGE (5yrs; age of mortality) cloned to the next gen  #
    # The remainder of the population is made up by random mating (parents chosen #
    #   randomly with replacement) of individuals >= MATING-AGE (3yrs)            #
    # Pairwise Fst calculated for all pre-evolution and post-evolution pops       #
    #                                                                             #
    ###############################################################################
    for i in range(replications):
        #print repetition number to stdout
        print "rep %s" % (i+1)

        #set random seed
        getRNG().set()
        random.seed(getRNG().seed())

        #remove 20 random Waipa fish from "rotowai"
        pop = pops.clone()
        while (pop.subPopSize(subPop=2) > 980): 
            pop.removeIndividuals(random.randint(2000,3000))

        #add 20 Rotorua fish to "rotowai"
        pop.mergeSubPops(subPops=[2,3])
        pop.setSubPopName('Waipa',0)
        pop.setSubPopName('Waipa2',1)
        pop.setSubPopName('Rotowai',2)

        #print pre-drift rotowai population to FSTAT file
        rwai = pop.extractSubPops(subPops=[2])
        #rwai.setSubPopName('Rotowai',0)
        fname = "ROTOWAI-REP_%s.gen" % (i+1)
        gpout(fname, rwai)

        #calculate pair-wise Fst on pre-evolution Waipa, Waipa2, Rotowai
        stat(pop, structure=ALL_AVAIL, subPops=[0, 1], suffix='_01pre', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[0, 2], suffix='_02pre', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[1, 2], suffix='_12pre', vars='structure_sp')

        #print (pre-evolution) Fst values to file
        fileout.write('%f\t%f\t%f\t' % (pop.dvars().F_st_01pre,
                                        pop.dvars().F_st_02pre,
                                        pop.dvars().F_st_12pre))


        #run simulation: Random-mating for 100 generations. 
        #pairwise Fst calculated at the end of the simulation.
        simu1=Simulator(pop, stealPops=False)
        simu1.evolve(initOps = [
                                PyExec('fstvals01=0'),
                                PyExec('fstvals02=0'),
                                PyExec('fstvals12=0')
                               ],
                    preOps = InfoExec('age += 1'), 
                    matingScheme = HeteroMating(
                        [
                        # age <= maxAge, copy to the next generation
                        CloneMating(subPops=[(x,y) for x in range(3) 
                                                   for y in range(2)], weight=-1),
                        # random mating for individuals in mating ages
                        RandomMating(subPops=[(x,1) for x in range(3)])
                        ],
                        subPopSize=[1000,1000,1000]
                    ),
                    gen=gens,
                    finalOps=  [
                                #calculate and store pair-wise FST
                                Stat(structure=ALL_AVAIL, subPops=[0, 1], suffix='_01'),
                                Stat(structure=ALL_AVAIL, subPops=[0, 2], suffix='_02'),
                                Stat(structure=ALL_AVAIL, subPops=[1, 2], suffix='_12'),
                                PyExec('fstvals01 = F_st_01'),
                                PyExec('fstvals02 = F_st_02'),
                                PyExec('fstvals12 = F_st_12')

                              ]
                    )

        #print genotypes to file (FSTAT format)
        fname = 'ROTO_WPA-REP_%s.gen' % (i+1)
        gpout(fname, simu1.population(0))


        #print (post-evolution) Fst values to file
        fileout.write('%f\t%f\t%f\t' % (simu1.dvars(0).fstvals01,
                                        simu1.dvars(0).fstvals02,
                                        simu1.dvars(0).fstvals12))


        #calculate and print pairwise Fst for pre vs post evolution pops
        pop.addIndFrom(simu1.population(0))
        stat(pop, structure=ALL_AVAIL, subPops=[0, 3], suffix='_03pp', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[0, 4], suffix='_04pp', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[0, 5], suffix='_05pp', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[1, 3], suffix='_13pp', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[1, 4], suffix='_14pp', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[1, 5], suffix='_15pp', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[2, 3], suffix='_23pp', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[2, 4], suffix='_24pp', vars='structure_sp')
        stat(pop, structure=ALL_AVAIL, subPops=[2, 5], suffix='_25pp', vars='structure_sp')
        fileout.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (pop.dvars().F_st_03pp,
                                                                pop.dvars().F_st_04pp,
                                                                pop.dvars().F_st_05pp,
                                                                pop.dvars().F_st_13pp,
                                                                pop.dvars().F_st_14pp,
                                                                pop.dvars().F_st_15pp,
                                                                pop.dvars().F_st_23pp,
                                                                pop.dvars().F_st_24pp,
                                                                pop.dvars().F_st_25pp))
        
    fileout.close()


"""
print genotypes to file (GENEPOP format)

this function adapted from fstatUtil.py by Bo Peng, downloaded from:
http://simupop.sourceforge.net/cookbook/pmwiki.php/Cookbook/FstatUtil
"""
def gpout(fname, pop):
    #get num subpops, num loci
    np = pop.numSubPop()
    nl = pop.totNumLoci()
    loci = range(nl)

    #write to file
    f = open(fname, 'w')
    #write header
    f.write('%s\n' % (fname))
    #write loci names
    for loc in loci:
        f.write(pop.locusName(loc)+'\n')
    #write genotypes
    for sp in range(0, np):
        f.write("POP\n")
        spn = pop.subPopName(sp)
        for ind in range(0, pop.subPopSize(sp)):
            f.write("%s ," % (spn))
            geno = pop.individual(ind+sp*pop.subPopSize(sp)).genotype()
            for al in range(nl):
                a = str(geno[al]).rjust(2,'0')
                b = str(geno[al+nl]).rjust(2,'0')
                f.write(" %s%s" % (a,b))
            f.write('\n')      
    f.close()



if __name__ == '__main__':
    if (len(sys.argv) > 1):
        gen = int(sys.argv[1])
        simuVanilla(gens=gen)

    else: simuVanilla()


