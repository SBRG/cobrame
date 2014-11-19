
# coding: utf-8

## Build a basic ME model

# We will try to build an ME model from the NC_000913.2 Genbank file, the iJO1366 M model, and the complex reconstruction from iJL1650-ME

# In[1]:

from sympy import init_printing, Basic, lambdify
from collections import defaultdict
import re

import cobra.test

from minime import *

init_printing()  # to make things pretty


# In[2]:

me = MEmodel("iJO1366-ME")
me.id = "iJO1366-ME"


#### Start with iJO1366

# In[3]:

iJO1366 = cobra.test.create_test_model(cobra.test.ecoli_pickle)

# add in all the metabolites
for met in iJO1366.metabolites:
    metab=Metabolite(met.id)
    metab.name=met.name
    metab.charge=met.charge
    metab.compartment=met.compartment
    metab.formula=met.formula
    me.add_metabolites(metab)


# In[4]:

for r in iJO1366.reactions:
    # We will add all catalyzed reactions in using the MetabolicReactionData class.
    # Boundary reactions and the biomass reaction are not enzyme catalyzed.
    if not r.boundary and not r.id.startswith("Ec_biomass"):
        if r.gene_reaction_rule != 's0001':
            reaction = MetabolicReactionData(r.id, me)
            reaction._stoichiometry = {met.id: value for met, value in r._metabolites.iteritems()}      
            reaction.lower_bound = r.lower_bound
            reaction.upper_bound = r.upper_bound
        else:
            # Spontaneous reactions added directly
            me.add_reaction(r)
# Boundary reactions added directly.
me.add_reactions([i for i in iJO1366.reactions if i.id.startswith("EX_") and i.lower_bound < 0])
me.add_reactions([i for i in iJO1366.reactions if i.id.startswith("DM_")])


# Also make a dummy reaction

# In[5]:

dummy = MetabolicReactionData("dummy_reaction", me)
dummy.lower_bound = 0
dummy.upper_bound = 1000
dummy._stoichiometry = {}


#### Add Transcription/Translation from genBank

# In[6]:

from Bio import SeqIO
gb_file = SeqIO.read('NC_000913.2.gb', 'gb')
full_seq = str(gb_file.seq)


# In[7]:

def add_transcription_reaction(me_model, TU_name, locus_ids, sequence, update=True):
    """add a transcription reaction"""
    transcription = TranscriptionReaction("transcription_" + TU_name)
    transcription.transcription_data = TranscriptionData(TU_name, me)
    transcription.transcription_data.nucleotide_sequence = sequence
    transcription.transcription_data.RNA_products = {"RNA_" + i for i in locus_ids}
    me.add_reaction(transcription)
    if update:
        transcription.update()


# In[8]:

for feature in gb_file.features:
    if feature.type == "CDS":
        bnum = feature.qualifiers["locus_tag"][0]
        seq = full_seq[feature.location.start:feature.location.end]
        if feature.strand == -1:
            seq = util.dogma.reverse_transcribe(seq)
        add_transcription_reaction(me, "mRNA_" + bnum, {bnum}, seq)
        try:
            amino_acid_sequence = feature.qualifiers["translation"][0]
        except KeyError:
            continue
            #translaion.translation_data.compute_sequence_from_DNA(dna_sequence)
        translation = TranslationReaction("translation_" + bnum)
        me.add_reaction(translation)
        translation.translation_data = TranslationData(bnum, me, "RNA_" + bnum, "protein_" + bnum)
        translation.translation_data.amino_acid_sequence = amino_acid_sequence.replace("U", "C")  # TODO make selenocystine
        translation.update()  # wasteful to update twice
        
    elif feature.type == "rRNA":
        bnum = feature.qualifiers["locus_tag"][0]
        seq = full_seq[feature.location.start:feature.location.end]
        if feature.strand == -1:
            seq = util.dogma.reverse_transcribe(seq)
        add_transcription_reaction(me, "rRNA_" + bnum, {bnum}, seq)
    

    elif feature.type == "tRNA":
        # TODO account for modifications
        bnum = feature.qualifiers["locus_tag"][0]
        seq = full_seq[feature.location.start:feature.location.end]
        if feature.strand == -1:
            seq = util.dogma.reverse_transcribe(seq)
        add_transcription_reaction(me, "tRNA_" + bnum, {bnum}, seq)    


# Build ribosome and RNA Polylmerase

# In[9]:

ribosome_formation = Reaction("ribosome_formation")
ribosome_components = {}
for i in ["b3851", "b3854", "b3855"]:
    ribosome_components[me.metabolites.get_by_id("RNA_" + i)] = -1
# 30S Listed as [rpsA -rpsU], sra 
# [rplA-rplF],  rplI, [rplK-rplY],
# [rpmA-rpmJ]

for i in ["b0911", "b0169", "b3314", "b3296", "b3303", "b4200", "b3341",
          "b3306", "b3230", "b3321", "b3297", "b3342", "b3298", "b3307", 
          "b3165", "b2609", "b3311", "b4202", "b3316", "b0023", "b3065",
          "b1480"]:
    ribosome_components[me.metabolites.get_by_id("protein_" + i)] = -1
    
# 50S listed as [rplA-rplF],rplJ, rplI, rplK [rplM-rplY], [rpmA-rpmJ]    
for i in ["b3984", "b3317", "b3320", "b3319", "b3308", "b3305", 
          "b3958", "b4203", "b3983", "b3231", "b3310", "b3301",
          "b3313", "b3294", "b3304", "b2606", "b1716", "b3186",
          "b3315", "b3318", "b3309", "b2185", 
          "b3185", "b3637", "b3312", "b3302", "b3936", "b1089", 
          "b3636", "b3703", "b1717", "b3299"]:
    ribosome_components[me.metabolites.get_by_id("protein_" + i)] = -1
# [rplJ, 2(2[rplL])]
ribosome_components[me.metabolites.get_by_id("protein_" + "b3986")] = -4
ribosome_components[Ribosome("ribosome")] = 1
ribosome_formation.add_metabolites(ribosome_components)
me.add_reaction(ribosome_formation)


# In[10]:

RNAP_formation=Reaction("RNAP_formation")
RNAP_components = {}
# Core RNA Polymerase Enzyme
for i in {"b3295" : "rpoA", "b3988" : "rpoC", "b3987" : "rpoB"}:
    if i == "b3295": RNAP_components[me.metabolites.get_by_id("protein_" + i)] = -2
    else: RNAP_components[me.metabolites.get_by_id("protein_" + i)] = -1
RNAP_components[RNAP("RNA_Polymerase")] = 1
RNAP_formation.add_metabolites(RNAP_components)
me.add_reaction(RNAP_formation)


# In[11]:

for r in me.reactions:
    if isinstance(r, TranslationReaction):
        r.update()
    if isinstance(r, TranscriptionReaction):
        r.update()


# Add a dummy protein in as well

# In[12]:

dna_sequence = "ATG" + "TTT" * 5 + "TAT"*5+ "ACG"*5 + "GAT" *5 + "AGT"*5 + "TAA"
add_transcription_reaction(me, "dummy", {"dummy"}, dna_sequence)
me.add_metabolites(TranslatedGene("protein_" + "dummy"))
translation = TranslationReaction("translation_" + "dummy")
me.add_reaction(translation)
translation.translation_data = TranslationData("dummy", me, "RNA_dummy", "protein_dummy")
translation.translation_data.compute_sequence_from_DNA(dna_sequence)
translation.update()

complex_data = ComplexData("CPLX_dummy", me)
complex_data._stoichiometry = {}
complex_data._stoichiometry["protein_" + "dummy"] = 1


#### Add in the complex formation

# In[13]:

ME_complex = open('protein_complexes.txt')
ME_complex_dict={}

for line in ME_complex:
    line=line.rstrip('\tM_protein_recon\n')
    line=line.rstrip('\t2011_Updated_E_recon\n')
    line=re.split('\t| AND |',line)
    ME_complex_dict[line[0]]=line[2:]

for cplx, value in ME_complex_dict.iteritems():
    complex_data = ComplexData(cplx, me)
    complex_data._stoichiometry = {}
    for gene in value:
        stoichiometry = gene[6]
        bnum = gene[0:5]
        try:
            complex_data._stoichiometry["protein_" + bnum] = float(stoichiometry)
        except:
            complex_data._stoichiometry["protein_" + bnum] = float(1)
ME_complex.close()


# Only use teddy's complex information file

# In[14]:

enzRxn=open('enzyme_reaction_association.txt','r')
rxnToModCplxDict={}
for line in enzRxn:
    line=line.rstrip('\n')
    line=re.split('\t| OR ',line)
    rxnToModCplxDict[line[0]]=line[1:]


# In[15]:

enzMod=open('protein_modification.txt','r')
modToCompDict={}
for line in enzMod:
    line=line.rstrip('\tM_protein_recon\n')
    line=re.split('\t| AND |',line)
    modToCompDict[line[0]]=line[1:]


# In[16]:

rxnToCompDict={}
for reaction in me.metabolic_reaction_data:
    compList=[]
    try:
        modcomplex=rxnToModCplxDict[reaction.id]
        for modcplx in modcomplex:
            try:
                compList.append(modToCompDict[modcplx])
            except KeyError:
                compList.append(modcplx)
        rxnToCompDict[reaction]=compList
    except KeyError:
        rxnToCompDict[reaction]=None


# Associate complexes with reactions

# In[17]:

for reaction, complexes in rxnToCompDict.items():
    # use a dummy complex for reactions with no assigned complex
    if complexes is None:
        complexes = ["CPLX_dummy"]
    for cpx in complexes:
        if isinstance(cpx, list):
            cpx = cpx[0]
        try:
            complex_data = me.complex_data.get_by_id(cpx)
        except KeyError:
            print cpx, "not found"
            continue
        if reaction.upper_bound > 0:
            r = MetabolicReaction(reaction.id + "_FWD_" + cpx)
            r.keff = 65
            r.complex_data = complex_data
            r.metabolic_reaction_data = reaction
            
            # Keep complexes with multiple modifications from
            # being added twice to me
            try:
                me.add_reaction(r)
            except:
                continue
            r.update()
        if reaction.lower_bound < 0:
            r = MetabolicReaction(reaction.id + "_REV_" + cpx)
            r.keff = 65.
            r.complex_data = complex_data
            r.metabolic_reaction_data = reaction
            r.reverse = True
            
            # Keep complexes with multiple modifications from
            # being added twice to me
            try:
                me.add_reaction(r)
            except:
                continue
            r.update()


# In[18]:

# associate the dummy reaction with a dummy complex as well
r = MetabolicReaction("dummy")
r.objective_coefficient = 1
r.keff = 1
r.complex_data = me.complex_data.get_by_id("CPLX_dummy")
r.metabolic_reaction_data = me.metabolic_reaction_data.get_by_id("dummy_reaction")
me.add_reaction(r)
r.update()


# Remove unused protein and mRNA to make the model solve faster

# In[19]:

for p in me.metabolites.query("protein"):
    if len(p._reaction) == 1:
        list(p._reaction)[0].delete(remove_orphans=True)
for m in me.metabolites.query("RNA"):
    if len(m._reaction) == 1:
        list(m._reaction)[0].delete(remove_orphans=True)


# This gives the total number of genes included

# In[20]:

len(me.reactions.query("transcription"))


### Solve

# In[21]:

from time import time

from cobra.solvers.soplex import Soplex


# In[22]:

compiled_expressions = {}
for i, r in enumerate(me.reactions):
    for met, stoic in r._metabolites.iteritems():
        if isinstance(stoic, Basic):
            compiled_expressions[(me.metabolites.index(met), i)] = lambdify(mu, stoic)


# In[23]:

lp = Soplex.create_problem(me)


# In[24]:

def solve_lp(lp, mu):
    for index, expr in compiled_expressions.iteritems():
        lp.change_coefficient(index[0], index[1], expr(mu))
    lp.solve_problem()
    status = lp.get_status()
    if status == "optimal":
        return (status, lp.get_objective_value())
    else:
        return (status,)


# Binary search and Newton's method code

# In[25]:

def try_mu(lp, mu, feasible_mu=[], infeasible_mu=[], productions=[]):
    r = solve_lp(lp, mu)
    if r[0] == "optimal" and r[1] > 0:
        print mu, r[1]
        feasible_mu.append(mu)
        productions.append(r[1])
        return True
    else:
        print mu
        infeasible_mu.append(mu)
        return False

def binary_search(lp, feasible_mu, infeasible_mu, productions, fraction=0.5):
    mu = fraction * infeasible_mu[-1] + (1 - fraction) * feasible_mu[-1]
    return try_mu(lp, mu, feasible_mu, infeasible_mu, productions)

def newtons_method(lp, feasible_mu, infeasible_mu, productions):
    m1, m2 = feasible_mu[-2:]
    p1, p2 = productions[-2:]
    mu = m1 - p1 * (m2 - m1) / (p2 - p1) if p2 != p1 else infeasible_mu[-1] * 1.5  # detect divide by 0
    if mu > infeasible_mu[-1]:
        print "guess will be infeasible"
        return binary_search(lp, feasible_mu, infeasible_mu, productions, fraction=0.25)
    return try_mu(lp, mu, feasible_mu, infeasible_mu, productions)


# In[26]:

feasible_mu = []
infeasible_mu = []
productions = []
start = time()
# try the edges of binary search
if not try_mu(lp, 0, feasible_mu, infeasible_mu, productions):
    print "0 needs to be feasible"
max_mu = 2
while try_mu(lp, max_mu, feasible_mu, infeasible_mu, productions):
    max_mu += 1
while infeasible_mu[-1] - feasible_mu[-1] > 1e-9:
    binary_search(lp, feasible_mu, infeasible_mu, productions, fraction=0.5)
print "completed in %.1f seconds and %d iterations" % (time() - start, len(feasible_mu) + len(infeasible_mu))

