#!/usr/bin/env python
# coding: utf-8

# In[1]:


import Bio
import pydna
import pydna.dseq
import Bio.SeqIO
import Bio.Restriction
import glob
import re


# In[2]:


def reducedPlasmidNamesToCanonical(splitInput):
    canonicalNames = []
    for s in splitInput:
        result = re.search("\(([0-9][0-9][0-9])\)",s)
        if result != None:
            canonicalNames.append("pYTK"+result.group(1))
        result = re.search("(pYTK[0-9][0-9][0-9])",s)
        if result != None:
            canonicalNames.append(result.group(1))
        result = re.search("pJKW[ _-]{0,2}([0-9]{1,6})",s)
        if result != None:
            canonicalNames.append("pJKW "+result.group(1))
    return(canonicalNames)


# In[3]:


def extractPlasmidSeqsAndFeatures(plasmidNames,files):
    featuresToApply = dict() ##Key is the sequence (uppered) of the feature, value is the feature itself
    plasmid_Seqs = dict()
    assert len(files) > 0
    for p in plasmidNames:
        for f in files:
            if p in f:
                #print(f)
                bioRecord = list(Bio.SeqIO.parse(f,"gb"))[0]
                pyDNASeq = pydna.dseq.Dseq(bioRecord.seq,bioRecord.seq.reverse_complement(),ovhg=0,circular=True)
                plasmid_Seqs[f] = pyDNASeq

                ##Store the features for reapplication later
                for t in bioRecord.features:
                    if t.location.strand == 1:
                        featSeq = str(bioRecord.seq[t.location.start:t.location.end]).upper()         
                    elif t.location.strand == -1:
                        featSeq = str(bioRecord.seq[t.location.start:t.location.end].reverse_complement()).upper()
                    featuresToApply[featSeq] = t
                ##Did the below just to show that all the annotation gets written back correctly.
                ##Bio.SeqIO.write(lSeq, "example"+lSeq.name+".gbk", "gb")
                ##Thinking of ways to transfer the annotation to the resulting sequence...
                ##(1) Just export all the annotations as sequence/annotation pairs, then reannotate.
    return plasmid_Seqs,featuresToApply


# In[4]:


def reapply_features_to_record(record,featuresDict):
    for k in featuresDict.keys():
        v = featuresDict[k]
        result = record.seq.upper().find(k)
        if result != -1:
            #print("found this +1 strand feature at index:",result)
            #print(v)
            span = v.location.end - v.location.start
            my_start_pos = Bio.SeqFeature.ExactPosition(result)
            my_end_pos = Bio.SeqFeature.ExactPosition(result+span)
            my_feature_location = Bio.SeqFeature.FeatureLocation(my_start_pos,my_end_pos)
            my_feature = Bio.SeqFeature.SeqFeature(my_feature_location,type=v.type,strand=1,qualifiers=v.qualifiers)
            record.features.append(my_feature)
                
        result_rc = record.seq.upper().reverse_complement().find(k)
        result_rc_fixed = len(record.seq)-result_rc
        if result_rc != -1:
            #print("found this -1 strand feature at index:",result_rc)
            #print(v)
            span = v.location.end - v.location.start
            my_start_pos = Bio.SeqFeature.ExactPosition(result_rc_fixed-span)
            my_end_pos = Bio.SeqFeature.ExactPosition(result_rc_fixed)
            my_feature_location = Bio.SeqFeature.FeatureLocation(my_start_pos,my_end_pos)
            my_feature = Bio.SeqFeature.SeqFeature(my_feature_location,type=v.type,strand=-1,qualifiers=v.qualifiers)
            record.features.append(my_feature)
            
    return record


# In[5]:


def digestBsaI(plasmid_Seqs,enzymeName):
    fragments = dict()
    for p in plasmid_Seqs:
        if enzymeName == "BsaI":
            cutSeqs = list(plasmid_Seqs[p].cut(Bio.Restriction.BsaI))
        elif enzymeName == "BsmBI":
            cutSeqs = list(plasmid_Seqs[p].cut(Bio.Restriction.BsmBI))
        ##print(p,len(cutSeqs),cutSeqs)
        fragments[p] = cutSeqs
    return fragments


# In[6]:


def ligate_fragments(theFrags):
    
    ##print("Input fragments num:",len(theFrags))
    all_frags = theFrags

    all_products = []
    for i in range(0,len(all_frags)):
        
        ##Also allow intramolecular ligation
        five_prime_overhang = Bio.Seq.Seq(all_frags[i].five_prime_end()[1].upper())
        three_prime_overhang = Bio.Seq.Seq(all_frags[i].three_prime_end()[1].upper())
        if five_prime_overhang == three_prime_overhang.reverse_complement():
            circSeq = all_frags[i].looped()
            all_products.append(circSeq)

        ##Do the intermolecular ligations
        for j in range(0,len(all_frags)):
            result = ""
            product = None
            
            try:
                product = all_frags[i]+all_frags[j]
                all_products.append(product)
                result += "1_worked"
            except Exception as e:
                pass
            try:
                product = all_frags[i]+all_frags[j].reverse_complement()
                all_products.append(product)
                result += "2_worked"
            except Exception as e:
                pass
            try:
                product = all_frags[i].reverse_complement()+all_frags[j]
                all_products.append(product)
                result += "3_worked"
            except Exception as e:
                pass
            try:
                product = all_frags[i].reverse_complement()+all_frags[j].reverse_complement()
                all_products.append(product)
                result += "4_worked"
            except Exception as e:
                pass
            try:
                product = all_frags[j]+all_frags[i]
                all_products.append(product)
                result += "5_worked"
            except Exception as e:
                pass
            try:
                product = all_frags[j]+all_frags[i].reverse_complement()
                all_products.append(product)
                result += "6_worked"
            except Exception as e:
                pass
            try:
                product = all_frags[j].reverse_complement()+all_frags[i]
                all_products.append(product)
                result += "7_worked"
            except Exception as e:
                pass
            try:
                product = all_frags[j].reverse_complement()+all_frags[i].reverse_complement()
                all_products.append(product)
                result += "8_worked"
            except Exception as e:
                pass
            ##print(i,j)
            #if result != "":
            #    print(result,len(product))
            #else: 
            #    print("None")
    ##print("Total ligation products:"+str(len(all_products)))
    dereplicated_products = []

    for i in all_products:
        ##This is how to do the "if not in" conditional
        itemIn = False
        for j in dereplicated_products:
            if i == j or i == j.reverse_complement():
                itemIn = True
        if itemIn == False:
            dereplicated_products.append(i)

    ##print("Dereplicated ligation products:"+str(len(dereplicated_products)))
    return(dereplicated_products)


# In[7]:


def subtract_frags_by_seq(frags):
    passed_frags = []
    GFP_30bp = "ataccctggtaaaccgcattgagctgaaag".upper()
    CamR_35bp = "agaagttgtccatattggccacgtttaaatcaaaa".upper()
    seq_to_screen = [GFP_30bp,CamR_35bp]
    for f in frags:
        all_pass = True
        for s in seq_to_screen:
            if s in str(f).upper() or s in str(f.reverse_complement()).upper():
                all_pass = False
        if all_pass == True:
            passed_frags.append(f)
    return(passed_frags)


# In[8]:


def filterToCircularDSeqs(combinedFrags):
    circularPlasmids = []
    i=0
    for f in combinedFrags:
        if f.circular:
            circularPlasmids.append(f)
            i+=1
    return(circularPlasmids)


# In[9]:


def dereplicate_circular_sequences(plasmids):
    dereplicated_plasmids = []
    for i in plasmids:
        ##This is how to do the "if not in" conditional
        itemIn = False
        for j in dereplicated_plasmids:
            ##If lengths don't match, can't be equal
            if len(i) != len(j):
                break
                
            ##Simple comparison
            if i == j or i == j.reverse_complement():
                itemIn = True
                
            ##Exhaustive shifting comparison
            for k in range(1,len(i)):
                shiftedPlasmid = i.shifted(k)
                if shiftedPlasmid == j or shiftedPlasmid == j.reverse_complement():
                    itemIn = True
                    
        if itemIn == False:
            dereplicated_plasmids.append(i)
    return(dereplicated_plasmids)


# In[14]:


if __name__ == '__main__':
    ##input="(003),(010),pJKW 2370,(052),(072),(095)"
    ##input="(002),(030),pJKW 2263,(051),(072),(095)"
    ##input="(002),(010),pJKW 2262,(053),(067),(095)"
    ##input="(002),(010),pJKW 2033,(053),(067),(095)"
    ##input="(002),(010),pJKW 2262,(053),(072),(095)"
    input="pYTK002,(010),pJKW 2033,(053),(072),(095)"

    splitInput=input.split(",")

    plasmidNames = reducedPlasmidNamesToCanonical(splitInput)
    plasmidSeqs,featuresToApply = extractPlasmidSeqsAndFeatures(plasmidNames)
    fragments = digestBsaI(plasmidSeqs)

    initial_frags = [item for sublist in fragments.values() for item in sublist] ##Get rid of the dictionary structure
    initial_frags_filtered = subtract_frags_by_seq(initial_frags) ##Subtract away GFP and chloramphenicol resistance

    pastFrags = []
    keepLooping = True
    i=0
    while keepLooping == True:
        i+=1
        print("Ligation loop:",i)
        combinedFrags = ligate_fragments(pastFrags+initial_frags_filtered)
        print("Number of resulting fragments:",len(combinedFrags))

        for f in combinedFrags:
            if f.circular:
                keepLooping = False
                print("Found circular plasmid. Exiting ligation loop.")
                break
        pastFrags = combinedFrags

    print("Number of ligation products:",len(combinedFrags))
    circularPlasmids = filterToCircularDSeqs(combinedFrags)
    print("Number of ligation products after fitering to circular only:",len(circularPlasmids))
    dereplicatedPlasmids = dereplicate_circular_sequences(circularPlasmids)
    print("Number of ligation products after dereplicating:",len(dereplicatedPlasmids))
    assert len(dereplicatedPlasmids) == 1


# In[16]:


if __name__ == '__main__':
    print("Reapplying feature annotations to sequence and exporting...")
    from Bio.Alphabet.IUPAC     import IUPACAmbiguousDNA
    import datetime
    theSeq = Bio.Seq.Seq(str(dereplicatedPlasmids[0]),IUPACAmbiguousDNA()).upper()
    theRecord = Bio.SeqRecord.SeqRecord(theSeq,id="testID", name="testName",description="testDescription")
    theRecord.annotations["topology"] = "circular"
    theRecord.annotations["molecule_type"] = "ds-DNA"
    today = datetime.datetime.today()
    today_string = today.strftime("%d-%b-%Y").upper() ##%d-%b-%Y -> uppercase
    theRecord.annotations["date"] = today_string

    theRecord = reapply_features_to_record(theRecord,featuresToApply)

    Bio.SeqIO.write(theRecord, "example.gbk", "gb")
    print("finished.")


# In[ ]:




