import Bio.Phylo.PAML.codeml
import Bio.AlignIO
import pandas as pd
import os
from mergedeep import merge
from IPython.display import HTML
from pypdb import *
from Bio.PDB import MMCIFParser
from Bio.PDB.mmtf import MMTFIO
from Bio.PDB.mmtf import MMTFParser
from Bio.PDB.mmtf import DefaultParser
from Bio.PDB.Atom import Atom
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
from pypdb import *
import time
import itertools
import random
import numpy as np


class data:
    '''
    Given list of dndsvalues and list of TAED protein family ids 
    returns a dictionary with the following structure
    "TAED id":
        "Species name" : string
        "Protein name" : string
        "Selection information":
            "N" : float
            "S" : float
            "omega" : float
            "dN" : float
            "dS" : float
            "N*dN" : float
            "S*dS" : float
            "Protein length" : float
        "Alignment sequence" : string
        "Substituted Sites":
            "Positions" : list with strings
            "Ancestral Nucleotide State" : list with strings
            "Ancestral Amino Acid State" : list with strings
            "Current States" : list with strings
        "PDB id" : string
    
    need to save the following TAED files in "files/{proteinfamilyid}/{filename}" format:
        newick tree files
        multiple sequence alignment files 
        paml files
        ancestral sequence files (codeml)
    '''
    def __init__(self, proteinfamilyidlist):
        self.posdict={}
        self.negdict={}
        self.proteinfamilyidlist=proteinfamilyidlist
        self.pamllist=[]
        self.ancestrallist=[]
        self.treelist=[]
        self.alignmentlist=[]

    def getpamlfilenamelist(self):
        '''
        takes list of protein family ids and returns a list with text file names with selection information
        '''
        for id in self.proteinfamilyidlist:
            for i in range (1,10):
                try:
                    filename="files/"+id+"/"+id+f".subTree_tree_{i}.paml_rooted.txt"
                    f=open(filename,"r")
                    self.pamllist.append(filename)
                    f.close()
                except:
                    pass

    def getancestralfilenameslist(self):
        '''
        takes list of protein ids and returns a list with text file names with substituion position information
        '''
        for id in self.proteinfamilyidlist:
            for i in range (1,10):
                try:
                    filename="files/"+id+"/"+id+f"_{i}.RST"
                    f=open(filename,"r")
                    self.ancestrallist.append(filename)
                    f.close()
                except:
                    pass
                


    def gettreefilenamelist(self):
        '''
        takes list of protein family ids and returns a list with text file names with newick trees
        '''
        for id in self.proteinfamilyidlist:
            self.treelist.append("files/"+id+"/"+id+".nhx.txt")


    def getalignmentfilenamelist(self):
        '''
        takes list of protein family ids and returns a list with text file names with newick trees
        '''
        for id in self.proteinfamilyidlist:
            self.alignmentlist.append("files/"+id+"/"+id+".interleaved.txt")


    def getdicts(self):
        print("Creating dictionaries\n")
        '''
        given a list with newick file names, 
        returns a dictionary with species names and proteins names of positively selected and negatively selected proteins
        '''
        templist=[]
        for filename in self.treelist:
                data=open(filename, 'r')
                for line in data:
                    line = line.replace("(", "").replace(")", ",").replace(";","").split(",")
                    for item in line:
                        item=item.split("#")
                        templist.append(item)
                        for i in item:
                            m=i.split(":")
                            for n in m:
                                n=n.strip("]")
                                if "DNDS" in n:
                                    if float(n[5:])>1:
                                        posindex=templist.index(item)
                                        if ("BSV" not in templist[posindex][0]) and ("root" not in templist[posindex][0]):
                                            self.posdict[templist[posindex][0]]={"Species name":templist[posindex][1], "Protein name":m[0], "Omega":float(n[5:])}
                                    elif float(n[5:])<0.5:
                                        negindex=templist.index(item)
                                        if ("BSV" not in templist[negindex][0]) and ("root" not in templist[negindex][0]):
                                            self.negdict[templist[negindex][0]]={"Species name":templist[negindex][1], "Protein name":m[0], "Omega":float(n[5:])}
        return self.posdict, self.negdict

    def getseq(self):
        print("Getting sequences\n")
        '''
        given a dictionary with TAED ids as keys, laods protein alignment sequence related to the TAED ids
        TAED id:
            "Alignment sequence" : string
        '''
        for filename in self.alignmentlist:
            alignments=Bio.AlignIO.read(open(filename),"phylip")
            for record in alignments:
                for id in self.posdict.keys():
                    if record.id==id:
                        self.posdict[id]["Alignment sequence"]=record.seq
                for id1 in self.negdict.keys():
                    if record.id==id1:
                        self.negdict[id1]["Alignment sequence"]=record.seq
        
        return self.posdict, self.negdict

    def getdndsinfo(self):
        print("Getting selection information\n")
        '''
        given a dictionary with TAED ids as keys, loads information about selection on the protein related to TAED id
        "TAED id":
            "Selection information":
                "N" : float
                "S" : float
                "omega" : float
                "dN" : float
                "dS" : float
                "N*dN" : float
                "S*dS" : float
                "Protein length" : float
        '''
        posdict={}
        negdict={}
        for filename in self.pamllist:
            data=Bio.Phylo.PAML.codeml.read(filename)
            omegatree=data["NSsites"][0]["omega tree"].split(" ")
            branches=data["NSsites"][0]["parameters"]["branches"]
            for i in omegatree:
                for id in self.posdict.keys():
                    if id in i:
                        dnds=self.posdict[id]["Omega"]
                        for m in branches:
                            for n in branches[m]:
                                if branches[m][n]==dnds and branches[m]["N*dN"]>=2:
                                    temp=branches[m].copy()
                                    temp["Protein length"] = (branches[m]["N"]+branches[m]["S"])/3
                                    posdict[id]={"Selection information": temp}
                for id1 in self.negdict.keys(): 
                    if id1 in i:
                        dnds=self.negdict[id1]["Omega"]
                        for m in branches:
                            for n in branches[m]:
                                if branches[m][n]==dnds and branches[m]["N*dN"]>=2:
                                    temp=branches[m].copy()
                                    temp["Protein length"] = (branches[m]["N"]+branches[m]["S"])/3
                                    negdict[id1]={"Selection information": temp}
        merge(self.posdict, posdict)
        merge(self.negdict, negdict)
        return self.posdict, self.negdict


    def getpossubspositions(self):
        print("Getting the subs positions\n")
        '''
        given a dictionary with TAED ids as keys, loads information about substituted sites on protein related to TAED id
        "TAED id":
            "Substituted Sites":
                "Positions" : list with strings
                "Ancestral Nucleotide State" : list with strings
                "Ancestral Amino Acid State" : list with strings
                "Current States" : list with strings
        '''
        dictcopy=self.posdict.copy()
        for filename in self.ancestrallist:
            for id in self.posdict.keys():
                data=open(filename,"r")
                datalines=data.readlines()
                data.close()
                newdatafile=open(str(filename)[:-4]+"_new.csv",'w')
                for line in datalines:
                    search = "(T" + id + ")  (n="
                    if search in line:
                        index=datalines.index(line)
                        for line1 in datalines[(index+1):]:
                            if "Branch" in line1:
                                index1=datalines.index(line1)
                                temp1=datalines[(index+1):(index1-2)] 
                                for templine in temp1:
                                    templine=templine.lstrip()
                                    newdatafile.write(templine)
                                break
                df=pd.read_csv(str(filename)[:-4]+"_new.csv", sep=" ", header=None, names=["Substituted Positions", "Ancestral Nucleotide State","Ancestral Amino Acid State", "Probability of Change", "Arrows", "Current Nucleotide State", "Current Amino Acid State"])
                if df.empty:
                    #print("empty")
                    os.remove(str(filename)[:-4]+"_new.csv")
                else:
                    df.drop(df.loc[df["Probability of Change"]=="->"].index, inplace=True)
                    df = df.astype({"Probability of Change":'float'})
                    df.drop(df.loc[df["Probability of Change"]<=0.5].index, inplace=True)
                    df.drop(df.loc[df["Current Nucleotide State"]=="---"].index, inplace=True)
                    df.drop("Probability of Change", inplace=True, axis=1)
                    df.drop("Arrows", inplace=True, axis=1)
                    dict1=df.to_dict("list")
                    dictcopy[id]["Substituted Sites"]=dict1
                    os.remove(str(filename)[:-4]+"_new.csv")
        for id1 in self.posdict.keys():
            try:
                dictcopy[id1]["Substituted Sites"]
            except:
                del dictcopy[id1]
        self.posdict=dictcopy
        return self.posdict


    def getnegsubspositions(self):
        print("Getting the subs positions\n")
        '''
        given a dictionary with TAED ids as keys, loads information about substituted sites on protein related to TAED id
        "TAED id":
            "Substituted Sites":
                "Positions" : list with strings
                "Ancestral Nucleotide State" : list with strings
                "Ancestral Amino Acid State" : list with strings
                "Current Nucleotide State" : list with strings
                "Current Amino Acid State" : list with strings
        '''
        dictcopy=self.negdict.copy()
        for filename in self.ancestrallist:
            for id in self.negdict.keys():
                data=open(filename,"r")
                datalines=data.readlines()
                data.close()
                newdatafile=open(str(filename)[:-4]+"_new.csv",'w')
                for line in datalines:
                    search = "(T" + id + ")  (n="
                    if search in line:
                        index=datalines.index(line)
                        for line1 in datalines[(index+1):]:
                            if "Branch" in line1:
                                index1=datalines.index(line1)
                                temp1=datalines[(index+1):(index1-2)] 
                                for templine in temp1:
                                    templine=templine.lstrip()
                                    newdatafile.write(templine)
                                break
                df=pd.read_csv(str(filename)[:-4]+"_new.csv", sep=" ", header=None, names=["Substituted Positions", "Ancestral Nucleotide State","Ancestral Amino Acid State", "Probability of Change", "Arrows", "Current Nucleotide State", "Current Amino Acid State"])
                if df.empty:
                    os.remove(str(filename)[:-4]+"_new.csv")
                else:
                    df.drop(df.loc[df["Probability of Change"]=="->"].index, inplace=True)
                    df = df.astype({"Probability of Change":'float'})
                    df.drop(df.loc[df["Probability of Change"]<=0.5].index, inplace=True)
                    df.drop(df.loc[df["Current Nucleotide State"]=="---"].index, inplace=True)
                    df.drop("Probability of Change", inplace=True, axis=1)
                    df.drop("Arrows", inplace=True, axis=1)
                    dict1=df.to_dict("list")
                    dictcopy[id]["Substituted Sites"]=dict1
                    os.remove(str(filename)[:-4]+"_new.csv")
        for id1 in self.negdict.keys():
            try:
                dictcopy[id1]["Substituted Sites"]
            except:
                del dictcopy[id1]
        self.negdict=dictcopy
        return self.negdict

    def searchpdbid(self):
        print("Matching PDB ids\n")
        posdict=self.posdict.copy()
        negdict=self.negdict.copy()
        '''
        given a dictionary with TAED ids as keys,
        loads PDB ids related to TAED ids to the dictionary
        loads xyz coordinates of each alpha carbon atom in the protein
        '''
        for id in self.posdict.keys():
            seq=''
            try:
                for aminoacid in self.posdict[id]["Alignment sequence"]:
                    if not aminoacid=="-":
                        seq+=aminoacid
                q = Query(seq, 
                query_type="sequence", 
                return_type="polymer_entity")
                posdict[id]["PDB id"]=q.search()["result_set"][0]["identifier"]
                posdict[id]["PDB aligned seq"]=q.search()["result_set"][0]["services"][0]["nodes"][0]["match_context"][0]["query_aligned_seq"]
                structure = MMTFParser.get_structure_from_url(posdict[id]["PDB id"][:4])
                coordinates={}
                for atom in structure.get_atoms():
                    if "Atom CA" in str(atom):
                        respos=atom.get_parent().get_id()[1]
                        coordinates[respos]=atom.coord.tolist()
                posdict[id]["CA coordinates"]=coordinates
            except:
                del posdict[id]
                pass
        for id1 in self.negdict.keys():
            seq1=''
            try:
                for aminoacid1 in self.negdict[id1]["Alignment sequence"]:
                    if not aminoacid1=="-":
                        seq1+=aminoacid1
                q1 = Query(seq1, 
                query_type="sequence", 
                return_type="polymer_entity")
                negdict[id1]["PDB id"]=q1.search()["result_set"][0]["identifier"]
                negdict[id1]["PDB aligned seq"]=q.search()["result_set"][0]["services"][0]["nodes"][0]["match_context"][0]["query_aligned_seq"]
                structure = MMTFParser.get_structure_from_url(negdict[id1]["PDB id"][:4])
                coordinates={}
                for atom in structure.get_atoms():
                    if "Atom CA" in str(atom):
                        respos=atom.get_parent().get_id()[1]
                        coordinates[respos]=atom.coord.tolist()
                negdict[id1]["CA coordinates"]=coordinates
            except:
                del negdict[id1]
                pass
        self.posdict=posdict
        self.negdict=negdict
        return self.posdict, self.negdict

start = time.time()
proteinfamilyid=["18975", '9781', '6283', '19419', '6543', '7131']
templist=["18975"]
data1=data(proteinfamilyid)
data.gettreefilenamelist(data1)
posdict,negdict=data.getdicts(data1)
data.getalignmentfilenamelist(data1)
posdict1,negdict1=data.getseq(data1)
data.getpamlfilenamelist(data1)
posdict2,negdict2=data.getdndsinfo(data1)
data.getancestralfilenameslist(data1)
posdict3=data.getpossubspositions(data1)
negdict3=data.getnegsubspositions(data1)
posdict4, negdict4=data.searchpdbid(data1)

#writing the dictionary information to file
with open("negatively_selected_proteins.txt", 'w') as f: 
    f.write("Information about negatively selected proteins from TAED and PDB databases.\n")
    f.write("Number of proteins:\n")
    f.write(str(len(negdict4)))
    f.write("\n\n")
    for key, value in negdict4.items():
        f.write('%s:\n' % (key))
        for key1, value1 in negdict4[key].items():
            f.write('\t%s:%s\n' % (key1, value1))
        f.write("\n\n")
    f.close()

#writing the dictionary information to file
with open("positively_selected_proteins.txt", 'w') as f: 
    f.write("Information about positively selected proteins from TAED database.\n\n")
    f.write("Number of proteins:\n")
    f.write(str(len(posdict4)))
    f.write("\n\n")
    for key, value in posdict4.items(): 
        f.write('%s:\n' % (key))
        for key1, value1 in posdict4[key].items():
            f.write('\t%s:%s\n' % (key1, value1))
        f.write("\n\n")
    f.close()