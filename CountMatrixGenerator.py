import sys
import os
import glob
import gffutils
from collections import defaultdict
import csv

###variables:
# for filename attributes
filebeginning="merge_graphs_"
fileending="gff3"
charac="genenames"
csv_fill_attr="_for_"

#for all possible types of alternative splicing events
splicing=["alt_3prime","alt_5prime","exon_skip","intron_retention","mult_exon_skip"]


### paths
# path to gff3 files (named S250Nr)
gff3_filepath="C:/Users/Rosa/Documents/Uni/7.Semester/Bachelorarbeit/alternative_splicing_analysis/data/gff_data/01-fastq/"

#csv result path for event count tables
csv_path="C:/Users/Rosa/Documents/Uni/7.Semester/Bachelorarbeit/alternative_splicing_analysis/countTables/"

test_result="C:/Users/Rosa/Documents/Uni/7.Semester/Bachelorarbeit/"

metasheet_path="C:/Users/Rosa/Documents/Uni/7.Semester/Bachelorarbeit/alternative_splicing_analysis/MetaSheet_RNA_Seq_Samples.csv"
test="C:/Users/Rosa/Documents/Uni/7.Semester/Bachelorarbeit/alternative_splicing_analysis/TEST.csv"

# return a counttable for one splicing type, that contains the eventcounts of every patient for a correponding gene.
# the columns represent the patients
# the rows repesent the gene
# the entries of the matrix are the eventcounts for a specific gene, calculated for every patient


def readInMetaSheet(metasheet_path):
    csv_dict=defaultdict(list)
    #keys=getKeys(metasheet_path)
    try:                                                                    
        with open(metasheet_path, 'r') as csvfile:                                # write the event count matrix to the csv file
            reader = csv.DictReader(csvfile)
            for rows in reader:
                for key in rows.keys():
                    csv_dict[key].append(rows[key])

    except IOError:
        print("I/O error") 
    checkMetasheet(csv_dict)
    return(csv_dict)
   

def checkMetasheet(data):
    keys=list(data.keys())
    if(len(keys)>=3):
        for key in keys:
            if(None in data[key]):
                raise ValueError("Metasheet is wrong: the number of entries in the columns differ")
        for val in data[keys[2]]:
            if(val!="normal" and val!="tumour"):
                raise ValueError("Metasheet is wrong: group is neither normal nor tumour" )
    else:
        raise ValueError("Metasheet is wrong: columns of the metasheet is smaller than 2")
    

def getCountTables(infilepath,outfilepath,splicetype,metasheet):
    """Function for event count calculation of a gene.

    `Calculation of event counts for one type of splicing
     correponding to a specific gene and a sample of patients from two cohorts (cancer-related and healthy)
     The function returns a count matrix for one splicing type, that contains the eventcounts
     of every patient for a specific gene.
     Explanation of the matrix:
     rows    : repesent the gene
     columns : represent the patients
     entries : represent the event count (of a patient and a gene correpsonding to the row)´

    Args:
        infilepath  (str): The first parameter is the path to the patient folders (named S250NrX),
                           that containing the gff3 files of every splice type in a folder called spladdrout.
        outfilepath (str): The second parameter is the path to the directory, where the count matrix should be
                           saved in csv format.
        splicetype  (str): The third parameter is the type of the alternative splice event.

    Returns:
        defaultdict dict: The return matrix (named data), that contains all event counts for all genes of a splice type found in the gff3 file.

    .. spladder documentation:
        https://github.com/ratschlab/spladder/wiki

    .. ggf3 file explanation:
        https://github.com/ratschlab/spladder/wiki/File-Format-Descriptions


    """
    data=defaultdict(dict)                                          # return matrix
    patients=[]                                                     # list of all patients
    events=defaultdict(list)  
    metaData=readInMetaSheet(metasheet)     
    metaKeys=list(metaData.keys() )                     # defaultdict for all events , containing a list of counts per gene              
    patientfolders=metaData[metaKeys[0]]                         # access to the patient folders
    print("Patientfolders found for %s: %s"%(splicetype,patientfolders))
    for patient in patientfolders:
        print("Calculation of eventcounts for patient: %s" %(patient))
        if(not(patient in patients)):                         
            patients.append(patient)                                # add new patients to patientlist      
        sourcebase="{}{}/spladdrout/".format(infilepath,patient)      # access to the gff3 files
        os.chdir(sourcebase)
        for file in glob.glob("%s%s*" %(filebeginning,splicetype)): # access to the gff3 file of a certain splicing type
            db=gffutils.create_db(file,':memory:')                  
            geneEvents=dict()                                       # init dictionary for the event counts of the current gff3 file
            for gene in db.features_of_type('gene'):
                if(gene.source!=splicetype):                        # check if the splicing type is contained in the file of interest
                    print("Warning: wrong splicing event found")
                genename=gene.attributes['GeneName'][0]             # get the genename of the current event

                if(genename in geneEvents.keys()):                  # sum up all events of one patient found for the same gene
                    geneEvents[genename]+=1
                else:
                    geneEvents[genename]=1                          # otherwise initialize a new count for a gene which was not seen before
            for gene in geneEvents.keys():
                if(gene in events.keys()):                          # if the gene is in the event matrix already
                    events[gene].append(geneEvents[gene])           # extend the matrix with the count of the current patient for this gene
                else:                                               
                    pat=0                                           # else add a new row to the matrix for the "new" gene 
                    while(pat<len(patients)-1):                     # and fill all counts of the passed patients with zero
                        events[gene].append(0)                      
                        pat+=1
                    events[gene].append(geneEvents[gene])           # last, append the calculated counts for the current patient

            for gene in events.keys():
                if(not(gene in geneEvents.keys())):                 # if a gene is missed in the current patient, add a zero for the patientscolumn
                    events[gene].append(0)

    dict_data=list()                                                # initialize the dictionary for the csv writer
    csv_columns=[charac]+patients                                  # initialize the csv columns
    
    for gene in events.keys():                                      # create the write structure for the csv writer
        dict_gene=dict()
        dict_gene[charac]=gene
       
        num=0
        for pat in patients:
            dict_gene[pat]=events[gene][num]        
            num+=1
        dict_data.append(dict_gene)

    csv_file = "%scountTables%s%s.csv"%(outfilepath,csv_fill_attr,splicetype)  # name the csv file corresponding to the current splicing type
    try:                                                                    
        with open(csv_file, 'w') as csvfile:                                # write the event count matrix to the csv file
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in dict_data:
                writer.writerow(data)
    except IOError:
        print("I/O error") 
                    
    return data                                                            # return the count matrix

def getCountTables_mod(infilepath,outfilepath,splicetype,metasheet):
    """Function for event count calculation of a gene.

    `Calculation of event counts for one type of splicing
     correponding to a specific gene and a sample of patients from two cohorts (cancer-related and healthy)
     The function returns a count matrix for one splicing type, that contains the eventcounts
     of every patient for a specific gene.
     Explanation of the matrix:
     rows    : repesent the gene
     columns : represent the patients
     entries : represent the event count (of a patient and a gene correpsonding to the row)´

    Args:
        infilepath  (str): The first parameter is the path to the patient folders (named S250NrX),
                           that containing the gff3 files of every splice type in a folder called spladdrout.
        outfilepath (str): The second parameter is the path to the directory, where the count matrix should be
                           saved in csv format.
        splicetype  (str): The third parameter is the type of the alternative splice event.

    Returns:
        defaultdict dict: The return matrix (named data), that contains all event counts for all genes of a splice type found in the gff3 file.

    .. spladder documentation:
        https://github.com/ratschlab/spladder/wiki

    .. ggf3 file explanation:
        https://github.com/ratschlab/spladder/wiki/File-Format-Descriptions


    """
    data=defaultdict(dict)                                          # return matrix
    patients=[]                                                     # list of all patients
    events=defaultdict(list)  
    metaData=readInMetaSheet(metasheet)     
    metaKeys=list(metaData.keys() )                     # defaultdict for all events , containing a list of counts per gene              
    patientfolders=metaData[metaKeys[0]]                         # access to the patient folders
    patientIDs=metaData[metaKeys[1]]
    print("Patientfolders found for %s: %s"%(splicetype,patientfolders))
    print("Patients found for %s: %s"%(splicetype,patientIDs))
   # print(metaData)
    for i in range(0,len(patientfolders)):
       # print(i)
        patientfolder=patientfolders[i]
        print("Patientfolder: %s" %(patientfolder))
        patient=patientIDs[i]
        print("Calculation of eventcounts for patient: %s" %(patient))
        if(not(patient in patients)):                         
            patients.append(patient)                                # add new patients to patientlist      
        sourcebase="{}{}/spladdrout/".format(infilepath,patientfolder)      # access to the gff3 files
        os.chdir(sourcebase)
        for file in glob.glob("%s%s*" %(filebeginning,splicetype)): # access to the gff3 file of a certain splicing type
            db=gffutils.create_db(file,':memory:')                  
            geneEvents=dict()                                       # init dictionary for the event counts of the current gff3 file
            for gene in db.features_of_type('gene'):
                if(gene.source!=splicetype):                        # check if the splicing type is contained in the file of interest
                    print("Warning: wrong splicing event found")
                genename=gene.attributes['GeneName'][0]             # get the genename of the current event

                if(genename in geneEvents.keys()):                  # sum up all events of one patient found for the same gene
                    geneEvents[genename]+=1
                else:
                    geneEvents[genename]=1                          # otherwise initialize a new count for a gene which was not seen before
            for gene in geneEvents.keys():
                if(gene in events.keys()):                          # if the gene is in the event matrix already
                    events[gene].append(geneEvents[gene])           # extend the matrix with the count of the current patient for this gene
                else:                                               
                    pat=0                                           # else add a new row to the matrix for the "new" gene 
                    while(pat<len(patients)-1):                     # and fill all counts of the passed patients with zero
                        events[gene].append(0)                      
                        pat+=1
                    events[gene].append(geneEvents[gene])           # last, append the calculated counts for the current patient

            for gene in events.keys():
                if(not(gene in geneEvents.keys())):                 # if a gene is missed in the current patient, add a zero for the patientscolumn
                    events[gene].append(0)

    dict_data=list()                                                # initialize the dictionary for the csv writer
    csv_columns=[charac]+patients                                  # initialize the csv columns
    
    for gene in events.keys():                                      # create the write structure for the csv writer
        dict_gene=dict()
        dict_gene[charac]=gene
       
        num=0
        for pat in patients:
            dict_gene[pat]=events[gene][num]        
            num+=1
        dict_data.append(dict_gene)

    csv_file = "%scountTables%s%s.csv"%(outfilepath,csv_fill_attr,splicetype)  # name the csv file corresponding to the current splicing type
    try:                                                                    
        with open(csv_file, 'w') as csvfile:                                # write the event count matrix to the csv file
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in dict_data:
                writer.writerow(data)
    except IOError:
        print("I/O error") 
                    
    return data                                                            # return the count matrix


def main(gff3_filepath,csv_path):
    """Main function.

    `Calculation of event counts for all types of alterntive splicing using the function getCountTables´

    Args:
        gff3_filepath (str): The first parameter is the path to the patient folders (named S250NrX),
                             that containing the gff3 files of every splice type in a folder called spladdrout.
        csv_path      (str): The second parameter is the path to the directory, where the count matrix should be
                             saved in csv format.

    Returns:
        nothing

    .. getCounttable:
        getCounttable.__doc__

    """
    for splice in splicing:
        getCountTables_mod(gff3_filepath,csv_path,splice,metasheet_path)


#data=getCountTables_mod(gff3_filepath,test_result,splicing[0],test)
main(gff3_filepath,csv_path)
#print(getCountTables.__doc__)
#print(main.__doc__)
#dic=readInMetaSheet(metasheet_path)
#print(dic)
