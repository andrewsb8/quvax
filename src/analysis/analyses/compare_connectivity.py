#Script to compare connectivity tables from output of design.py with known connectivity tables from a database
import pandas as pd

#Step 0 is to create a function that can input a ct table and converts it into a dataframe. 

def ct_to_dataframe(ct_file):
    df = pd.read_csv(ct_file, delim_whitespace=True, skiprows=1, header=None)
    df.columns = ["Index", "Nucleotide", "Previous", "Next", "Paired With", "Counter"]
    return df

#accessing the files from my desktop, will have to be changed when it enters the codebase
ct_file_path_target = '/Users/jakeabraham/Desktop/5s_Acetobacter-sp.-1_archiveII.ct'
ct_file_path_quvax = '/Users/jakeabraham/Desktop/5s_Acetobacter-sp.-1_quvax.ct'


#Step 1 is to load in the target connectivities, as well as the ones from the Quvax codebase. Indeces is just the column that shows which nucleotide index that one is paired with.

targetct = ct_to_dataframe(ct_file_path_target)
targetpairings = targetct["Paired With"]


quvaxct = ct_to_dataframe(ct_file_path_quvax)
quvaxpairings = quvaxct["Paired With"]


#Step 2 is to write an overall function called "compare_tables" that outputs sensitivity, PPV, F1 score, and the specificity


def compare_tables(pairings1, pairings2):
    print("Bases correctly paired: ", truepos(pairings1, pairings2))
    print("Bases with misidentified pairings (either missed or incorrect pairing): ", falsepos(pairings1, pairings2)+falseneg(pairings1, pairings2))
    print("Sensitivity: ", sensitivity(pairings1, pairings2), " ; PPV: ", PPV(pairings1, pairings2)," ; F1: ", F1(pairings1, pairings2), " ; Specificity: ", specificity(pairings1, pairings2))
    

#Step 3 is to create functions that calculate the sensitivity, PPV, F1 score, and the specificity given the inputs 

def sensitivity(pairings1, pairings2):
    sens = truepos(pairings1, pairings2)/(truepos(pairings1, pairings2)+falseneg(pairings1, pairings2))
    return sens

def PPV(pairings1, pairings2):
    ppv = truepos(pairings1, pairings2)/(truepos(pairings1, pairings2)+falsepos(pairings1, pairings2))
    return ppv

def F1(pairings1, pairings2):
    f = 2*sensitivity(pairings1, pairings2)*PPV(pairings1, pairings2)/(sensitivity(pairings1, pairings2)+PPV(pairings1, pairings2))
    return f

def specificity(pairings1, pairings2):
    spec = trueneg(pairings1, pairings2)/(trueneg(pairings1, pairings2)+falsepos(pairings1, pairings2))
    return spec

#Step 4 is to zoom even further to create functions that gather data directly from the tables. 
#Want to know how many true positives, true negatives, false positives, and false negatives there are. 

def truepos(pairings1, pairings2):
    numtruepos = 0
    for i in range(len(pairings1)):
        if pairings1[i] == 0:
            numtruepos += 0
        if pairings1[i] != 0:
            if pairings1[i] == pairings2[i]:
                numtruepos += 1
            else:
                numtruepos += 0
    return numtruepos

def trueneg(pairings1, pairings2):
    numtrueneg = 0
    for i in range(len(pairings1)):
        if pairings1[i] == 0:
            if pairings2[i] == 0:
                numtrueneg += 1
            else:
                numtrueneg += 0
        if pairings1[i] != 0:
            numtrueneg += 0
    return numtrueneg

def falsepos(pairings1, pairings2):
    numfalsepos = 0
    for i in range(len(pairings1)):
        if pairings1[i] == 0:
            numfalsepos += 0
        if pairings1[i] != 0:
            if pairings1[i] != pairings2[i]:
                numfalsepos += 1
            else:
                numfalsepos += 0
    return numfalsepos

def falseneg(pairings1, pairings2):
    numfalseneg = 0
    for i in range(len(pairings1)):
        if pairings1[i] == 0:
            if pairings2[i] != 0:
                numfalseneg += 1
            else:
                numfalseneg += 0
        if pairings1[i] != 0:
            numfalseneg += 0
    return numfalseneg

compare_tables(quvaxpairings,targetpairings) #first entry should be our data, and the second should be the database
