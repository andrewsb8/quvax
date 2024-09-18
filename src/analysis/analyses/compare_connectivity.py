import pandas as pd


def ct_to_dataframe(ct_file):
    """
    This function takes base connectivity data in the form of .ct files and converts it to
    a data frame.
    
    """
    df = pd.read_csv(ct_file, delim_whitespace=True, skiprows=1, header=None)
    df.columns = ["Index", "Nucleotide", "Previous", "Next", "Paired With", "Counter"]
    return df

#accessing the files from my desktop, will have to be changed when it enters the codebase
#ct_file_path_target = '/Users/jakeabraham/Desktop/5s_Acetobacter-sp.-1_archiveII.ct'
#ct_file_path_quvax = '/Users/jakeabraham/Desktop/5s_Acetobacter-sp.-1_quvax.ct'

targetct = ct_to_dataframe(ct_file_path_target)
targetpairings = targetct["Paired With"]

quvaxct = ct_to_dataframe(ct_file_path_quvax)
quvaxpairings = quvaxct["Paired With"]


def compare_tables(test_pairings, ref_pairings):
    """
    Function that inputs the base pairing data from both our code and from a reference database, 
    and outputs the sensitivity, PPV, F1 score, specificity, bases correctly paired, and .

    """
#    truepos, trueneg, falsepos, falseneg = truthvalues(test_pairings,ref_pairings)
    print("Bases correctly paired: ", truepos)
    print("Bases with misidentified pairings (either missed or incorrect pairing): ", falsepos + falseneg)
    print("Sensitivity: ", sensitivity(test_pairings, ref_pairings), " ; PPV: ", PPV(test_pairings, ref_pairings)," ; F1: ", F1(test_pairings, ref_pairings), " ; Specificity: ", specificity(test_pairings, ref_pairings))

def sensitivity(test_pairings, ref_pairings):
    return truepos/(truepos+falseneg)

def PPV(test_pairings, ref_pairings):
    return truepos/(truepos+falsepos)
    
def F1(test_pairings, ref_pairings):
    return 2*truepos/(2*truepos + falsepos + falseneg)

def specificity(test_pairings, ref_pairings):
    return trueneg/(trueneg+falsepos)
    
  
def truthvalues(test_pairings, ref_pairings):
    """
    This function gathers data from the tables to find how many true positives, true negatives, 
    false positives, and false negatives there are given both data sets . 
    
    """
    numtruepos = 0
    numtrueneg = 0
    numfalsepos = 0
    numfalseneg = 0
    for i in range(len(test_pairings)):
        if test_pairings[i] != 0 and test_pairings[i] == ref_pairings[i]:
            numtruepos += 1
    for i in range(len(test_pairings)):
        if test_pairings[i] == 0 and ref_pairings[i] == 0:
            numtrueneg += 1
    for i in range(len(test_pairings)):
        if test_pairings[i] != 0 and test_pairings[i] != ref_pairings[i]:
            numfalsepos += 1
    for i in range(len(test_pairings)):
        if test_pairings[i] == 0 and ref_pairings[i] != 0:
            numfalseneg += 1
    return numtruepos, numtrueneg, numfalsepos, numfalseneg

truepos, trueneg, falsepos, falseneg = truthvalues(quvaxpairings,targetpairings)
compare_tables(quvaxpairings,targetpairings) #first entry should be our data, and the second should be the database
