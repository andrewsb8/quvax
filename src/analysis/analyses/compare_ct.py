import pandas as pd

def ct_to_dataframe(ct_file):
    """
    Takes base connectivity data in the form of .ct files and converts it to a data frame.
    
    """
    df = pd.read_csv(ct_file, delim_whitespace=True, skiprows=1, header=None)
    df.columns = ["Index", "Nucleotide", "Previous", "Next", "Paired With", "Counter"]
    return df

#accessing the files from my desktop, will have to be changed when it enters the codebase
ct_file_path_target = '/Users/jakeabraham/Desktop/5s_Acetobacter-sp.-1_archiveII.ct'
ct_file_path_quvax = '/Users/jakeabraham/Desktop/5s_Acetobacter-sp.-1_quvax.ct'

targetct = ct_to_dataframe(ct_file_path_target)
targetpairings = targetct["Paired With"]

quvaxct = ct_to_dataframe(ct_file_path_quvax)
quvaxpairings = quvaxct["Paired With"]

def compare_tables(test_pairings, ref_pairings):
    """
    Inputs the data frame, outputs the sensitivity, PPV, F1 score, 
    specificity, bases correctly paired, and bases with incorrect pairing.

    """
    print("Bases correctly paired: ", truepos)
    print("Bases with misidentified pairings (either missed or incorrect pairing): ", falsepos + falseneg)
    print("Sensitivity: ", sensitivity(truepos, falseneg), " ; PPV: ", PPV(truepos, falsepos)," ; F1: ", F1(truepos, falsepos, falseneg), " ; Specificity: ", specificity(trueneg, falsepos))

def sensitivity(truepos, falseneg):
    return truepos/(truepos+falseneg)

def PPV(truepos, falsepos):
    return truepos/(truepos+falsepos)
    

def F1(truepos, falsepos, falseneg):
    return 2*truepos/(2*truepos + falsepos + falseneg)

def specificity(trueneg, falsepos):
    return trueneg/(trueneg+falsepos)
    
def truthvalues(test_pairings, ref_pairings):
    """
    Gathers data from the tables to find how many true positives, true negatives, 
    false positives, and false negatives there are given both data sets . 
    
    """
    numtruepos = 0
    numtrueneg = 0
    numfalsepos = 0
    numfalseneg = 0
    for i in range(len(test_pairings)):
        if test_pairings[i] != 0 and test_pairings[i] == ref_pairings[i]:
            numtruepos += 1
        if test_pairings[i] == 0 and ref_pairings[i] == 0:
            numtrueneg += 1
        if test_pairings[i] != 0 and test_pairings[i] != ref_pairings[i]:
                numfalsepos += 1
        if test_pairings[i] == 0 and ref_pairings[i] != 0:
            numfalseneg += 1
    return numtruepos, numtrueneg, numfalsepos, numfalseneg

truepos, trueneg, falsepos, falseneg = truthvalues(quvaxpairings,targetpairings)
compare_tables(quvaxpairings,targetpairings) 
