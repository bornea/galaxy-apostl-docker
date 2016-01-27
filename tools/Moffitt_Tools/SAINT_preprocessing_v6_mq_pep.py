#See other note about standardizing header and adding all authors. We can say "original script by Brent Kuenzi" to give Brent with coming up with the original version of this file.
#######################################################################################
# Python-code: SAINT pre-processing from maxquant "Samples Report" output
# Author: Brent Kuenzi
#######################################################################################
# This program reads in a raw maxquant "Samples Report" output and a user generated
# bait file and autoformats it into prey and interaction files for SAINTexpress 
# analysis
#######################################################################################
import sys
import urllib2
import os
#######################################################################################
## REQUIRED INPUT ##

# 1) infile: maxquant "Samples Report" output
# 2) baitfile: SAINT formatted bait file generated in Galaxy
# 3) fasta_db: fasta database for use (defaults to SwissProt_HUMAN_2014_08.fasta)
# 4) prey: Y or N for generating a prey file 
# 5) make_bait: String of bait names, assignment, and test or control boolean 
#######################################################################################
mq_file = sys.argv[1]
ins_path = "/galaxy-apostl-docker/tools/Moffitt_Tools/"
names_path = str(ins_path) + r"uniprot_names.txt"
cmd = r"Rscript "+ str(ins_path) +"pre_process_protein_name_set.R " + str(mq_file) + " " + str(names_path)
os.system(cmd)

infile = "./tukeys_output.txt" #maxquant "Samples Report" output
prey = sys.argv[2] # Y or N
fasta_db = sys.argv[3]
if fasta_db == "None":
    fasta_db = str(ins_path)  + "SwissProt_HUMAN_2014_08.fasta"
make_bait= sys.argv[6]
bait_bool = sys.argv[9]

def bait_create(baits, infile):
    #Takes the Bait specified by the user and makes them into a Bait file and includes a check to make sure they are using valid baits.
    baits = make_bait.split()
    i = 0
    bait_file_tmp = open("bait.txt", "wr")
    order = [] 
    bait_cache = []
    while i < len(baits):
        if baits[i+2] == "true":
            T_C = "C"
        else:
            T_C = "T"
        line1 = baits[i] + "\t" + baits[i+1] + "\t" + T_C + "\n"
        q = open(infile,"r")
        for line2 in q:
           line2 = line2.replace("\"", "")
           line2 = line2.replace(r"Intensity.", "") #R coerces "-" into "." this changes them back and remove Intensity from the Bait names.
           line2 = line2.replace(r".", r"-")
           temp = line2.split()
           if "mapped_protein" in str(temp):
                #If the bait is in the original file then write to cache it if not exit.
                if baits[i] in temp:
                    number_bait = temp.index(str(baits[i]))
                    number_bait = number_bait - 9 #What's this?
                    bait_cache.append((number_bait, str(line1)))
                else:
                    print "Error: bad bait " + str(baits[i]) #Standardize messages to the user (see other comment about stdout vs stderr). Can this be more helpful like suggesting what the problem is based on what the script can't find?
                    # AB: Sure we can do lots of stuff, this is just telling them what they typed as a bait is not found in the file as expected.
                    sys.exit()
           else: 
                pass                    
        i = i + 3
    #Writes cache to file.
    bait_cache.sort()
    for line in bait_cache:
        bait_file_tmp.write(line[1])            
        
    bait_file_tmp.close()  


if bait_bool == 'false':
    bait_create(make_bait, infile)
    baitfile = "bait.txt" 
else:
    bait_temp_file = open(sys.argv[10], 'r')
    bait_cache = bait_temp_file.readlines()
    bait_file_tmp = open("bait.txt", "wr")
    for line in bait_cache:
        bait_file_tmp.write(line)                    
    bait_file_tmp.close()
    baitfile = "bait.txt" 


class ReturnValue1(object):
    def __init__(self, sequence, gene):
     self.seqlength = sequence
     self.genename = gene
class ReturnValue2(object):
    def __init__(self, getdata, getproteins, getheader):
        self.data = getdata
        self.proteins = getproteins
        self.header = getheader

def main(maxquant_input, make_bait):  
    #bait_check(baitfile, maxquant_input)
    make_inter(maxquant_input)
    if prey == 'true':
        make_prey(maxquant_input)
        no_error_inter(maxquant_input)
        os.rename('prey.txt', sys.argv[5])
    elif prey == 'false':
        if os.path.isfile('error proteins.txt') == True:
            no_error_inter(maxquant_input)
        pass
    elif prey != 'true' or 'false':
        sys.exit("Invalid Prey Argument: Y or N")
    os.rename('inter.txt', sys.argv[4])
    os.rename("bait.txt", sys.argv[7])


def get_info(uniprot_accession_in): # get aa lengths and gene name
    error = open('error proteins.txt', 'a+')
    data = open(fasta_db,'r')
    lines = data.readlines()
    db_len = len(lines)
    seqlength = 0
    count = 0
    for i in lines:           
        if ">sp" in i:
            if uniprot_accession_in == i.split("|")[1]:
                match = count+1
                if 'GN=' in i:
                    lst = i.split('GN=')
                    lst2 = lst[1].split(' ')
                    genename = lst2[0]
                if 'GN=' not in i:
                    genename = 'NA'
                while ">sp" not in lines[match]:
                    if match <= db_len:
                        seqlength = seqlength + len(lines[match].strip())
                        match = match + 1
                    else:
                        break
                return ReturnValue1(seqlength, genename)
        count = count + 1
        

    if seqlength == 0:
        error.write(uniprot_accession_in + '\t' + "Uniprot not in Fasta" + '\n')
        error.close
        seqlength = 'NA'
        genename = 'NA'
        return ReturnValue1(seqlength, genename)


def readtab(infile):
    with open(infile,'r') as x: # read in tab-delim text
        output = []
        for line in x:
            line = line.strip()
            temp = line.split('\t')
            output.append(temp)
    return output
def read_maxquant(maxquant_input): # Get data, proteins and header from maxquant output
    dupes = readtab(maxquant_input)
    header_start = 0 
    header = dupes[header_start]
    for i in header:
        i = i.replace(r"\"", "")
        i = i.replace(r"Intensity.", r"")
        i = i.replace(r".", r"-")
    data = dupes[header_start+1:len(dupes)] #cut off blank line and END OF FILE
    proteins = []
    for protein in data:
        proteins.append(protein[0])
    return ReturnValue2(data, proteins, header)
def make_inter(maxquant_input):
    bait = readtab(baitfile)
    data = read_maxquant(maxquant_input).data
    header = read_maxquant(maxquant_input).header
    proteins = read_maxquant(maxquant_input).proteins
    bait_index = []
    for i in bait:
        bait_index.append(header.index("mapped_protein") + 1) # Find just the baits defined in bait file
    with open('inter.txt', 'w') as y:
            a = 0; l=0
            for bb in bait:
                for lst in data:
                    y.write(header[bait_index[l]] + '\t' + bb[1] + '\t' + proteins[a] + '\t' + lst[bait_index[l]] + '\n')
                    a+=1
                    if a == len(proteins):
                        a = 0; l+=1
def make_prey(maxquant_input):
    proteins = read_maxquant(maxquant_input).proteins
    output_file = open("prey.txt",'w')
    for a in proteins:
        a = a.replace("\n","") # remove \n for input into function
        a = a.replace("\r","") # ditto for \r
        seq = get_info(a).seqlength
        GN = get_info(a).genename
        if seq != 'NA':
            output_file.write(a+"\t"+str(seq)+ "\t" + str(GN) + "\n")
    output_file.close()
def no_error_inter(maxquant_input): # remake inter file without protein errors from Uniprot
    err = readtab("error proteins.txt")
    bait = readtab(baitfile)
    data = read_maxquant(maxquant_input).data
    header = read_maxquant(maxquant_input).header #I think this can be done in one line with multiple .replace() or maybe else statements within the comprehension?
    # AB: remember trying and it simply would not work with one line. 
    header = [i.replace(r"\"", "") for i in header]
    header = [i.replace(r"Intensity.", r"") for i in header]
    header = [i.replace(r".", r"-") for i in header]
    bait_index = []
    for i in bait:
        bait_index.append(header.index(i[0]))
    proteins = read_maxquant(maxquant_input).proteins
    errors = []
    for e in err:
        errors.append(e[0])
    with open('inter.txt', 'w') as y:
        l = 0; a = 0
        for bb in bait:
            for lst in data:
                if proteins[a] not in errors:
                    y.write(header[bait_index[l]] + '\t' + bb[1] + '\t' + proteins[a] + '\t' + lst[bait_index[l]] + '\n')
                a+=1
                if a == len(proteins):
                    l += 1; a = 0
def bait_check(bait, maxquant_input): # check that bait names share header titles
    bait_in = readtab(bait)
    header = read_maxquant(maxquant_input).header
    for i in bait_in:
        if i[0] not in header:
            sys.exit("Bait must share header titles with MaxQuant output")

if __name__ == '__main__':
    main(infile, make_bait)
