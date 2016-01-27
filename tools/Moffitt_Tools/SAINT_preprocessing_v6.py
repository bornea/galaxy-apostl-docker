#######################################################################################
# Python-code: SAINT pre-processing from Scaffold "Samples Report" output
# Author: Brent Kuenzi
#######################################################################################
# This program reads in a raw Scaffold "Samples Report" output and a user generated
# bait file and autoformats it into prey and interaction files for SAINTexpress 
# analysis
#######################################################################################
import sys
import urllib2
import os.path
#######################################################################################
## REQUIRED INPUT ##

# 1) infile: Scaffold "Samples Report" output
# 2) baitfile: SAINT formatted bait file generated in Galaxy
# 3) fasta_db: fasta database for use (defaults to SwissProt_HUMAN_2014_08.fasta)
# 4) prey: Y or N for generating a prey file 
# 5) make_bait: String of bait names, assignment, and test or control boolean 
#######################################################################################
infile = sys.argv[1] #Scaffold "Samples Report" output
prey = sys.argv[2] # Y or N
fasta_db = sys.argv[3]
tool_path = r"/galaxy-apostl-docker/tools/Moffitt_Tools/bubblebeam"
if fasta_db == "None":
    fasta_db = str(tool_path)  + "/SwissProt_HUMAN_2014_08.fasta"
make_bait = sys.argv[6]
bait_bool = sys.argv[9]




#What are C and T? 
def bait_create(baits, infile):
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
            line2 = line2.strip()
            temp = line2.split('\t')
            if "Quantitative Variance" in str(temp):
                if baits[i] in temp:
                    number_bait = temp.index(str(baits[i]))
                    number_bait = number_bait - 9 #What's this?
                    # AB: Counter reset?
                    bait_cache.append((number_bait, str(line1)))
                else:
                    print "Error: bad bait " + str(baits[i])
                    sys.exit()
            else: 
                pass
        i = i + 3 #Why iterate by 3? Is there anything that can break this? Rearranged file rows or columns? Can this be soft coded (regex magic) instead of hard coded?
        # AB: It is by 3 because of how the trade off is handled. Galaxy provides this information as a string and I have to split it so the first 3 are the first bait and 4-6 are the second.

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
    print bait_cache
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

def main(scaffold_input, baits): 
    bait_check(baitfile, scaffold_input)
    make_inter(scaffold_input)
    if prey == 'true':
        make_prey(scaffold_input)
        no_error_inter(scaffold_input)
        os.rename('prey.txt', sys.argv[5])
    elif prey == 'false':
        if os.path.isfile('error proteins.txt') == True:
            no_error_inter(scaffold_input)
        pass
    elif prey != 'true' or 'false':
        sys.exit("Invalid Prey Argument: Y or N")

def get_info(uniprot_accession_in): # get aa lengths and gene name
    error = open('error proteins.txt', 'a+')
    data = open(fasta_db,'r')
    lines = data.readlines()
    db_len = len(lines)
    seqlength = 0
    count = 0
    for i in lines:           
        if ">sp" in i:
            namer = i.split("|")[2]
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
            elif uniprot_accession_in == namer.split(" ")[0]:
                match = count+1 #consistent spacing throughout
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
def read_scaffold(scaffold_input): # Get data, proteins and header from scaffold output
    dupes = readtab(scaffold_input)
    cnt = 0
    for i in dupes:
        cnt += 1
        if i[0] == '#': # finds the start of second header
            header_start = cnt-1
    header = dupes[header_start]
    prot_start = header.index("Accession Number")
    data = dupes[header_start+1:len(dupes)-2] # cut off blank line and END OF FILE
    #A lot of text editors and programs will automatically add a blank line at the end of an output file. Do we check or correct for them anywhere else? Maybe want to standardize this.
    proteins = []
    for i in data:
        i[4] = i[4].split()[0] # removes the (+##) that sometimes is attached
    for protein in data:
        proteins.append(protein[prot_start])
    return ReturnValue2(data, proteins, header)
def make_inter(scaffold_input):
    bait = readtab(baitfile)
    data = read_scaffold(scaffold_input).data
    header = read_scaffold(scaffold_input).header
    proteins = read_scaffold(scaffold_input).proteins
    bait_index = []
    for i in bait:
        bait_index.append(header.index(i[0])) # Find just the baits defined in bait file
    with open('inter.txt', 'w') as y:
            a = 0; l=0
            for bb in bait:
                for lst in data:
                    y.write(header[bait_index[l]] + '\t' + bb[1] + '\t' + proteins[a] + '\t' + lst[bait_index[l]] + '\n')
                    a+=1
                    if a == len(proteins):
                        a = 0; l+=1
def make_prey(scaffold_input):
    proteins = read_scaffold(scaffold_input).proteins
    output_file = open("prey.txt",'w')
    for a in proteins:
        a = a.replace("\n","") # remove \n for input into function
        a = a.replace("\r","") # ditto for \r
        seq = get_info(a).seqlength
        GN = get_info(a).genename
        if seq != 'NA':
            output_file.write(a+"\t"+str(seq)+ "\t" + str(GN) + "\n")
    output_file.close()
def no_error_inter(scaffold_input): # remake inter file without protein errors from Uniprot
    err = readtab("error proteins.txt")
    bait = readtab(baitfile)
    data = read_scaffold(scaffold_input).data
    header = read_scaffold(scaffold_input).header
    bait_index = []
    for i in bait:
        bait_index.append(header.index(i[0]))
    proteins = read_scaffold(scaffold_input).proteins
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
def bait_check(bait, scaffold_input): # check that bait names share header titles
    bait_in = readtab(bait)
    header = read_scaffold(scaffold_input).header
    for i in bait_in:
        if i[0] not in header:
            sys.exit("Bait must share header titles with Scaffold output")

if __name__ == '__main__':
    main(infile, baitfile)

os.rename("inter.txt", sys.argv[4])
os.rename("bait.txt", sys.argv[7])
