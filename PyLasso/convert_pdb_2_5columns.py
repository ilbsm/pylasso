#!/usr/bin/python3
import sys
import numpy as np
import argparse
import re
from os import rename, remove
from shutil import copyfile

date = "05.06.2017"
parser = argparse.ArgumentParser(prog="convert_columns", formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description="#################################################################\n\
#	convert_column - script converting PDB to XYZ files.	#\n\
#	Date: 02.10.2015, version from " + date + "		#\n\
#       Author: Pawel Dabrowski-Tumanski			#\n\
#	p.dabrowski [at] cent.uw.edu.pl				#\n\
#	version 2.1						#\n\
#################################################################")
parser.add_argument('input_file', action="store", help="The input PDB file")
parser.add_argument('-t', '--trajectory', action="store_true", dest="traj", default=False,
                    help="Declare, that the input file is a trajectory")
parser.add_argument('-f', '--fourcolumn', action="store_true", dest="fourcolumn", default=False,
                    help="Print XYZ output in 4-column format (default 5-column)")
parser.add_argument('-r', action="store_true", dest="romek", default=False,
                    help=argparse.SUPPRESS)
parser.add_argument('--version', action='version', version='%(prog)s 2.1')

args = parser.parse_args()

global find_index, amino_acids

################################ Possible amino acids ################################
amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO',
               'SER', 'THR', 'TRP', 'TYR', 'VAL', 'BTC', 'FCY', 'GGL']


################################ Functions ################################
def find_index(number, arr):
    for k in range(len(arr)):
        if (arr[k][0] == number):
            return k


def parse_traj(name, out, four):
    f = open(name, 'r')
    got_chain = 0
    nextchain = 0
    new_chain = 1
    art_time = 0
    new_model = 0
    time_changed = 0
    oldtime = 0
    time = 'unset'
    chains = []
    names = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnoprstuvwxyz"
    output = open(out + '_' + str(nextchain) + '.xyz', 'w')

    for line in f:
        if ((line[0:6] == "TITLE ") or (line[0:6] == "REMARK")) and len(
                re.findall("t=[ ]*([0-9]+\.[0-9]+|[0-9]+)", line)) > 0 and got_chain == 0:
            if not output.closed:
                output.close()
            nextchain = 0
            output = open(out + '_' + str(nextchain) + '.xyz', 'a')
            if re.findall("(t=[ ]*[0-9]+\.[0-9]+)|(t=[ ]*[0-9]+)", line)[0][0] != '':
                time = re.findall("(t=[ ]*[0-9]+\.[0-9]+)|(t=[ ]*[0-9]+)", line)[0][0][2:].strip()
                time = "{0:.5f}".format(float(time))
            else:
                time = re.findall("(t=[ ]*[0-9]+\.[0-9]+)|(t=[ ]*[0-9]+)", line)[0][1][2:].strip()
                time = "{0:.5f}".format(float(time))
            got_chain = 1
            new_chain = 1
        if (line[0:6] == "MODEL ") and got_chain == 0:
            if not output.closed: output.close()
            nextchain = 0
            output = open(out + '_' + str(nextchain) + '.xyz', 'a')
            if re.findall("([0-9]+\.[0-9]+)|([0-9]+)", line)[0][0] != '':
                time = re.findall("([0-9]+\.[0-9]+)|([0-9]+)", line)[0][0]
                time = "{0:.5f}".format(float(time))
            else:
                time = re.findall("([0-9]+\.[0-9]+)|([0-9]+)", line)[0][1]
                time = "{0:.5f}".format(float(time))
            got_chain = 1
            new_chain = 1
            new_model = 1
        if (line[0:6] == "TER   "):
            output.close()
            new_chain = 1
            nextchain = nextchain + 1
            output = open(out + '_' + str(nextchain) + '.xyz', 'a')
        if ((line[0:6] == "ATOM  ") or (line[0:6] == "HETATM")) and (line[12:16].strip() == "CA"):
            if time == 'unset':
                art_time = 1
                time = 0
                time = "{0:.5f}".format(float(time))
            if (art_time == 1) and (new_model == 1):
                oldtime = oldtime + 1
                time = oldtime
                time = "{0:.5f}".format(float(time))
                new_model = 0
                time_changed = 1
            if (new_chain == 1):
                if (line[21] in chains) and (art_time == 1) and (time_changed == 0):
                    time = time + 1
                    time = "{0:.5f}".format(float(time))
                output.write("t " + str(time) + "\n")
                new_chain = 0
            if four == False:
                output.write(
                    str(int(line[22:26])) + "  " + str(float(line[30:38])) + " " + str(float(line[38:46])) + " " + str(
                        float(line[46:54])) + " " + str(line[17:20]) + "\n")
            if four == True:
                output.write(
                    str(int(line[22:26])) + "  " + str(float(line[30:38])) + " " + str(float(line[38:46])) + " " + str(
                        float(line[46:54])) + "\n")
            got_chain = 0
            time_changed = 0
            if line[21] not in chains: chains.append(line[21])
        if (line[0:6] == "ENDMDL"): new_model = 1
    if not output.closed: output.close()
    f.close()

    for k in range(len(chains)):
        if chains[k] != ' ':
            rename(out + '_' + str(k) + '.xyz', out + '_' + chains[k] + '.xyz')
        else:
            rename(out + '_' + str(k) + '.xyz', out + '_' + names[k] + '.xyz')
    if chains[0] != ' ':
        copyfile(out + '_' + chains[0] + '.xyz', out + '.xyz')
    else:
        copyfile(out + '_A.xyz', out + '.xyz')
    for i in range(k + 1, nextchain + 1):
        remove(out + '_' + str(i) + '.xyz')


class Chain:
    def __init__(self, name):
        self.name = name
        self.residues = []
        self.coordinates = []
        self.bridges = []
        self.missing = []
        self.gaps = []
        self.helix = []
        self.sheet = []

    def find_length(self):
        self.length = max(len(self.residues), len(self.coordinates))

    def add_residue(self, line):
        line = line[19:70]
        for k in range(len(line.split())):
            if (line.split()[k] in amino_acids):
                self.residues.append(line.split()[k])
        self.find_length()

    def add_bridge(self, bridge):
        self.bridges.append(bridge)

    def add_helix(self, helix):
        self.helix.append(helix)

    def add_sheet(self, sheet):
        self.sheet.append(sheet)

    def add_missing(self, missing, residue):
        self.missing.append([missing, residue])

    def add_coordinate(self, index, coordinate, residue):
        if (len(self.coordinates) == 0):
            for k in range(len(self.missing)):
                if (index > self.missing[k][0]):
                    self.coordinates.append([self.missing[k][0], [], self.missing[k][1]])
            self.coordinates.append([index, coordinate, residue])
        else:
            diff = index - self.coordinates[len(self.coordinates) - 1][0]
            if (diff > 1):
                x = self.coordinates[len(self.coordinates) - 1][1][0]
                y = self.coordinates[len(self.coordinates) - 1][1][1]
                z = self.coordinates[len(self.coordinates) - 1][1][2]
                vec_diff = [(coordinate[0] - x) / diff, (coordinate[1] - y) / diff, (coordinate[2] - z) / diff]
                self.gaps.append([diff, self.coordinates[len(self.coordinates) - 1][0], index])
            for k in range(diff - 1):
                x = round(x + vec_diff[0], 3)
                y = round(y + vec_diff[1], 3)
                z = round(z + vec_diff[2], 3)
                #    res=self.missing[find_index(index-diff+k+1,self.missing)][1]
                #    if (res==None):
                res = "XXX"
                self.coordinates.append([index - diff + k + 1, [x, y, z], res])
            self.coordinates.append([index, coordinate, residue])

    ######### characterize the bond type
    def bond_type(self, res1, atom1, Nend, res2, atom2, Cend):
        if ((atom1[0], atom2[0]) == ("C", "N")):
            if ((((res1, atom1) == ("GLU", "CD")) or ((res1, atom1) == ("ASP", "CG"))) and (
                    (res2, atom2) == ("LYS", "NZ"))):
                return "AMIDE"
            else:
                return "AMIDE-like"
        if ((atom1[0], atom2[0]) == ("N", "C")):
            if ((((res1, atom1) == ("LYS", "NZ")) or (Nend and atom1 == "N")) and (
                    ((res2, atom2) == ("GLU", "CD")) or ((res2, atom2) == ("ASP", "CG")) or (
                    Cend and atom2 == "C"))):
                return "AMIDE"
            else:
                return "AMIDE-like"
        if ((atom1[0], atom2[0]) == ("C", "O")):
            if ((((res1, atom1) == ("GLU", "CD")) or ((res1, atom1) == ("ASP", "CG"))) and (
                    ((res2, atom2) == ("SER", "OG")) or ((res2, atom2) == ("THR", "OG1")))):
                return "ESTER"
            else:
                return "ESTER-like"
        if ((atom1[0], atom2[0]) == ("O", "C")):
            if ((((res1, atom1) == ("SER", "OG")) or ((res2, atom2) == ("THR", "OG1"))) and (
                    ((res2, atom2) == ("GLU", "CD")) or ((res2, atom2) == ("ASP", "CG")) or (
                    Cend and atom2 == "C"))):
                return "ESTER"
            else:
                return "ESTER-like"
        if ((atom1[0], atom2[0]) == ("C", "S")):
            if ((((res1, atom1) == ("GLU", "CD")) or ((res1, atom1) == ("ASP", "CG"))) and (
                    (res2, atom2) == ("CYS", "SG"))):
                return "THIOESTER"
            else:
                return "THIOESTER-like"
        if ((atom1[0], atom2[0]) == ("S", "C")):
            if (((res1, atom1) == ("CYS", "SG")) and (
                    ((res2, atom2) == ("GLU", "CD")) or ((res2, atom2) == ("ASP", "CG")) or (
                    Cend and atom2 == "C"))):
                return "THIOESTER"
            else:
                return "THIOESTER-like"
        else:
            return "OTHER"

    ######## check, whether it is N- or C-end
    def N_end(self, index):
        if (index == self.coordinates[0][0]):
            return True
        else:
            return False

    def C_end(self, index):
        if (index == self.coordinates[len(self.coordinates) - 1][0]):
            return True
        else:
            return False

    ######### cleaning data
    def clean(self):
        ### adding C-end
        if (len(self.coordinates) != 0):
            for k in range(len(self.missing)):
                if (self.coordinates[len(self.coordinates) - 1][0] < self.missing[k]):
                    self.coordinates.append([self.missing[k][0], [], self.missing[k][1]])
        else:
            for k in range(len(self.missing)):
                self.coordinates.append([self.missing[k][0], [], self.missing[k][1]])
                ### clearing double CA atoms
        for k in range(len(self.coordinates) - 1, 0, -1):
            if (self.coordinates[k][0] == self.coordinates[k - 1][0]):
                self.coordinates.pop(k)
                #   if (self.coordinates_residue[k][0]==self.coordinates[k-1][0]):
                #    self.coordinates_residue.pop(k)
                ### clearing non-protein links
        for k in range(len(self.bridges) - 1, -1, -1):
            if ((self.bridges[k][1] not in amino_acids) or (self.bridges[k][4] not in amino_acids)):
                self.bridges.pop(k)
                ### clearing to small loops
        for k in range(len(self.bridges) - 1, -1, -1):
            if (abs(self.bridges[k][3] - self.bridges[k][6]) < 5):
                self.bridges.pop(k)
                ### defining bond type
        for k in range(len(self.bridges)):
            if (self.bridges[k][0] == "LINK"):
                self.bridges[k][0] = self.bond_type(self.bridges[k][1], self.bridges[k][2],
                                                    self.N_end(self.bridges[k][3]), self.bridges[k][4],
                                                    self.bridges[k][5], self.C_end(self.bridges[k][6]))
                ### removing bonds and links which do not exist!
        for k in range(len(self.bridges) - 1, -1, -1):
            if (find_index(self.bridges[k][3], self.coordinates) == None or find_index(self.bridges[k][6],
                                                                                       self.coordinates) == None):
                self.bridges.pop(k)
                ### removing amide bonds used to extend the backbone
        for k in range(len(self.bridges) - 1, -1, -1):
            if ((self.bridges[k][0] == "AMIDE-like") and (abs(
                    find_index(self.bridges[k][3], self.coordinates) - find_index(self.bridges[k][6],
                                                                                  self.coordinates)) == 1)):
                self.bridges.pop(k)
                ### dealing with different number of residues in different parts of PDB file
        k = len(self.residues) - len(self.coordinates)
        if (k > 0):
            check = 1  # check=1, if residues in N-terminus are same
            for i in range(min(k, len(self.residues), len(self.coordinates))):
                if (self.residues[i] != self.coordinates[i][2]):
                    check = 0
            if (check == 0):
                for i in range(k):
                    self.coordinates.insert(0, [self.coordinates[0][0] - 1, [], ""])
            else:
                if (len(self.residues) != 0 and len(self.coordinates) != 0):
                    for i in range(len(self.residues) - k, len(self.residues), 1):
                        self.coordinates.append([self.coordinates[len(self.coordinates) - 1][0] + 1, [], ""])
        if (k < 0):
            check = 1  # check=1, if residues in N-terminus are same
            for i in range(k):
                if (self.residues[i] != self.coordinates[i][2]):
                    check = 0
            if (check == 0):
                for i in range(k):
                    self.residues.insert(0, "UNK")
            else:
                for i in range(len(self.residues), len(self.residues) + k, 1):
                    self.residues.append("UNK")
                    ### length check, just for sure
        self.find_length()
        ######### checking gaps

    def check_gaps(self):
        communicate = "\n"
        for k in range(len(self.gaps)):
            if (self.gaps[k][0] > 2):
                communicate += "WARNING!!! In chain " + self.name + " there is a gap of length " + str(
                    self.gaps[k][0] - 1) + " between residues " + str(self.gaps[k][1]) + " and " + str(
                    self.gaps[k][2]) + "\n"
        communicate = communicate[:-1]
        return communicate

    ######### printing data
    def chain_print(self, PDB, four):
        output_file = open(PDB + "_" + self.name + ".xyz", 'w')
        print(self.check_gaps())
        if (len(self.residues) == self.length):
            for k in range(min(self.length, len(self.coordinates))):
                if (self.coordinates[k][1] != []) and (four == False):
                    output_file.write(str(self.coordinates[k][0]) + " " + str(self.coordinates[k][1][0]) + " " + str(
                        self.coordinates[k][1][1]) + " " + str(self.coordinates[k][1][2]) + " " + str(
                        self.residues[k]) + "\n")
                if (self.coordinates[k][1] != []) and (four == True):
                    output_file.write(str(self.coordinates[k][0]) + " " + str(self.coordinates[k][1][0]) + " " + str(
                        self.coordinates[k][1][1]) + " " + str(self.coordinates[k][1][2]) + "\n")
        else:
            for k in range(self.length):
                if (self.coordinates[k][1] != []) and (four == False):
                    output_file.write(str(self.coordinates[k][0]) + " " + str(self.coordinates[k][1][0]) + " " + str(
                        self.coordinates[k][1][1]) + " " + str(self.coordinates[k][1][2]) + " " + str(
                        self.coordinates[k][2]) + "\n")
                if (self.coordinates[k][1] != []) and (four == True):
                    output_file.write(str(self.coordinates[k][0]) + " " + str(self.coordinates[k][1][0]) + " " + str(
                        self.coordinates[k][1][1]) + " " + str(self.coordinates[k][1][2]) + "\n")
        output_file.close()
        output_file = open(PDB + "_" + self.name + ".pdb", 'w')
        for k in range(len(self.helix)):
            output_file.write(self.helix[k])
        for k in range(len(self.sheet)):
            output_file.write(self.sheet[k])
        if (len(self.residues) == self.length):
            for k in range(min(self.length, len(self.coordinates))):
                if (self.coordinates[k][1] != []):
                    output_file.write(
                        "ATOM  %(atom)5s  CA  %(resname)3s A%(res)4s    %(x)8s%(y)8s%(z)8s  1.00  1.00           C\n" % {
                            "atom": self.coordinates[k][0], "res": self.coordinates[k][0], "resname": self.residues[k],
                            "x": self.coordinates[k][1][0], "y": self.coordinates[k][1][1],
                            "z": self.coordinates[k][1][2]})
        else:
            for k in range(self.length):
                if (self.coordinates[k][1] != []):
                    output_file.write(
                        "ATOM  %(atom)5s  CA  %(resname)3s A%(res)4s    %(x)8s%(y)8s%(z)8s  1.00  1.00           C\n" % {
                            "atom": self.coordinates[k][0], "res": self.coordinates[k][0],
                            "resname": self.coordinates[k][2], "x": self.coordinates[k][1][0],
                            "y": self.coordinates[k][1][1], "z": self.coordinates[k][1][2]})
        output_file.write("END\n")
        output_file.close()

    ######### printing commands to program
    def commands_print(self, PDB, flag=0):
        for k in range(len(self.bridges)):
            if (flag == 1):
                print(self.bridges[k][0] + " ./surfacesMyOrient " + PDB + "_" + self.name + ".xyz " + str(
                    self.bridges[k][3]) + " " + str(self.bridges[k][6]) + " 0 0")
            if (flag == 2):
                print(self.bridges[k][0] + " " + PDB + "_" + self.name + " " + str(self.bridges[k][3]) + " " + str(
                    self.bridges[k][6]))
            else:
                print("./surfacesMyOrient " + PDB + "_" + self.name + ".xyz " + str(self.bridges[k][3]) + " " + str(
                    self.bridges[k][6]) + " 0 0")
                ######### finding distance between residues

    def find_distance(self, res1, res2):
        for k in range(len(self.coordinates)):
            if (self.coordinates[k][0] == res1):
                vec1 = self.coordinates[k][1]
            if (self.coordinates[k][0] == res2):
                vec2 = self.coordinates[k][1]
        print(np.linalg.norm(np.asarray(vec1) - np.asarray(vec2)))

    ######### finding index of first residue with coordinates
    def find_first(self):
        for k in range(len(self.coordinates)):
            if (self.coordinates[k][1] != []):
                return self.coordinates[k][0]


if (args.romek == True):
    import webbrowser

    webbrowser.open("https://www.youtube.com/watch?v=niiYv09hHOI")
    sys.exit(0)

################################ Main part ################################
### search for chains and build chain classes

if args.traj:
    parse_traj(args.input_file, args.input_file, args.fourcolumn)

else:
    chains = []
    input_file = open(args.input_file, 'r')
    ternum = 0
    terfound = 1
    tercount = 1
    art_chains = 0
    names = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnoprstuvwxyz"
    for line in input_file:
        if (line[0:6] == "SEQRES"):
            if (line[11] not in chains):
                chains.append(line[11])
        if (line[0:4] == "ATOM"):
            terfound = 0
            if (line[21] not in chains):
                chains.append(line[21])
        if (line[0:3] == "TER") and (terfound == 0) and (tercount == 1):
            ternum = ternum + 1
            terfound = 1
        if (line[0:3] == "END"):
            if (terfound == 0) and (tercount == 1): ternum = ternum + 1
            tercount = 0
    if (tercount == 1) and (terfound == 0): ternum = ternum + 1
    input_file.close()
    if len(chains) == 1 and chains[0] == ' ' and ternum > 0:
        chains = []
        art_chains = 1
        for k in range(ternum):
            chains.append(names[k])

    ### fill the chain class information
    cross_chain_bridges = []
    chain_data = {}
    ternum = 0
    for k in range(len(chains)):
        chain_data[chains[k]] = Chain(chains[k])
        input_file = open(args.input_file, 'r')
        for line in input_file:
            if ((line[0:10] == "REMARK 465") and (line[15:18] in amino_acids) and (
                    line[19] == chains[k]) and isinstance(
                line[21:26], int)):
                chain_data[chains[k]].add_missing(int(line[21:26]), line[15:18])
            if ((line[0:6] == "SEQRES") and (line[11] == chains[k])):
                chain_data[chains[k]].add_residue(line)
            if ((line[0:5] == "HELIX") and (line[19] == chains[k]) and (line[31] == chains[k])):
                chain_data[chains[k]].add_helix(line)
            if ((line[0:5] == "SHEET") and (line[21] == chains[k]) and (line[32] == chains[k])):
                chain_data[chains[k]].add_sheet(line)
            if ((line[0:6] == "SSBOND") and (line[15] == chains[k]) and (line[29] == chains[k])):
                chain_data[chains[k]].add_bridge(
                    ["SS", line[11:14], "S", int(line[17:21]), line[25:28], "S", int(line[31:35])])
            if ((line[0:6] == "SSBOND") and (line[15] != line[29])):
                cross_chain_bridges.append(
                    ["SS", line[15], line[11:14], "S", int(line[17:21]), line[29], line[25:28], "S", int(line[31:35])])
            if ((line[0:4] == "LINK") and (line[21] == chains[k]) and (line[51] == chains[k])):
                chain_data[chains[k]].add_bridge(
                    ["LINK", line[17:20], line[12:16].strip(), int(line[22:26]), line[47:50], line[42:46].strip(),
                     int(line[52:56])])
            if ((line[0:4] == "LINK") and (line[21] != line[51])):
                cross_chain_bridges.append(
                    ["LINK", line[21], line[17:20], line[12:16].strip(), int(line[22:26]), line[51], line[47:50],
                     line[42:46].strip(), int(line[52:56])])
            if ((line[0:6] == "HETATM") and (line[13:15] == "CA") and (line[21] == chains[k])):
                chain_data[chains[k]].add_coordinate(int(line[22:26]),
                                                     [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                                                     line[17:20])
            if ((line[0:6] == "HETATM") and (line[13:15] == "CA") and (k == ternum) and (art_chains == 1)):
                chain_data[chains[k]].add_coordinate(int(line[22:26]),
                                                     [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                                                     line[17:20])
            if ((line[0:4] == "ATOM") and (line[13:15] == "CA") and (line[21] == chains[k])):
                chain_data[chains[k]].add_coordinate(int(line[22:26]),
                                                     [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                                                     line[17:20])
            if ((line[0:4] == "ATOM") and (line[13:15] == "CA") and (k == ternum) and (art_chains == 1)):
                chain_data[chains[k]].add_coordinate(int(line[22:26]),
                                                     [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                                                     line[17:20])
            if ((line[0:3] == "TER") and (len(chain_data[chains[k]].coordinates) != 0)):
                break
            if ((line[0:3] == "TER") and (len(chain_data[chains[k]].coordinates) == 0)):
                ternum = ternum + 1
        input_file.close()
        chain_data[chains[k]].clean()  # clean bridges data - do not comment
        chain_data[chains[k]].chain_print(args.input_file, args.fourcolumn)  # save coordinates to .xyz and .pdb file
        chain_data[chains[k]].commands_print(args.input_file, 2)
