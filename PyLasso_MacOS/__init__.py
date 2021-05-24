# -*- coding: utf-8 -*-
# PyLasso: A PyMOL plugin to identify lassos, version for Mac platform
# Script/plugin by Aleksandra Gierut (a.gierut@cent.uw.edu.pl)
# Questions should be addressed to: Joanna Sulkowska (jsulkowska@cent.uw.edu.pl)
# Any technical difficulties and remarks please report to: Aleksandra Gierut (a.gierut@cent.uw.edu.pl)
#
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTUOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------


import re
import subprocess
import shutil
import textwrap
import imp
import decimal

import sys

from matplotlib.ticker import FormatStrFormatter

try:
    import Pmw
    import tkinter as tk
    import tkinter.filedialog
except:
    print("  ### Graphic libraries not found. Please install them (Tkinter and Pmw) and re-run the plugin.")
try:
    import matplotlib as mplt
    mplt.use('TKAgg')
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
except:
    print("  ### Matplotlib library not found. Please install it and re-run the plugin.")
from pymol.cgo import *
from pymol import cmd
from tkinter.font import Font
from math import ceil

errors = {
    1: "There 2 Ca atoms with identical indices and position.",
    2: "There are two consecutive atoms in the non-natural distance lower than 2 angstrom or larger than 4.2 angstrom.",
    3: "There are at least two Ca atoms with higher index occurring before a lower one.",
    4: "The input file is empty.",
    5: "The indices of the bridge forming residues were not given, or only one of them was given.",
    6: "The indices of the bridge are wrong.",
    7: "There were to few points in the closed loop.",
    8: "There is no bridge - the distance between bridge forming Ca atoms does not lie in the range 3.1-10 angstrom",
    9: "One of the parameters given to the program is wrong."
}

lasso_description = {
    "L0": " Trivial loop, the closed loop for which there is no essential (not artificial) piercing.",
    "L1": "Single lasso, a covalent loop is pierced once by a tail.",
    "L2": "Double lasso, a covalent loop is pierced twice by the same tail.",
    "L3": "Triple lasso, a covalent loop is pierced three times by the same tail.",
    "L4": "Quadruple lasso, a covalent loop is pierced three times by the same tail.",
    "L5": "Quintuple lasso, a covalent loop is pierced three times by the same tail.",
    "L6": "Sextuplet lasso, a covalent loop is pierced three times by the same tail.",
    "LL1,1": "Two-sided lasso, a covalent loop is pierced by two tails, 1 time by one tail and 1 time by another tail",
    "LL1,2": "Two-sided lasso, a covalent loop is pierced by two tails, 1 time by one tail and 2 times by another tail",
    "LL1,4": "Two-sided lasso, a covalent loop is pierced by two tails, 1 time by one tail and 4 times by another tail",
    "LL2,1": "Two-sided lasso, a covalent loop is pierced by two tails, 2 times by one tail and 1 time by another tail",
    "LL2,2": "Two-sided lasso, a covalent loop is pierced by two tails, 2 times by one tail and 2 times by another tail",
    "LL4,1": "Two-sided lasso, a covalent loop is pierced by two tails, 4 times by one tail and 1 time by another tail",
    "LL4,2": "Two-sided lasso, a covalent loop is pierced by two tails, 4 times by one tail and 2 times by another tail",
    "LL4,3": "Two-sided lasso, a covalent loop is pierced by two tails, 4 times by one tail and 3 times by another tail",
    "LS2": "Supercoiling, one tail pierces the covalent loop, then winds around the protein chain comprising the loop, "
           "and pierces it again.",
    "LS3": "Supercoiling, one tail pierces the covalent loop, then winds around the protein chain comprising the loop, "
           "and pierces it again.",
    "LS4": "Supercoiling, one tail pierces the covalent loop, then winds around the protein chain comprising the loop, "
           "and pierces it again.",
    "LS5": "Supercoiling, one tail pierces the covalent loop, then winds around the protein chain comprising the loop, "
           "and pierces it again.",
    "LS7": "Supercoiling, one tail pierces the covalent loop, then winds around the protein chain comprising the loop, "
           "and pierces it again.",
    "ERR": "No picture for such lasso type. Please report their absence."
}

lassos = ['L0', 'L+1N', 'L-1N', 'L+1C', 'L-1C', 'L+2N', 'L-2N', 'L+2C', 'L-2C', 'L+3N', 'L-3N', 'L+3C', 'L-3C',
          'LL+1,+1', 'LL+1,-1', 'LL-1,+1', 'LL-1,-1', 'LL+1,+2', 'LL+1,-2', 'LL-1,+2', 'LL-1,-2', 'LL+2,+1',
          'LL+2,-1', 'LL-2,+1', 'LL-2,-1', 'LL+4,+2', 'LL+4,-2', 'LL-4,+2', 'LL-4,-2', 'LL+2,+4', 'LL+2,-4',
          'LL-2,+4', 'LL-2,-4', 'LS2++N', 'LS2--N', 'LS2++C', 'LS2--C', 'LS3++-N', 'LS3+--N', 'LS3-++N',
          'LS3--+N', 'LS3++-C', 'LS3+--C', 'LS3-++C', 'LS3--+C', 'ERR', 'Other']

colors = ['#FF4000', '#FF8000', '#FFBF00', '#FFFF00', '#BFFF00', '#80FF00', '#40FF00', '#00FF00',
          '#00FF40', '#00FF80', '#00FFBF', '#00FFFF', '#00BFFF', '#0080FF', '#0040FF', '#0000FF', '#4000FF',
          '#8000FF', '#BF00FF', '#FF00FF', '#FF00BF', '#FF0080', '#FF0040', '#CC9900', '#999966', '#000000',
          '#FF9999', '#FFCC99', '#FFFF99', '#B3FF99', '#99FFB3', '#99FFCC', '#9999FF', '#FF99FF', '#FF99B3',
          '#FF3399', '#990033', '#CC3300', '#FFFFFF', '#9900FF', '#0033CC', '#339966', '#336600', '#cc3300',
          '#CC6699', 'red', '#FF9966']

plugin_path = os.path.dirname(__file__)
sys.path.append(plugin_path)
system_working_directory = os.getcwd()


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command', 'PyLasso',
                             label='PyLasso',
                             command=lambda s=self: PyLasso(s))


class PyLasso:
    def __init__(self, app):

        self.parent = app.root
        if cmd.get_version()[1] < 2.0:
            self.is_pymol_2 = False
            #self.parent.option_add("*Font", "Arial 10")
            self.bold_font = Font(family="Helvetica", size=11, weight="bold")
            self.subscript_font = Font(family="Helvetica", size=9)
            self.row_name_size = "Helvetica 13"
            self.hull_width = 1000
        else:
            self.is_pymol_2 = True
            self.parent.option_add("*Font", "Helvetica 10")
            self.bold_font = Font(family="Helvetica", size=9, weight="bold")
            self.subscript_font = Font(family="Helvetica", size=9)
            self.row_name_size = "Helvetica 10"
            self.hull_width = 1150

        # Trajectory and single structure variables
        self.screen_height = self.parent.winfo_screenheight()
        self.screen_width = self.parent.winfo_screenwidth()
        self._filename = ""
        self._file_extension = ""
        self._img_extension = ".gif"
        self.program_execution = plugin_path + os.sep + "detect_lassos "
        self.python_compiler = [plugin_path + os.sep + 'convert_pdb_2_5columns.py']
        self.img_button_height = 22
        self.img_button_width = 75
        self.view_btn_width = 3
        self.extend_clear_btns_width = 4
        self.retrieve_btns_width = 9
        self.retrieve_btns_pad = 1
        self.table_name_width = 14
        self.gln_figsize = 5
        self.gln_dpi = 98
        self.lassos = []
        self.hint_width = 60
        # Advanced variables
        self.is_stable = tk.IntVar()
        self.is_bad_caca_enabled = tk.IntVar()
        self.is_gln_checkbutton_selected = tk.IntVar()

        self.previous_bond_in_view = ["", ""]

        self.load_file()

    ####################################################################################################################
    #                             CHECK FILE EXTENSION & ADJUST POLYMER REPRESENTATION
    ####################################################################################################################

    def load_file(self):
        self.open_file_window = tkinter.filedialog.askopenfile(initialdir=os.getcwd(), title="PyLasso",
                                                         filetypes=(("PDB", "*.pdb"), ("XYZ", "*.xyz")))

        if self.open_file_window == None:
            print("  ### No file was chosen. PyLasso was shut down.")
            return

        self._full_path_to_file = self.open_file_window.name

        self._full_path_to_dir = os.sep.join(self._full_path_to_file.split(os.sep)[:-1])
        self.is_original_pdb = False
        self.is_trajectory = False
        self.is_artifact = False

        cmd.reinitialize()
        cmd.load(filename=self._full_path_to_file)

        self._file_extension = self._full_path_to_file[-3:]
        self._filename = self._full_path_to_file.split(os.sep)[-1]
        
        self.chains = ["A"] if (len(cmd.get_chains()) < 2) else cmd.get_chains()
        
        self.check_if_trajectory()
        if not self.is_trajectory:
            self.is_original_pdb = self.contains_bridge_information()

        if self._file_extension == "xyz":
            self._filename = self._filename[:-4] + "_xyz2.pdb"
            self.convert_trajectory_xyz_to_pdb() if self.is_trajectory else self.convert_protein_xyz_to_pdb()
            self._full_path_to_file = self._full_path_to_dir + os.sep + self._filename

        self.initialise_plugin_interface()
        self.adjust_object_representation()
        self.delete_pymol_objects()
        print("  (Bio)Polymer reloaded...")

    def raise_popup_menu(self, error_message):
        if hasattr(self, "error_pop_menu") and self.error_popup.winfo_exists():
            self.error_popup.withdraw()

        self.error_popup = Pmw.MessageDialog(self.parent, title='Error!', defaultbutton=0,
                                             message_text=textwrap.fill(error_message, 80))
        self.error_popup.geometry("+%d+%d" % (self.screen_width / 2 - 150, self.screen_height / 2))

        self.error_popup.focus_force()
        self.error_popup.wait_window()

    def check_if_trajectory(self):
        with open(self._full_path_to_file, "r") as traj_file:
            if not self._file_extension == "pdb":
                first_line = traj_file.readline()
                if first_line.__contains__("t"):
                    self.initialize_trajectory_variables()
                    return
            else:
                for line in traj_file:
                    if line.__contains__("SOLUTION NMR"):
                        self.initialize_nmr_messagebox()
                        return
                    if line.__contains__("ENDMDL"):
                        self.initialize_trajectory_variables()
                        return
            self.initialize_protein_variables()

    def initialize_nmr_messagebox(self):
        nmr_message = "Multichain NMR has been detected. Would you like to treat it as a single structure or a " \
                      "trajectory?"
        self.nmr_messagebox = Pmw.MessageDialog(self.parent, title='NMR detected!', defaultbutton=0,
                                                buttons=["A single structure", "Trajectory"],
                                                message_text=textwrap.fill(nmr_message, 80),
                                                command=self._invoke_nmr)
        self.nmr_messagebox.geometry("+%d+%d" % (self.screen_width / 2 - 150, self.screen_height / 2))
        self.nmr_messagebox.focus_force()
        self.nmr_messagebox.wait_window()

    def _invoke_nmr(self, btn):
        if btn == "Trajectory":
            self.initialize_trajectory_variables()
        elif btn == "A single structure":
            self.initialize_protein_variables()
        self.nmr_messagebox.destroy()

    def initialize_trajectory_variables(self):
        self.is_trajectory = True
        self.is_detailed_alg = tk.IntVar()
        self.is_detailed_out = tk.IntVar()
        self.is_detailed_out_frame = tk.IntVar()

        if self._file_extension == "pdb":
            chain = self.chains[0]
            atoms = {'atoms': []}
            cmd.iterate_state(state=1, selection="all and chain " + chain, expression="atoms.append(resi)", space=atoms)
            atoms = atoms['atoms']
            self.marginal_atoms = [int(atoms[0]), int(atoms[-1])]

    def initialize_protein_variables(self):
        self.is_trajectory = False
        self.number_of_own_loops = 3
        self.type_bridge_exists = False
        self.bridge_selected = None
        self.bridge_button_list = []
        self.type_closing_choice = [0, 0]  # choice between defined type of bridge and own loops
        self.number_of_own_loops = 3

    def contains_bridge_information(self):
        input_file = open(self._full_path_to_file, "r").read()
        ss_bonds = list(re.findall('SSBOND|LINK', input_file, flags=re.M | re.S))
        return len(ss_bonds) is not 0

    def convert_trajectory_xyz_to_pdb(self):
        traj_xyz2pdb = open(self._full_path_to_dir + os.sep + self._filename, "w")
        num_frames = 1
        atom_idx = 1

        with open(self._full_path_to_file) as f:
            for idx, line in enumerate(f):
                if line.__contains__("t"):
                    frame_time = line.split(" ")[1]
                    if num_frames == 1:
                        traj_xyz2pdb.write(
                            "REMARK    GENERATED BY PYLASSO\nTITLE    t= %s\nMODEL        %d\n" %
                            (frame_time, num_frames))
                    else:
                        traj_xyz2pdb.write("TER\nENDMDL\nREMARK    GENERATED BY PYLASSO\nTITLE    t= %s\nMODEL        "
                                           "%d\n" % (frame_time, num_frames))
                    num_frames += 1
                    atom_idx = 1
                else:
                    elems = list(filter(len, line.split(" ")))
                    line = "ATOM%7d  CA  GLY A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n" % \
                           (atom_idx, atom_idx, float(elems[1]), float(elems[2]), float(elems[3][:-1]))
                    atom_idx += 1
                    traj_xyz2pdb.write(line)
        traj_xyz2pdb.write("TER\nENDMDL")
        traj_xyz2pdb.close()

    def convert_protein_xyz_to_pdb(self):
        output_file = open(self._full_path_to_dir + os.sep + self._filename, 'w')

        xyz_atm = re.compile(
            r"^\s*(?P<resid>[0-9]*)\s+(?P<x>-?\d+\.\d*)\s+(?P<y>-?\d+\.\d*)\s+(?P<z>-?\d+\.\d*).*$")
        with open(self._full_path_to_file) as f:
            for line in f:
                data = xyz_atm.match(line)
                if data:
                    resid = int(data.groups()[0])
                    x = float(data.groups()[1])
                    y = float(data.groups()[2])
                    z = float(data.groups()[3])
                    output_file.write("ATOM%7d  CA  GLY A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n" % (
                        resid, resid, x, y, z))
        output_file.close()

    def delete_pymol_objects(self):
        if not self.is_trajectory:
            chain = self.chain_index.get()
            cmd.delete(name="CHAIN_" + chain)
        cmd.delete(name="TMP_*")
        cmd.delete(name="DIST*")
        cmd.delete(name="BR_*")
        cmd.delete(name="SEQ")
        cmd.delete(name="SMOOTH_CHAIN_*")
        cmd.delete(name="PIERC")
        cmd.delete(name="TRIANG")
        cmd.delete(name="NEG_*")
        cmd.delete(name="POS_*")
        cmd.delete(name="SHALLOW*")
        cmd.hide(representation="spheres", selection="all")

        if self.previous_bond_in_view[0] is not "" and self.previous_bond_in_view[1] is not "":
            cmd.unbond(atom1=self.previous_bond_in_view[0], atom2=self.previous_bond_in_view[1])

    def adjust_object_representation(self):
        self.get_marginal_atoms()
        if self._file_extension == "xyz" or self.is_trajectory or not self.is_original_pdb:
            cmd.delete(name=self._filename[:-9] + "*")
            cmd.load(filename=self._full_path_to_dir + os.sep + self._filename)
            self.connect_xyz_points()
            cmd.show(representation="lines", selection=self._filename[:-9])
            cmd.set(name="line_width", value="4")
            if self.is_trajectory and len(self.chains) >= 2:
                self.display_first_chain_in_trajectory()
        else:
            cmd.remove(selection="solvent")
            cmd.hide(representation="lines", selection="all")
            cmd.show(representation="cartoon", selection="all")
            cmd.cartoon(type="tube", selection="all")
        cmd.spectrum(palette="rainbow", selection="all")

    def get_marginal_atoms(self):
        atoms = []
        reg = re.compile('ATOM\s\s+\d+')

        with open(self._full_path_to_dir + os.sep + self._filename, "r") as f:
            for line in f:
                clear_line = list(filter(len, line.split(" ")))
                if reg.match(line) and self.is_trajectory:
                    try:
                        atoms.append(int(clear_line[4]))
                    except:
                        atoms.append(clear_line[5])
                elif reg.match(line) and not self.is_trajectory:
                    if clear_line[4].isalpha():
                        atoms.append(clear_line[5])
                    else:
                        idx = clear_line[4][1:]
                        atoms.append(idx)

        self.marginal_atoms = [int(atoms[0]), int(atoms[-1])]

    def connect_xyz_points(self):
        """
            Due to the fact, that file with .xyz extension represents a set of points in space that are by no means
            connected, a plugin iteratively creates a bond between atom(i) and atom(i+1). Thereby one can clearly see a
            chain, a lasso, etc. Method does not apply to file with .pdb extension.
        """
        for i in range(self.marginal_atoms[0], self.marginal_atoms[1]):
            cmd.bond(atom1="id " + str(i) + " and name ca",
                     atom2="id " + str(i + 1) + " and name ca")

    def display_first_chain_in_trajectory(self):
        hide_chains = "all and "
        for i in self.chains[1:]:
            hide_chains += "chain " + i + " and "
        hide_chains = hide_chains[:-5]
        cmd.hide(representation="lines", selection=hide_chains)

    def _invoke_plugin_action(self, clicked_btn):
        if clicked_btn == "Proceed":
            print("  PyLasso is running...")
            self._invoke_program()
        else:
            self.dialog.withdraw()
            if hasattr(self, "error_pop_menu") and self.error_popup.winfo_exists():
                self.error_popup.withdraw()
            if hasattr(self, "win_lasso_info") and self.win_lasso_info.winfo_exists():
                self.win_lasso_info.withdraw()
            if hasattr(self, "win_trajectory_") and self.win_trajectory_.winfo_exists():
                self.win_trajectory_.withdraw()
            if hasattr(self, "artifact_found") and self.artifact_found.winfo_exists():
                self.artifact_found.withdraw()
            if hasattr(self, "bridge_found") and self.bridge_found.winfo_exists():
                self.bridge_found.withdraw()
            print("  PyLasso has been shut down...")

    ####################################################################################################################
    #                                            INITIALISE INTERFACE
    ####################################################################################################################

    def initialise_plugin_interface(self):
        self.dialog = Pmw.Dialog(self.parent, buttons=('Proceed', 'Exit'), title='PyLasso [' + self._filename + ']',
                                 command=self._invoke_plugin_action)
        self.dialog.geometry("+%d+%d" % (self.screen_width / 4, self.screen_height / 4))

        self.main_window = tk.Label(self.dialog.interior())
        self.main_window.pack(fill='both', expand=1, padx=5, pady=5)

        self.group_advanced = Pmw.Group(self.main_window, tag_text='Advanced options')
        self.group_advanced.grid(sticky='eswn', column=1, row=0, padx=5, pady=5, rowspan=3)

        self.create_trajectory_interface() if self.is_trajectory else self.create_protein_interface()
        self.create_smooth_interior()
        self.create_correct_mode_interior()
        if not self.is_trajectory:
            self.create_gln_interior()
        self.create_advanced_frame_hints()

        Pmw.setbusycursorattributes(self.dialog.component('hull'))
        self.dialog.resizable(0, 0)
        self.dialog.show()

    ####################################################################################################################
    #                                          1) TRAJECTORY INTERFACE
    ####################################################################################################################

    def create_trajectory_interface(self):
        print("  Trajectory detected. Setting correct interface...")
        if hasattr(self, "fr_selected_loop") and hasattr(self, "fr_detailed_output") and hasattr(self, "fr_step"):
            self.fr_selected_loop.grid_forget()
            self.fr_detailed_output.grid_forget()
            self.fr_step.grid_forget()

        cmd.set(name="mouse_selection_mode", value=6)
        cmd.set(name="seq_view", value=1)
        self.distance_error = {}

        self.fr_selected_loop = Pmw.Group(self.main_window, tag_text="Selected loop")
        self.fr_selected_loop.grid(sticky='eswn', column=0, row=0, padx=5, pady=5)

        self.fr_detailed_output = Pmw.Group(self.main_window, tag_text="Accuracy of calculations")
        self.fr_detailed_output.grid(sticky='eswn', column=0, row=1, padx=5, pady=5)

        self.fr_step = Pmw.Group(self.main_window, tag_text="Step")
        self.fr_step.grid(sticky='eswn', column=0, row=2, padx=5, pady=5)

        self.label_trajectory_loop = tk.Label(self.fr_selected_loop.interior(), justify="left",
                                              text=textwrap.fill('Define the first and the last index of residue for '
                                                                 'loop closing. Only one loop is available in the '
                                                                 'trajectory analysis mode.', 50))
        self.label_trajectory_loop.grid(column=0, columnspan=5, row=0, pady=2, padx=5)

        self.loops_list = []
        self.loops_list.append([Pmw.EntryField(self.fr_selected_loop.interior(), labelpos='w', entry_width=8,
                                               validate={'validator': 'integer'}),
                                Pmw.EntryField(self.fr_selected_loop.interior(), labelpos='w', entry_width=8,
                                               validate={'validator': 'integer'})])
        self.loops_list[0][0].grid(sticky='w', column=1, row=2)
        self.loops_list[0][1].grid(sticky='w', column=2, row=2)

        self.btns_view = []
        self.btns_view.append(tk.Button(self.fr_selected_loop.interior(), text="View", width=self.view_btn_width,
                                        command=lambda: self.view_possible_trajectory()))
        self.btns_view[0].grid(sticky='w', column=3, row=2, padx=2)

        self.btns_extend_clear = Pmw.ButtonBox(self.fr_selected_loop.interior(), orient="horizontal")
        self.btns_extend_clear.add("Clear", width=self.extend_clear_btns_width + 4,
                                   command=lambda: self.clear_trajectory())
        self.btns_extend_clear.grid(sticky='swn', column=1, row=3, columnspan=2, padx=28, pady=2)

        tmp_frame = tk.Frame(self.fr_selected_loop.interior(), padx=0, pady=0)
        self.btns_get_data = Pmw.ButtonBox(tmp_frame, orient="vertical", pady=self.retrieve_btns_pad)
        self.btns_get_data.add("Get from\nstructure", width=self.retrieve_btns_width, font=self.bold_font,
                               command=lambda: self.get_loop_from_trajectory_structure())
        self.btns_get_data.add("Get from\nsequence", width=self.retrieve_btns_width, font=self.bold_font,
                               command=lambda: self.get_loop_from_trajectory_sequence())
        self.btns_get_data.grid(sticky='wen', column=0, row=0)
        tmp_frame.grid(column=0, row=2, rowspan=4)

        self.more_detailed_output = tk.Checkbutton(self.fr_detailed_output.interior(), text='More detailed output file',
                                                   variable=self.is_detailed_out)
        self.more_detailed_output.grid(sticky='w', column=0, row=0, padx=2, pady=2)

        self.more_detailed_algorithm = tk.Checkbutton(self.fr_detailed_output.interior(),
                                                      text='More detailed algorithm',
                                                      variable=self.is_detailed_alg)
        self.more_detailed_algorithm.grid(sticky='w', column=0, row=1, padx=2, pady=2)

        self.step = Pmw.EntryField(self.fr_step.interior(), labelpos='w', entry_width=8,
                                   validate={'validator': 'integer'})
        self.step.grid(sticky='eswn', column=0, row=0, padx=2, pady=2)

        self.create_trajectory_interface_hints()

    def view_possible_trajectory(self):
        self.validate_trajectory()

        for i in cmd.get_names(type="all"):
            if str(i).startswith("TMP") or str(i).startswith("DIST_") or str(i).startswith("SMOOTH_CHAIN_*") or \
                    str(i).startswith("SEQ") or str(i).startswith("PIERC") or str(i).startswith("TRIANG") \
                    or str(i).startswith("BR") or str(i).startswith("SMOOTH"):
                cmd.hide(representation="everything", selection=i)
                cmd.delete(name=i)

        if self.previous_bond_in_view[0] is not "" and self.previous_bond_in_view[1] is not "":
            cmd.unbond(atom1=self.previous_bond_in_view[0], atom2=self.previous_bond_in_view[1])
        self.simplify_polymer_representation()

        atom = "ca"
        chain = self.chains[0]
        cmd.set(name="seq_view", value=1)
        cmd.set(name="mouse_selection_mode", value=6)
        cmd.hide(representation="sticks", selection="all")
        cmd.spectrum(palette="rainbow", selection="all")
        cmd.deselect()

        cmd.select("TMP_BR_" + self.loops_list[0][0].getvalue() + "_" + self.loops_list[0][1].getvalue(),
                   selection="(chain " + chain + " and residue " + self.loops_list[0][0].getvalue() + " and name "
                             + atom + ")+(chain " + chain + " and residue " + self.loops_list[0][1].getvalue() +
                             " and name " + atom + ")")
        cmd.select("TMP_SEQ",
                   selection="chain " + chain + " and residue " + self.loops_list[0][0].getvalue() + "-" +
                             self.loops_list[0][1].getvalue())

        cmd.bond(atom1="chain " + chain + " and residue " + self.loops_list[0][0].getvalue() + " and name " + atom,
                 atom2="chain " + chain + " and residue " + self.loops_list[0][1].getvalue() + " and name " + atom)
        cmd.color(color="gray", selection="TMP_SEQ")
        cmd.show(representation='sphere', selection="TMP_BR_" + self.loops_list[0][0].getvalue() + "_" +
                                                    self.loops_list[0][1].getvalue())
        cmd.show(representation="sticks", selection="TMP_BR_" + self.loops_list[0][0].getvalue() + "_" +
                                                    self.loops_list[0][1].getvalue())
        cmd.set(name="sphere_color", value="orange", selection="all")
        cmd.set(name="stick_color", value="orange", selection="all")
        cmd.set("sphere_scale", value=0.5)

        self.check_distance_between_atoms([self.loops_list[0][0], self.loops_list[0][1]], 0)
        self.previous_bond_in_view = ["chain " + chain + " and residue " + self.loops_list[0][0].getvalue()
                                      + " and name " + atom,
                                      "chain " + chain + " and residue " + self.loops_list[0][1].getvalue()
                                      + " and name " + atom]

    def validate_trajectory(self):
        if (len(self.loops_list[0][0].getvalue()) is 0 and len(self.loops_list[0][1].getvalue()) is not 0) or (
                        len(self.loops_list[0][0].getvalue()) is not 0 and len(self.loops_list[0][1].getvalue()) is 0):
            self.raise_popup_menu('There are missing data in fields.')
        if len(self.loops_list[0][0].getvalue()) is 0 or len(self.loops_list[0][1].getvalue()) is 0:
            self.raise_popup_menu('Please fill in the fields in the Selected loop.')
        if (int(self.loops_list[0][0].getvalue()) < self.marginal_atoms[0]) or \
                (int(self.loops_list[0][1].getvalue()) > self.marginal_atoms[1]):
            self.raise_popup_menu('There are no such indices in the structure. Please choose other (correct) ones.')

    def simplify_polymer_representation(self):
        if self._file_extension == "xyz" or self.is_trajectory or not self.is_original_pdb:
            cmd.delete(name=self._filename[:-9] + "*")
            cmd.load(filename=self._full_path_to_file)
            self.connect_xyz_points()
            cmd.show(representation="lines", selection=self._filename[:-9])
            cmd.set(name="line_width", value="4", selection=self._filename[:-9])
            if self.is_trajectory and len(self.chains) >= 2:
                self.display_first_chain_in_trajectory()
        else:
            cmd.load(filename=self._full_path_to_file)
            cmd.remove(selection="solvent")
            cmd.hide(representation="lines", selection="all")

    def clear_trajectory(self):
        self.loops_list[0][0].setentry("")
        self.loops_list[0][1].setentry("")

        for i in list(self.distance_error.keys()):
            self.distance_error[i].grid_forget()

        self.label_trajectory_loop.configure(text=textwrap.fill(self.label_trajectory_loop.cget("text"), 50))
        self.label_trajectory_loop.grid_configure(columnspan=5)
        self.distance_error = {}

    def get_loop_from_trajectory_structure(self):
        atoms = {'atoms': []}
        atom = "ca"

        if len(cmd.get_names(type="selections")) == 0:
            self.raise_popup_menu('No atoms has been chosen.')

        cmd.iterate_state(state=0, selection="sele and name " + atom, expression='atoms.append(resi)', space=atoms)
        atoms = {'atoms': list(set(atoms['atoms']))}
        atoms['atoms'].sort(key=int)

        if len(atoms['atoms']) > 2:
            self.raise_popup_menu('Too many selection (' + str(len(atoms['atoms'])) + ').')
        elif len(atoms['atoms']) < 2:
            self.raise_popup_menu('Too little selections or there are selected atoms from wrong chain.')

        self.loops_list[0][0].setentry(atoms['atoms'][0])
        self.loops_list[0][1].setentry(atoms['atoms'][1])

        if len(cmd.get_names(type="selections")) > 0:
            cmd.delete(name="sele")

    def get_loop_from_trajectory_sequence(self):
        atoms = {'atoms': []}
        atom = "ca"

        if len(cmd.get_names(type="selections")) == 0:
            self.raise_popup_menu('No atoms has been chosen.')

        cmd.iterate_state(state=0, selection="sele and name " + atom, expression='atoms.append(resi)', space=atoms)
        atoms = {'atoms': list(set(atoms['atoms']))}
        atoms['atoms'].sort(key=int)

        if len(atoms['atoms']) < 2:
            self.raise_popup_menu('Too little selections or there are selected atoms from wrong chain.')

        self.loops_list[0][0].setentry(atoms['atoms'][0])
        self.loops_list[0][1].setentry(atoms['atoms'][-1])

        if len(cmd.get_names(type="selections")) > 0:
            cmd.delete(name="sele")

    def create_trajectory_interface_hints(self):
        hint_trajectory = Pmw.Balloon(self.main_window, relmouse="both")
        hint_trajectory.bind(self.fr_detailed_output, "Information about the specificity of data.")
        hint_trajectory.bind(self.more_detailed_output, "Use more detailed output file of trajectory.")
        hint_trajectory.bind(self.more_detailed_algorithm, "Use more detailed algorithm for each frame of trajectory.")
        hint_trajectory.bind(self.step, textwrap.fill("Choose how many frames from the original trajectory should "
                                                      "be skipped between two consecutive frames analyzed by the "
                                                      "PyLasso.", self.hint_width))
        hint_trajectory.bind(self.btns_get_data.button(0), textwrap.fill("Retrieve the selected tuple of C-alpha "
                                                                         "atoms from the (bio)polymer loaded in "
                                                                         "PyMOL.", self.hint_width))
        hint_trajectory.bind(self.btns_get_data.button(1), textwrap.fill("Retrieve the selected sequence of C-alpha"
                                                                         " atoms from the (bio)polymer loaded in "
                                                                         "PyMOL.", self.hint_width))
        hint_trajectory.bind(self.btns_extend_clear.button(0), "Clear the above data.")
        hint_trajectory.bind(self.btns_view[0], "View the loop closed by the the bridge chosen.")

    ####################################################################################################################
    #                                           2) PROTEIN INTERFACE
    ####################################################################################################################

    def create_protein_interface(self):
        print("  Single structure protein detected. Setting correct interface...")
        if hasattr(self, "fr_selected_loop"):
            self.fr_selected_loop.grid_forget()
            self.fr_detailed_output.grid_forget()
            self.fr_step.grid_forget()

        self.main_fr_bridge = tk.Label(self.main_window)
        self.main_fr_bridge.grid(sticky='eswn', column=0, row=0)

        self.fr_chain = Pmw.Group(self.main_fr_bridge, tag_text='Chain information')
        self.fr_chain.grid(sticky='eswn', column=0, row=0, padx=5, pady=5)

        self.fr_loop_closing_bridge = Pmw.Group(self.main_fr_bridge, tag_text='Method to set a loop closing')
        self.fr_loop_closing_bridge.grid(sticky='eswn', column=0, row=1, rowspan=3, padx=5, pady=5)

        tmp_fr_chain = tk.Frame(self.fr_chain.interior(), width=40)
        tmp_fr_chain.grid(column=0, row=0)

        self.chain_index = Pmw.ComboBox(tmp_fr_chain, labelpos='nw', label_text='Chain index',
                                        scrolledlist_items=tuple(self.chains),
                                        selectioncommand=lambda x: self.pymol_display_chain())
        self.chain_index.component("entryfield").insert(0, self.chains[0])
        self.chain_index.configure(entryfield_entry_state="readonly")
        self.chain_index.grid(sticky='swn', column=0, row=0, pady=5, padx=5)

        self.type_loop_closing_bridge = Pmw.ComboBox(self.fr_loop_closing_bridge.interior(), entry_width=30,
                                                     scrolledlist_items=('automatic detections of closed loops',
                                                                         'choose a type of loop closing',
                                                                         'choose two atoms to form a bridge'),
                                                     entryfield_entry_state="readonly",
                                                     selectioncommand=lambda x: self.set_bridge_closing())
        self.type_loop_closing_bridge.selectitem(0)
        self.type_loop_closing_bridge.grid(sticky='eswn', column=0, row=1, padx=5, pady=5)

        self.create_protein_interface_hints()

    def set_bridge_closing(self):
        if str(self.type_loop_closing_bridge.get()) == "automatic detections of closed loops":
            self.choose_automatic_closing()
        elif str(self.type_loop_closing_bridge.get()) == "choose a type of loop closing":
            self.choose_type_closing()
        elif str(self.type_loop_closing_bridge.get()) == "choose two atoms to form a bridge":
            self.choose_own_closing()

    def choose_automatic_closing(self):
        cmd.set(name="seq_view", value=0)

        if hasattr(self, "chosen_type_lc"):
            self.chosen_type_lc.grid_forget()
            self.type_closing_choice[0] = 0
            for elem in self.bridge_button_list:
                elem.deselect()
                elem.grid_forget()
            self.bridge_selected = None
        if hasattr(self, "chosen_own_lc"):
            self.chosen_own_lc.grid_forget()
            self.type_closing_choice[1] = 0

    ####################################################################################################################
    #                                       2b) CHOOSE TYPE OF CLOSING BRIDGE
    ####################################################################################################################

    def choose_type_closing(self):
        if hasattr(self, "chosen_own_lc"):
            self.chosen_own_lc.grid_forget()
            self.type_closing_choice[1] = 0
        if self.type_closing_choice[0] == 0:
            self.type_closing_choice[0] = 1
            self.chosen_type_lc = tk.Label(self.fr_loop_closing_bridge.interior())
            self.chosen_type_lc.grid(column=0, row=2, padx=3)
            cmd.set(name="seq_view", value=0)

            list_img = []
            full_path = plugin_path + os.sep + "img"
            for img in ("SS", "Amide", "Ester", "Thioester", "Other"):
                list_img.append(tk.PhotoImage(file=os.path.join(full_path, img + self._img_extension)))
            tmp_label = tk.Label(self.fr_loop_closing_bridge.interior())  # without this label, pics will disappear!
            tmp_label.image = list_img

            self.bridge_button_list = []
            for idx, text in enumerate(('SS', 'Amide', 'Ester', 'Thioester', 'Others')):
                self.bridge_button_list.append(tk.Radiobutton(self.chosen_type_lc, text=text, image=list_img[idx],
                                                              compound=tk.LEFT, anchor="w", value=idx, indicatoron=0,
                                                              width=self.img_button_width,
                                                              height=self.img_button_height,
                                                              padx=5, command=lambda x=text: self.set_bridge(x)))
                self.bridge_button_list[-1].grid(column=0, row=idx + 2, padx=2, pady=5)
                self.bridge_button_list[-1].deselect()
            self.create_bridge_type_hints()

    def set_bridge(self, idx):
        self.bridge_selected = idx

    def create_bridge_type_hints(self):
        for i in self.bridge_button_list:
            hint_choose_bridge_type = Pmw.Balloon(i, relmouse="both")
            hint_choose_bridge_type.bind(i, textwrap.fill("Choose the type of bridges (disulfide, amide, ester, "
                                                          "thioester, etc.) to be detected by the PyLasso. If no "
                                                          "bridge of the chosen type is found in a given protein "
                                                          "chain, the PyLasso identifies other bridges "
                                                          "automatically.", self.hint_width))

    ####################################################################################################################
    #                                       2c) CHOOSE OWN CLOSING BRIDGE
    ####################################################################################################################

    def choose_own_closing(self):
        if hasattr(self, "chosen_type_lc"):
            self.chosen_type_lc.grid_forget()
            self.type_closing_choice[0] = 0

            for elem in self.bridge_button_list:
                elem.deselect()
                elem.grid_forget()
            self.bridge_selected = None
        if self.type_closing_choice[1] == 0:
            cmd.set(name="mouse_selection_mode", value=6)
            cmd.set(name="seq_view", value=1)
            self.loops_list = []
            self.distance_error = {}
            self.type_closing_choice[1] = 1

            self.chosen_own_lc = tk.Frame(self.fr_loop_closing_bridge.interior())
            self.label_ll = tk.Label(self.chosen_own_lc, justify="left", text=textwrap.fill(
                'Residues can be entered in fields manually or retrieved from PyMOL via Get buttons.', 50))
            self.label_ll.grid(column=0, columnspan=4, row=1, pady=2)

            for i in range(3):
                field = [
                    Pmw.EntryField(self.chosen_own_lc, labelpos='w', entry_width=8, validate={'validator': 'integer'}),
                    Pmw.EntryField(self.chosen_own_lc, labelpos='w', entry_width=8, validate={'validator': 'integer'})]
                self.loops_list.append(field)
                self.loops_list[-1][0].grid(sticky='w', column=1, row=i + 2)
                self.loops_list[-1][1].grid(sticky='w', column=2, row=i + 2)

            self.btns_extend_clear = Pmw.ButtonBox(self.chosen_own_lc, orient="horizontal")
            self.btns_extend_clear.add("Extend", width=self.extend_clear_btns_width,
                                       command=lambda: self.extend_loop_list())
            self.btns_extend_clear.add("Clear", width=self.extend_clear_btns_width,
                                       command=lambda: self.clear_loop_list())
            self.btns_extend_clear.grid(sticky='swn', column=1, row=6, columnspan=2, padx=5, pady=2)

            tmp_frame = tk.Frame(self.chosen_own_lc, padx=0, pady=0)
            self.btns_get_data = Pmw.ButtonBox(tmp_frame, orient="vertical", pady=self.retrieve_btns_pad)
            self.btns_get_data.add("Get from\nstructure", width=self.retrieve_btns_width, font=self.bold_font,
                                   command=lambda: self.get_loop_from_protein_structure())
            self.btns_get_data.add("Get from\nsequence", width=self.retrieve_btns_width, font=self.bold_font,
                                   command=lambda: self.get_loop_from_protein_sequence())
            self.btns_get_data.grid(sticky='wen', column=0, row=0, rowspan=3)
            tmp_frame.grid(column=0, row=2, rowspan=4)
            self.chosen_own_lc.grid(sticky='swen', column=0, row=2)

            self.btns_view = []
            for i in range(3):
                self.btns_view.append(tk.Button(self.chosen_own_lc, text="View", width=self.view_btn_width,
                                                command=lambda x=self.loops_list[i], y=i: self.view_possible_loop(x,y)))
                self.btns_view[-1].grid(sticky='w', column=3, row=i + 2, padx=2)

            self.create_bridge_list_hints()

    def extend_loop_list(self):
        if (0 == self.loops_list[-1][1].getvalue()) or (0 == len(self.loops_list[-1][0].getvalue())):
            self.raise_popup_menu('All fields must be filled in order to load a new loop.')

        self.number_of_own_loops += 1
        self.loops_list.append(
            [Pmw.EntryField(self.chosen_own_lc, labelpos='w', entry_width=8, validate={'validator': 'integer'}),
             Pmw.EntryField(self.chosen_own_lc, labelpos='w', entry_width=8, validate={'validator': 'integer'})])
        self.loops_list[-1][0].grid(sticky='w', column=1, row=len(self.loops_list) + 1, pady=5)
        self.loops_list[-1][1].grid(sticky='w', column=2, row=len(self.loops_list) + 1, pady=5)

        self.btns_view.append(tk.Button(self.chosen_own_lc, text="View", width=self.view_btn_width,
                                        command=lambda x=self.loops_list[-1], y=len(self.loops_list):
                                        self.view_possible_loop(x, y)))
        self.btns_view[-1].grid(sticky='w', column=3, row=len(self.loops_list) + 1, padx=2)
        self.btns_extend_clear.grid_configure(row=len(self.loops_list) + 2)

    def clear_loop_list(self):
        for i in range(0, self.number_of_own_loops):
            self.loops_list[i][0].setentry("")
            self.loops_list[i][1].setentry("")

        for i in range(3, self.number_of_own_loops):
            self.loops_list[i][0].grid_forget()
            self.loops_list[i][1].grid_forget()
            self.btns_view[i].grid_forget()
            self.number_of_own_loops -= 1

        for i in list(self.distance_error.keys()):
            self.distance_error[i].grid_forget()

        self.loops_list[3:] = []
        self.btns_view[3:] = []
        self.distance_error = {}

        self.label_ll.configure(text=textwrap.fill(self.label_ll.cget("text"), 50))
        self.label_ll.grid_configure(columnspan=4)

        if len(cmd.get_names(type="selections")) > 0:
            cmd.delete(name="sele")

    def get_loop_from_protein_structure(self):
        chain = self.chain_index.get()
        atoms = {'atoms': []}

        if not cmd.get_names(type="selections").__contains__("sele"):
            self.raise_popup_menu('No atoms has been chosen.')

        cmd.iterate_state(state=0, selection="sele and chain " + chain + " and name ca",
                          expression='atoms.append(resi)', space=atoms)
        atoms = {'atoms': list(set(atoms['atoms']))}
        atoms['atoms'].sort(key=int)

        if len(atoms['atoms']) > 2:
            self.raise_popup_menu('Too many selection (' + str(len(atoms['atoms'])) + ').')
        elif len(atoms['atoms']) < 2:
            self.raise_popup_menu('Too little selections or there are selected atoms from wrong chain.')

        i = 0
        while i != self.number_of_own_loops and (len(self.loops_list[i][0].getvalue()) is not 0
                                                 and len(self.loops_list[i][0].getvalue()) is not 0):
            i += 1
        if i == self.number_of_own_loops:
            self.extend_loop_list()
        self.loops_list[i][0].setentry(atoms['atoms'][0])
        self.loops_list[i][1].setentry(atoms['atoms'][1])

        if len(cmd.get_names(type="selections")) > 0:
            cmd.delete(name="sele")

    def get_loop_from_protein_sequence(self):
        chain = self.chain_index.get()
        atoms = {'atoms': []}

        if not cmd.get_names(type="selections").__contains__("sele"):
            self.raise_popup_menu('No atoms has been chosen.')

        cmd.iterate_state(state=0, selection="sele and chain " + chain + " and name ca",
                          expression='atoms.append(resi)', space=atoms)
        atoms = {'atoms': list(set(atoms['atoms']))}
        atoms['atoms'].sort(key=int)

        if len(atoms['atoms']) < 2:
            self.raise_popup_menu('Too little selections or there are selected atoms from wrong chain.')

        i = 0
        while (i != self.number_of_own_loops) and (len(self.loops_list[i][0].getvalue()) is not 0
                                                   and len(self.loops_list[i][0].getvalue()) is not 0):
            i += 1
        if i == self.number_of_own_loops:
            self.extend_loop_list()
        self.loops_list[i][0].setentry(atoms['atoms'][0])
        self.loops_list[i][1].setentry(atoms['atoms'][-1])

        if len(cmd.get_names(type="selections")) > 0:
            cmd.delete(name="sele")

    def view_possible_loop(self, loop, pos):
        if (len(loop[0].getvalue()) is 0 and len(loop[1].getvalue()) is not 0) or (
                        len(loop[0].getvalue()) is not 0 and len(loop[1].getvalue()) is 0):
            self.raise_popup_menu('There are missing data in fields.')
        if len(loop[0].getvalue()) is 0 or len(loop[1].getvalue()) is 0:
            self.raise_popup_menu('No data in fields. There is nothing to show.')

        atom = "ca"
        chain = self.chain_index.get()
        atoms = {'atoms': []}
        cmd.iterate_state(state=0, selection="all and chain " + chain, expression="atoms.append(resi)", space=atoms)
        atoms = atoms['atoms']
        if (int(loop[0].getvalue()) < int(atoms[0])) or (int(loop[1].getvalue()) > int(atoms[-1])):
            self.raise_popup_menu('There are no such indices in the structure. Please choose other (correct) ones.')

        for i in cmd.get_names(type="all"):
            if str(i).startswith("TMP") or str(i).startswith("DIST_") or str(i).startswith("CHAIN_") or \
                    str(i).startswith("SEQ") or str(i).startswith("PIERC") or str(i).startswith("TRIANG") or \
                    str(i).startswith("BR*") or str(i).startswith("SMOOTH_CHAIN_"):
                cmd.hide(representation="everything", selection=i)
                cmd.delete(name=i)

        if self.previous_bond_in_view[0] is not "" and self.previous_bond_in_view[1] is not "":
            cmd.unbond(atom1=self.previous_bond_in_view[0], atom2=self.previous_bond_in_view[1])

        self.simplify_polymer_representation()

        cmd.set(name="seq_view", value=1)
        cmd.set(name="mouse_selection_mode", value=6)
        cmd.select(name="CHAIN_" + chain, selection="all and chain " + chain)
        cmd.hide(representation="sticks", selection="all")
        cmd.show(representation="cartoon", selection="CHAIN_" + chain)
        cmd.cartoon(type="tube", selection="CHAIN_" + chain)
        cmd.spectrum(palette="rainbow", selection="all")
        cmd.deselect()

        cmd.select("TMP_BR_" + loop[0].getvalue() + "_" + loop[1].getvalue(),
                   selection="(chain " + chain + " and residue " + loop[0].getvalue() + " and name " + atom + ")+"
                             "(chain " + chain + " and residue " + loop[1].getvalue() + " and name " + atom + ")")
        cmd.select("TMP_SEQ",
                   selection="chain " + chain + " and residue " + loop[0].getvalue() + "-" + loop[1].getvalue())
        cmd.bond(atom1="chain " + chain + " and residue " + loop[0].getvalue() + " and name " + atom,
                 atom2="chain " + chain + " and residue " + loop[1].getvalue() + " and name " + atom)
        cmd.color(color="gray", selection="TMP_SEQ")
        cmd.hide(representation="lines", selection="TMP_SEQ and sc.")
        cmd.show(representation='sphere', selection="TMP_BR_" + loop[0].getvalue() + "_" + loop[1].getvalue())
        cmd.show(representation="sticks", selection="TMP_BR_" + loop[0].getvalue() + "_" + loop[1].getvalue())
        cmd.set(name="sphere_color", value="orange", selection="all")
        cmd.set(name="stick_color", value="orange", selection="all")
        cmd.set("sphere_scale", value=0.5)

        self.check_distance_between_atoms(loop, pos)
        self.previous_bond_in_view = ["chain " + chain + " and residue " + loop[0].getvalue() + " and name " + atom,
                                      "chain " + chain + " and residue " + loop[1].getvalue() + " and name " + atom]

    def check_distance_between_atoms(self, loop, pos):
        atom = "ca"

        if self.is_trajectory:
            chain = self.chains[0]
            selection1 = "(chain " + chain + " and residue " + loop[0].getvalue() + " and name " + atom + ")"
            selection2 = "(chain " + chain + " and residue " + loop[1].getvalue() + " and name " + atom + ")"
        else:
            chain = self.chain_index.get()
            selection1 = "(chain " + chain + " and residue " + loop[0].getvalue() + " and name " + atom + ")"
            selection2 = "(chain " + chain + " and residue " + loop[1].getvalue() + " and name " + atom + ")"

        cmd.distance(name="DIST_" + loop[0].getvalue() + "_" + loop[1].getvalue(),
                     selection1=selection1, selection2=selection2)

        calc_distance = round(cmd.get_distance(atom1=selection1, atom2=selection2), 1)
        print("  ## DISTANCE between atoms " + loop[0].getvalue() + " and " + loop[1].getvalue() + " is " + \
              str(calc_distance))

        usr_dist = (self.min_dist_crossings.getvalue() if self.is_stable.get() == 0 else 10)
        if pos < 4:
            pos += 1
        self.create_view_loop_hints(calc_distance, usr_dist, pos)

    def create_view_loop_hints(self, calc_dist, usr_dist, pos):
        self.list_hint_distance = []

        interior = self.chosen_own_lc if hasattr(self, "chosen_own_lc") else self.fr_selected_loop.interior()

        if float(calc_dist) > float(usr_dist):
            if pos not in list(self.distance_error.keys()):
                if hasattr(self, "chosen_own_lc"):
                    self.label_ll.configure(text=textwrap.fill(self.label_ll.cget("text"), 55))
                    self.label_ll.grid_configure(columnspan=5)
                else:
                    self.label_trajectory_loop.configure(
                        text=textwrap.fill(self.label_trajectory_loop.cget("text"), 60))
                    self.label_trajectory_loop.grid_configure(columnspan=5)
                self.distance_error[pos] = (tk.Label(interior, text='too big!', foreground="red"))
                self.distance_error[pos].grid(column=4, row=pos + 1)
                self.list_hint_distance.append(Pmw.Balloon(interior, relmouse="both"))
                self.list_hint_distance[-1].bind(self.distance_error[pos], textwrap.fill("A distance between residues "
                                                "is too big to form a cysteine bridge. Another pair of residues "
                                                "forming a bridge should be chosen or an option 'Ignore bad length of "
                                                "bridge or Ca-Ca bonds' should be ticked.", self.hint_width))
        else:
            if pos in list(self.distance_error.keys()):
                self.distance_error[pos].grid_forget()
                del self.distance_error[pos]

    def create_bridge_list_hints(self):
        hint_choose_bridge_type = Pmw.Balloon(self.chosen_own_lc, relmouse="both")
        hint_choose_bridge_type.bind(self.label_ll,
                                     "Provide sequential numbers of two residues forming a bridge.")
        hint_choose_bridge_type.bind(self.btns_get_data.button(0), textwrap.fill("Retrieve the selected tuple of "
                                                                                 "C-alpha atoms from the (bio)polymer "
                                                                                 "loaded in PyMOL.", self.hint_width))
        hint_choose_bridge_type.bind(self.btns_get_data.button(1), textwrap.fill("Retrieve the selected sequence of "
                                                                                 "C-alpha atoms from the (bio)polymer "
                                                                                 "loaded in PyMOL.", self.hint_width))
        hint_choose_bridge_type.bind(self.btns_extend_clear.button(0), "Add one more loop.")
        hint_choose_bridge_type.bind(self.btns_extend_clear.button(1), "Clear the above data.")

        for i in self.btns_view:
            hint_choose_bridge_type.bind(i, "View the loop closed by the the bridge chosen.")

    ####################################################################################################################
    #                                              3) ADVANCED OPTIONS
    ####################################################################################################################

    def create_smooth_interior(self):
        self.group_smooth = Pmw.Group(self.group_advanced.interior(), tag_text='Smoothing')
        self.group_smooth.grid(sticky='eswn', column=1, row=0, padx=5, pady=5)

        self.smooth_val = Pmw.EntryField(self.group_smooth.interior(), labelpos='w', label_text='Level of smoothness',
                                         validate={'validator': 'integer', 'min': 1, "max": 100},
                                         value=2, entry_width=5)
        self.smooth_val.grid(sticky='w', column=0, row=0, padx=2, pady=2)

    def create_correct_mode_interior(self):
        self.group_correct_mode = Pmw.Group(self.group_advanced.interior(), tag_text='Correct mode')
        self.group_correct_mode.grid(sticky='eswn', column=1, row=1, padx=5, pady=5)

        self.stable_lasso = tk.Checkbutton(self.group_correct_mode.interior(), text='Stable lasso', onvalue=1,
                                           offvalue=0, variable=self.is_stable,
                                           command=lambda: self.enable_parametrization_of_algorithm())
        self.stable_lasso.grid(sticky='w', column=0, row=0, padx=2, pady=2)
        self.stable_lasso.select()

        self.bad_caca = tk.Checkbutton(self.group_correct_mode.interior(), variable=self.is_bad_caca_enabled,
                                       state=tk.DISABLED,
                                       text='Ignore an inappropriate length of a bridge and Ca-Ca bond')
        self.bad_caca.grid(sticky='w', column=0, row=1, padx=2, pady=2)

        self.min_dist_crossings = Pmw.EntryField(self.group_correct_mode.interior(), labelpos='w',
                                                 label_text='Minimal distance between crossings',
                                                 entry_state=tk.DISABLED,
                                                 entry_width=5, validate={'validator': 'integer'}, value='10')
        self.min_dist_crossings.grid(sticky='eswn', column=0, row=2, padx=2, pady=2)

        self.min_dist_cross_end = Pmw.EntryField(self.group_correct_mode.interior(), labelpos='w',
                                                 label_text='Minimal distance between a crossing and a tail end',
                                                 entry_state=tk.DISABLED, entry_width=5,
                                                 validate={'validator': 'integer'},
                                                 value='3')
        self.min_dist_cross_end.grid(sticky='eswn', column=0, row=3, padx=2, pady=2)

        self.min_dist_cross_loop = Pmw.EntryField(self.group_correct_mode.interior(), labelpos='w',
                                                  label_text='Minimal distance between a crossing and a loop',
                                                  entry_state=tk.DISABLED, entry_width=5,
                                                  validate={'validator': 'integer'},
                                                  value='3')
        self.min_dist_cross_loop.grid(sticky='eswn', column=0, row=4, padx=2, pady=2)
        self.stable_lasso.grid(column=0, row=0, padx=2, pady=2)

    def create_gln_interior(self):
        self.group_gln = Pmw.Group(self.group_advanced.interior(), tag_text='Gaussian Linking Number')
        self.group_gln.grid(sticky='eswn', column=1, row=2, padx=5, pady=5)

        self.gln_checkbutton = tk.Checkbutton(self.group_gln.interior(), text='Calculate GLN',
                                              variable=self.is_gln_checkbutton_selected)
        self.gln_checkbutton.grid(sticky='w', column=0, row=0, padx=2, pady=2)

    def enable_parametrization_of_algorithm(self):
        if self.is_stable.get():
            self.min_dist_crossings.configure(entry_state='disabled')
            self.min_dist_cross_loop.configure(entry_state='disabled')
            self.min_dist_cross_end.configure(entry_state='disabled')
            self.bad_caca.configure(state='disabled')
            self.bad_caca.deselect()
        else:
            self.min_dist_crossings.configure(entry_state='normal')
            self.min_dist_cross_loop.configure(entry_state='normal')
            self.min_dist_cross_end.configure(entry_state='normal')
            self.bad_caca.configure(state='normal')
            self.bad_caca.select()

    def create_advanced_frame_hints(self):
        hint_advanced = Pmw.Balloon(self.main_window, relmouse="both")
        if not self.is_trajectory:
            hint_advanced.bind(self.group_gln, textwrap.fill("A numerical invariant that describes the linking of two "
                                                             "closed curves in the three-dimensional space. "
                                                             "Intuitively, the linking number represents the number of "
                                                             "times that each curve winds around the other one.",
                                                             self.hint_width))
        hint_advanced.bind(self.group_smooth, textwrap.fill("Adjust the level of smoothness, by choosing the value of "
                                                            "this parameter in the range 2-100. The higher value "
                                                            "results in more refined triangulation of the minimal "
                                                            "surface. Above a certain value (characteristic for a "
                                                            "given configuration) there is no visible change in the "
                                                            "triangulation.", self.hint_width))
        hint_advanced.bind(self.group_correct_mode,
                           "Impose geometric conditions that define shallow lassos.")
        hint_advanced.bind(self.bad_caca,
                           "Ignores the validity check.")

    def create_protein_interface_hints(self):
        hint_protein = Pmw.Balloon(self.main_window, relmouse="both")
        hint_protein.bind(self.fr_chain, textwrap.fill("Choose a chain to be analyzed. If the field is empty, the "
                                                       "first chain (typically A in the case of proteins) will be "
                                                       "analyzed. For files in the XYZ format only one chain is "
                                                       "available.", self.hint_width))
        hint_protein.bind(self.type_loop_closing_bridge, textwrap.fill("Choose a method to identify closed loops. For "
                                                                       "files in XYZ format only the last method is "
                                                                       "available.", self.hint_width))

    ####################################################################################################################
    #                                               EXECUTE PROGRAM 1/2
    ####################################################################################################################

    def _invoke_program(self):
        self.displayed_lasso = None

        self.convert_to_5columns_format()

        self.delete_pymol_objects()
        cmd.spectrum(palette="rainbow", selection="all")

        if self._file_extension == "pdb":
            self.warning_gaps = [(i.split(" ")[3], i.split(" ")[-3], i.split(" ")[-1]) for i in self.pdb_bridges if
                                 i.__contains__("WARNING")]
            self.pdb_bridges = [num for num in self.pdb_bridges if not num.__contains__("WARNING")]

        self.user_data = self.generate_invoking_commands()

        if not self.is_trajectory:
            self.display_pymol_chain()

        self.call_lasso_detection()

        self.separate_smooth_crossings_from_output()

        if hasattr(self, "type_loop_closing_bridge") \
                and not str(self.type_loop_closing_bridge.get()) == 'choose two atoms to form a bridge' \
                and all(i.__contains__("ERROR") for i in self.output_data) and not self.is_trajectory:
            if self.is_stable.get() and not self.is_bad_caca_enabled.get():
                if hasattr(self, "artifact_found") and self.artifact_found.winfo_exists():
                    self.artifact_found.withdraw()
                greatest_gap = self.get_greatest_gap()
                artifact_message = "Detected lasso(s) may be artificial. The biggest gap in chain is " + greatest_gap +\
                                   " amino acids. An option ''Ignore an inappropriate length of a bridge and Ca-Ca " \
                                   "bond'' has been ticked in order to identify bridge(s) and to determine unbroken " \
                                   "segment of a chain (a protein backbone). Be careful with interpreting results!"
                self.artifact_found = Pmw.MessageDialog(self.parent, title=' ', defaultbutton=0,
                                                        message_text=textwrap.fill(artifact_message, self.hint_width))
                self.artifact_found.geometry("+%d+%d" % (self.screen_width / 2 - 150, self.screen_height / 2))
                self.call_program_for_artifacts()

        if self.is_gln_checkbutton_selected.get():
            print("  Generating GLN matrices...")
            self.call_gln_generator()

        self.move_files_to_polymer_directory()

        if hasattr(self, "win_lasso_info") and self.win_lasso_info.winfo_exists():
            self.win_lasso_info.destroy()
        if hasattr(self, "win_trajectory_analysis") and self.win_trajectory_analysis.winfo_exists():
            self.win_trajectory_analysis.destroy()

        self.create_trajectory_window() if self.is_trajectory else self.create_protein_window()

        if os.path.exists("niewaznypliczek.txt"):  # a file needed for proper executing program finding lasso, not
            os.remove("niewaznypliczek.txt")       # needed in further calculations

        if hasattr(self, "error_pop_menu") and self.error_popup.winfo_exists():
            self.error_popup.withdraw()

    def convert_to_5columns_format(self):
        self.pdb_bridges = None
        all_bridges = None

        if self.is_trajectory:
            all_bridges = subprocess.Popen(self.python_compiler + (self._full_path_to_file + " -f -t").split(" "),
                                           stdout=subprocess.PIPE).communicate()[0]
        else:
            all_bridges = subprocess.Popen(self.python_compiler + [self._full_path_to_file],
                                           stdout=subprocess.PIPE).communicate()[0].splitlines()
        if self.is_trajectory:
            self.pdb_bridges = all_bridges
        else:
            chain = self.chain_index.get()
            file_with_bridge = (self._full_path_to_file + "_" + chain)

            all_bridges = list(filter(len, all_bridges))
            self.pdb_bridges = [i for i in all_bridges if i.__contains__(file_with_bridge) or
                                                i.__contains__("WARNING")]

    ####################################################################################################################
    #                                               1) GET PLUGIN DATA
    ####################################################################################################################

    def generate_invoking_commands(self):
        print("  Collecting user's data...")
        arguments = []
        self.list_bridges = []

        arguments += self.get_trajectory_data() if self.is_trajectory else self.get_method_to_find_bridges()

        if self.is_gln_checkbutton_selected.get():
            arguments = [i + " -f" for i in arguments]
            arguments = [i + " 3" for i in arguments]
        else:
            arguments = [i + " -f" for i in arguments]
            arguments = [i + " 2" for i in arguments]
        if len(self.smooth_val.getvalue()) < 0:
            self.raise_popup_menu('No surface value given. Please insert appropriate value from range <2;100>.')
        if len(self.smooth_val.getvalue()) > 0:
            arguments = [i + " -sm_nr" for i in arguments]
            arguments = [i + " " + str(self.smooth_val.getvalue()) for i in arguments]
        if self.is_bad_caca_enabled.get():
            arguments = [i + " -cd" for i in arguments]
            arguments = [i + " 0" for i in arguments]
        if not self.is_stable.get():
            arguments = [i + " -redAC" for i in arguments]
            arguments = [i + " " + self.min_dist_crossings.getvalue() for i in arguments]
            arguments = [i + " -redEnd" for i in arguments]
            arguments = [i + " " + self.min_dist_cross_end.getvalue() for i in arguments]
            arguments = [i + " -redBr" for i in arguments]
            arguments = [i + " " + self.min_dist_cross_loop.getvalue() for i in arguments]
        return arguments

    def get_trajectory_data(self):
        args = self.program_execution
        new_filename = self._full_path_to_file + "_" + self.chains[0] + ".xyz"
        args += new_filename + " "
        args += self.get_trajectory_loop()
        args += self.get_accuracy_of_calculations()
        args += self.get_step()
        return [args]

    def get_trajectory_loop(self):
        self.validate_trajectory()
        return str(self.loops_list[0][0].getvalue()) + " " + str(self.loops_list[0][1].getvalue()) + " "

    def get_accuracy_of_calculations(self):
        args = ""

        args += "-traj " + ("2 " if self.is_detailed_alg.get() else "1 ")
        args += "-trajout " + ("1 " if self.is_detailed_out.get() else "0 ")
        return args

    def get_step(self):
        if len(self.step.getvalue()) == 0:
            return "-step 1"
        return "-step " + str(self.step.getvalue())

    def get_method_to_find_bridges(self):
        args = []

        if str(self.type_loop_closing_bridge.get()) == "automatic detections of closed loops":
            args += self.get_automatic_closing_data()
        elif str(self.type_loop_closing_bridge.get()) == "choose a type of loop closing":
            args += self.get_type_closing_data()
        elif str(self.type_loop_closing_bridge.get()) == "choose two atoms to form a bridge":
            args += self.get_own_closing_data()
        return args

    def get_automatic_closing_data(self):
        if self._file_extension == "xyz":
            self.raise_popup_menu('No automatic detection of closed loops for .xyz files.')
        if not self.is_original_pdb:
            self.raise_popup_menu('No bridge has been detected in the file. Please use na option "choose two atoms '
                                  'to form a bridge"')

        chain = self.chain_index.get()
        list_atoms = []
        list_args = []

        for i in list(self.pdb_bridges):
            elem = i.split(" ")
            if len(elem) > 2:
                if int(elem[-2]) > int(elem[-1]):
                    bridge = (elem[0], elem[-1], elem[-2])
                else:
                    bridge = (elem[0], elem[-2], elem[-1])
                list_atoms.append(bridge)
        for i in list_atoms:
            if not list_args.__contains__(self.program_execution + self._full_path_to_file + "_" + chain + ".xyz " +
                                                  i[1] + " " + i[2]):
                list_args.append(self.program_execution + self._full_path_to_file + "_" + chain + ".xyz " + i[1] +
                                 " " + i[2])
                self.list_bridges.append(i)
        return list_args

    def get_type_closing_data(self):
        if self._file_extension == "xyz":
            self.raise_popup_menu('Not possible to choose a type of loop closing for .xyz files.')
        if not self.is_original_pdb:
            self.raise_popup_menu('No bridge has been detected in the file. Please use na option "choose two atoms '
                                  'to form a bridge"')
        if hasattr(self, "bridge_selected") and None is self.bridge_selected:
            self.raise_popup_menu('No bridge type selected. Please select correct one.')

        chain = self.chain_index.get()
        list_atoms = []
        list_atoms_br = []
        list_args = []

        for i in list(self.pdb_bridges):
            elem = i.split(" ")

            if len(elem) > 2:
                if int(elem[-2]) > int(elem[-1]):
                    bridge = (elem[0], elem[-1], elem[-2])
                else:
                    bridge = (elem[0], elem[-2], elem[-1])
                if bridge[0] == (self.bridge_selected.upper()) or bridge[0] == (self.bridge_selected.upper() + "-like"):
                    list_atoms_br.append(bridge)
                else:
                    list_atoms.append(bridge)
        if len(list_atoms_br) > 0:  # A bridge - specified in bridge-type-selection - exists
            self.type_bridge_exists = True
            for i in list_atoms_br:
                if not list_args.__contains__(self.program_execution + self._full_path_to_file + "_" + chain + ".xyz " +
                                                      i[1] + " " + i[2]):
                    list_args.append(
                        self.program_execution + self._full_path_to_file + "_" + chain + ".xyz " + i[1] + " " + i[2])
                    self.list_bridges.append(i)
        else:
            if hasattr(self, "bridge_found") and self.bridge_found.winfo_exists():
                self.bridge_found.withdraw()
            self.bridge_found = Pmw.MessageDialog(self.parent, title=' ', defaultbutton=0,
                                                  message_text=textwrap.fill("Selected bridge was not found! Further "
                                                                             "calculations were conducted "
                                                                             "automatically.", self.hint_width))
            self.bridge_found.geometry("+%d+%d" % (self.screen_width / 2 - 150, self.screen_height / 2))
            self.type_bridge_exists = False

            for i in list_atoms:
                if not list_args.__contains__(self.program_execution + self._full_path_to_file + "_" + chain + ".xyz " +
                                                      i[1] + " " + i[2]):
                    list_args.append(
                        self.program_execution + self._full_path_to_file + "_" + chain + ".xyz " + i[1] + " " + i[2])
                    self.list_bridges.append(i)
        return list_args

    def get_own_closing_data(self):
        self.validate_loop_list()
        chain = self.chain_index.get()
        list_args = []

        for i in self.loops_list:
            if (len(i[0].getvalue()) == 0) or len(i[1].getvalue()) == 0:
                break
            if not list_args.__contains__(self.program_execution + self._full_path_to_file + "_" + chain + ".xyz " + str(
                    i[0].getvalue()) + " " + str(i[1].getvalue())):
                list_args.append(
                    self.program_execution + self._full_path_to_file + "_" + chain + ".xyz " + str(i[0].getvalue()) +
                    " " + str(i[1].getvalue()))
        return list_args

    def validate_loop_list(self):
        if len(self.loops_list[0][0].getvalue()) is 0 and len(self.loops_list[0][1].getvalue()) is 0:
            self.raise_popup_menu('At least one loop must be given in above fields.')

        for elem in self.loops_list:
            if (len(elem[0].getvalue()) is 0 and len(elem[1].getvalue()) is not 0) \
                    or (len(elem[0].getvalue()) is not 0 and len(elem[1].getvalue()) is 0):
                self.raise_popup_menu('There are missing data in fields.')

    ####################################################################################################################
    #                                               EXECUTE PROGRAM 2/2
    ####################################################################################################################

    def display_pymol_chain(self):
        if self._file_extension == "xyz":
            cmd.load(filename=self._full_path_to_file)
            self.connect_xyz_points()
            cmd.show(representation="lines", selection=self._filename[:-9])
            cmd.set(name="line_width", value="4")
        else:
            chain = self.chain_index.get()
            cmd.hide(representation="everything", selection="all")
            cmd.select(name="CHAIN_" + chain, selection="all and chain " + chain)
            cmd.show(representation="cartoon", selection="CHAIN_" + chain)
            cmd.cartoon(type="tube", selection="all")
        cmd.spectrum(palette="rainbow", selection="all")
        cmd.deselect()

    def call_lasso_detection(self):
        self.output_data = []
        try:
            for i in self.user_data:
                process = subprocess.Popen(i.split(" "), stdout=subprocess.PIPE).communicate()[0]
                self.output_data.append(process)
            self.output_data = list(filter(len, self.output_data))
        except Exception:
            print("Something went wrong with executable file. Please make sure you changed access permission to " \
                  "it (can be obtained by typing in console chmod a+x detect_lassos).")
        print("  Data passed to program and executed...")

    def call_gln_generator(self):
        iterate_list = self.user_data
        chain = self.chain_index.get()

        for i in iterate_list:
            command = i.split(" ")
            res = [command[2], command[3]]
            matrix1 = self.python_compiler[0] + " matrixGLNtoPNG.py matrixGLN_" + self._filename + "_" + chain + "_" + \
                      res[0] + "_" + res[1] + "_t1"
            matrix2 = self.python_compiler[0] + " matrixGLNtoPNG.py matrixGLN_" + self._filename + "_" + chain + \
                      "_" + res[0] + "_" + res[1] + "_t2"
            subprocess.Popen(matrix1.split(" "), stdout=subprocess.PIPE).communicate()[0]
            subprocess.Popen(matrix2.split(" "), stdout=subprocess.PIPE).communicate()[0]

    def separate_smooth_crossings_from_output(self):
        self.smooth_crossings = []

        for i in self.output_data:
            elem = i.split("\n")
            elem = list(filter(len, elem))
            smooth_cross = elem[-1]
            if smooth_cross.__contains__("ERROR"):
                self.smooth_crossings.append([])
            else:
                n_cross = self.get_n_crossings(smooth_cross.split(" "))
                c_cross = self.get_c_crossings(smooth_cross.split(" "))
                self.smooth_crossings.append((n_cross + c_cross).split(" "))

    def get_greatest_gap(self):
        chain = self.chain_index.get()
        max = 0
        for i in self.warning_gaps:
            if i[0] == chain and int(int(i[2]) - int(i[1])) > max:
                max = int(int(i[2]) - int(i[1]))
        return str(max)

    def call_program_for_artifacts(self):
        self.output_data = []
        self.is_artifact = True
        self.stable_lasso.deselect()
        self.enable_parametrization_of_algorithm()
        self.user_data = self.generate_invoking_commands()

        for i in self.user_data:
            self.output_data.append(subprocess.Popen(i.split(" "), stdout=subprocess.PIPE).communicate()[0])
        self.output_data = list(filter(len, self.output_data))
        print("  Modified data passed to program again and executed...")

    def move_files_to_polymer_directory(self):
        self.current_working_dir = system_working_directory
        directory = self.create_polymer_directory(self._filename.replace(".", "_"))
        directory_in_workspace = self._full_path_to_dir + os.sep + directory

        os.chdir(directory_in_workspace)

        if self.is_gln_checkbutton_selected.get():
            self.separate_files_to_directory(self.current_working_dir, os.getcwd() + os.sep + "_GLN", "matrixGLN_")
        self.separate_files_to_directory(self.current_working_dir, os.getcwd() + os.sep + "_barycentric",
                                         "barycentric_", "F_PYsvgBari_")
        self.separate_files_to_directory(self.current_working_dir, os.getcwd() + os.sep + "_surfaces", "surface_")
        self.separate_files_to_directory(self.current_working_dir, os.getcwd() + os.sep + "_smooth", "_smooth.pdb")
        self.separate_files_to_directory(self.current_working_dir,  os.getcwd() + os.sep, self._filename + "_")
        self.separate_files_to_directory(self._full_path_to_dir, os.getcwd() + os.sep, self._filename + "_")
        os.chdir(self.current_working_dir)
        print("  Resulting files moved to separate directories...")

    def create_polymer_directory(self, prot):
        direct = os.sep.join(self._full_path_to_file.split(os.sep)[:-1]) + os.sep + prot

        if not os.path.exists(direct):
            os.makedirs(direct)
            print("  Creating a new, separate directory for data...")
        else:
            print("  Directory found. Replacing existing files with new data...")
        return direct.split(os.sep)[-1]

    def separate_files_to_directory(self, source, target, *files):
        if not os.path.exists(target):
            os.makedirs(target)

        for filename in os.listdir(source):
            for arg in files:
                if filename.__contains__(arg):
                    shutil.move(os.path.join(source, filename),
                                os.path.join(target, filename))

    ####################################################################################################################
    #                                       1. CREATE TRAJECTORY RESULT WINDOW
    ####################################################################################################################

    def create_trajectory_window(self):
        if hasattr(self, "win_trajectory_") and self.win_trajectory_.winfo_exists():
            self.win_trajectory_.withdraw()

        self.win_trajectory_ = Pmw.Dialog(self.parent, title='Trajectory information window [' + self._filename + ']',
                                          buttons=["Show", "Exit"],
                                          command=self.invoke_trajectory_window_buttons)
        self.win_trajectory_.resizable(0, 0)

        self.win_trajectory_analysis = Pmw.ScrolledFrame(self.win_trajectory_.interior(), borderframe=0,
                                                         hscrollmode="none", usehullsize=1, hull_width=self.hull_width,
                                                         hull_height=600)
        self.win_trajectory_analysis.grid(sticky="swen", column=0, row=0)
        self.win_trajectory_.geometry("+%d+%d" % (self.screen_width / 2, self.screen_height / 4))

        self.win_trajectory_lasso_information = Pmw.Group(self.win_trajectory_analysis.interior(),
                                                          tag_text="Lasso information")
        self.win_trajectory_lasso_information.grid(sticky="swen", column=0, row=0, columnspan=2, padx=5, pady=5)

        self.win_lasso_type = Pmw.Group(self.win_trajectory_analysis.interior(), tag_text="Lasso(s) type")
        self.win_lasso_type.grid(sticky="swen", column=0, row=1, columnspan=2, padx=5, pady=5)

        self.win_atoms_piercing_lasso = Pmw.Group(self.win_trajectory_analysis.interior(),
                                                  tag_text="Atoms piercing the closed loop")
        self.win_atoms_piercing_lasso.grid(sticky="swen", column=0, row=2, columnspan=2, padx=5, pady=5)

        self.create_win_frame_interior()

        self.win_trajectory_analysis_log = Pmw.Group(self.win_trajectory_analysis.interior(),
                                                     tag_text="Trajectory analysis log")
        self.win_trajectory_analysis_log.grid(sticky="swen", column=1, row=3, padx=5, pady=5)

        self.delete_show_button()
        self.create_trajectory_window_hints()

    def create_win_frame_interior(self):
        self.win_frame_number = Pmw.Group(self.win_trajectory_analysis.interior(), tag_text='Frame number')
        self.win_frame_number.grid(sticky="swen", column=0, row=3, padx=5, pady=5)

        self.frame_list = []
        self.frame_list.append(Pmw.EntryField(self.win_frame_number.interior(), labelpos='w', entry_width=7,
                                              validate={'validator': 'real'}))
        self.frame_list[-1].grid(sticky='eswn', column=0, row=0, padx=5, pady=5)

        self.frame_button = Pmw.ButtonBox(self.win_frame_number.interior(), orient="vertical")
        self.frame_button.add("Calculate", command=lambda: self.calculate_trajectory_frame())
        self.frame_button.grid(sticky='e', column=2, columnspan=4, row=0, padx=5)

        self.more_detailed_frame_algorithm = tk.Checkbutton(self.win_frame_number.interior(),
                                                            text='More accurate algorithm - \na single frame approach',
                                                            variable=self.is_detailed_out_frame)
        self.more_detailed_frame_algorithm.grid(sticky='ne', column=2, columnspan=4, row=1, rowspan=5, padx=2, pady=2)

        self.frame_btns = Pmw.ButtonBox(self.win_frame_number.interior(), orient="horizontal")
        self.frame_btns.add("Extend", width=self.extend_clear_btns_width, command=lambda: self.extend_frame_list())
        self.frame_btns.add("Clear", width=self.extend_clear_btns_width, command=lambda: self.clear_frame_list())
        self.frame_btns.grid(sticky='swn', column=0, row=2, columnspan=2, padx=5, pady=2)

    def extend_frame_list(self):
        if len(self.frame_list) > 5:
            self.raise_popup_menu('Only six frames are possible to calculate.')
        if 0 == len(self.frame_list[-1].getvalue()):
            self.raise_popup_menu('All fields must be filled in order to load a frame field.')

        self.frame_list.append(Pmw.EntryField(self.win_frame_number.interior(), labelpos='w', entry_width=7,
                                              validate={'validator': 'real'}))
        self.frame_list[-1].grid(sticky='eswn', column=0, row=len(self.frame_list) + 1, padx=5, pady=5)

        self.frame_btns.grid_configure(row=len(self.frame_list) + 2)

    def clear_frame_list(self):
        for i in self.frame_list:
            i.grid_forget()

        self.frame_list = []
        self.frame_list.append(Pmw.EntryField(self.win_frame_number.interior(), labelpos='w', entry_width=7,
                                              validate={'validator': 'real'}))
        self.frame_list[-1].grid(sticky='eswn', column=0, row=0, padx=5, pady=5)
        self.frame_btns.grid_configure(row=2)

        if len(cmd.get_names(type="selections")) > 0:
            cmd.delete(name="sele")

    def delete_show_button(self):
        self.win_trajectory_.invoke("Show")
        self.win_trajectory_.component("buttonbox").delete(0)

    def create_trajectory_window_hints(self):
        hint_trajectory = Pmw.Balloon(self.win_trajectory_analysis.interior(), relmouse="both")
        hint_trajectory.bind(self.win_lasso_type, textwrap.fill("A chart presenting detected lasso types (the "
                                                                "vertical axis) at a given time frame (the horizontal "
                                                                "axis)", self.hint_width))
        hint_trajectory.bind(self.win_atoms_piercing_lasso, textwrap.fill("A chart presenting sequential numbers of "
                                                                          "piercings. The orange rectangle represents "
                                                                          "atoms which form a closed loop.",
                                                                          self.hint_width))
        hint_trajectory.bind(self.win_frame_number.interior(), textwrap.fill("Use the single frame version of the "
                                                                             "algorithm instead of the simplified "
                                                                             "algorithm for trajectories (the results "
                                                                             "may vary).", self.hint_width))
        if self.is_detailed_out.get():
            hint_trajectory.bind(self.traj_log.interior(), textwrap.fill(
                "Each line of the output consists of the frame number "
                "followed by: the number of piercings for tail N and "
                "tail C before reduction and information if surface is"
                "orientable X indices of all piercings (reduced ones "
                "as well) in tail N X indices of all piercings (reduced"
                " ones as well) in tail C XX the number of piercings for"
                " tail N and tail C after reduction X indices of all "
                "piercings (only not reduced ones) in tail N X indices "
                "of all piercings (only not reduced ones) in tail C XX "
                "the lasso type XX the radius of gyration. ",
                self.hint_width))
        else:
            hint_trajectory.bind(self.traj_log.interior(), textwrap.fill(
                "Each line of the output consists of the frame number "
                "followed by: the number of piercings for N-terminus "
                "tail and tail C-terminal tail| indices of piercings"
                " in tail N | indices of piercings in tail C | the "
                "lasso type.", self.hint_width))

    ####################################################################################################################
    #                                     FILL TRAJECTORY WINDOW WITH DATA
    ####################################################################################################################

    def invoke_trajectory_window_buttons(self, clicked_button):
        if clicked_button == "Show":
            self.get_chart_data_from_file()
            self.set_trajectory_analysis_log()
            self.calculate_lasso_in_trajectory()

            self.window_parent = self.win_trajectory_lasso_information.interior()

            self.display_lasso_information_table()
            if not all(elem.__contains__("ERROR") for elem in self.output_data):
                self.create_lasso_information_buttons()
                self.create_surface_hints()

            reversed_trajectory_lasso_set = list(set(self.retrieved_trajectory_lassos))
            if len(self.retrieved_frames) == 0:
                self.draw_error_charts("Given step is bigger than total number of frames. There is nothing to draw.",
                                       self.win_lasso_type.interior(), self.win_atoms_piercing_lasso.interior())
            elif all(lasso.__contains__("ERR") for lasso in reversed_trajectory_lasso_set):
                self.draw_error_charts("Please use ''Ignore an inappropriate length of a bridge and Ca-Ca bonds'' "
                                       "option to ignore\ninappropriate length between consecutive  atoms or the "
                                       "loop closing atoms - the bridge.",
                                       self.win_lasso_type.interior(), self.win_atoms_piercing_lasso.interior())
            else:
                self.draw_lassos_type_chart()
                self.draw_atoms_piercing_lasso_chart()
        else:
            self.given_frames = []

            if hasattr(self, "lasinf_smooth_button") and self.lasinf_smooth_button.winfo_exists():
                self.lasinf_smooth_button.grid_forget()
            if hasattr(self, "lasinf_surface_button") and self.lasinf_surface_button.winfo_exists():
                self.lasinf_surface_button.grid_forget()
            if hasattr(self, "lasinf_shallow_lasso_button") and self.lasinf_shallow_lasso_button.winfo_exists():
                self.lasinf_shallow_lasso_button.grid_forget()
            if hasattr(self, "win_trajectory_"):
                self.win_trajectory_analysis.grid_forget()
                self.win_trajectory_.withdraw()
            else:
                self.win_trajectory_analysis.withdraw()

    def get_chart_data_from_file(self):
        path_to_trajectory = self._full_path_to_dir + os.sep + self._filename.replace(".", "_")

        if not os.path.exists(path_to_trajectory):
            self.raise_popup_menu("A path to generated data of trajectory does not exist. There is no such directory.")

        self.file_with_trajectory = path_to_trajectory + os.sep + "traj_" + self._filename + "_" + self.chains[0] + \
                                    "_" + str(self.loops_list[0][0].getvalue()) + "_" + \
                                    self.loops_list[0][1].getvalue() + ".txt"

        self.trajectory_chain_range = []
        self.trajectory_chain_loop_indexes = []
        self.retrieved_frames = []
        self.retrieved_trajectory_lassos = []
        self.retrieved_trajectory_crossings = []

        idx_and_lasso = (0, -3) if self.is_detailed_out.get() else (0, -1)

        with open(self.file_with_trajectory) as f:
            for idx, line in enumerate(f):
                if idx == 1:
                    el = line.split()
                    self.trajectory_chain_range = [el[-4], el[-3]]
                    self.trajectory_chain_loop_indexes = [el[-2], el[-1]]
                elif idx > 7:
                    elem = line.split()
                    if elem.__contains__("ERROR"):
                        self.retrieved_frames.append(elem[0])
                        self.retrieved_trajectory_lassos.append("ERR")
                        self.retrieved_trajectory_crossings.append("ERR")
                    else:
                        self.retrieved_frames.append(elem[idx_and_lasso[0]])
                        self.retrieved_trajectory_lassos.append(elem[idx_and_lasso[-1]])
                        cross = self.get_crossings_from_trajectory(elem)
                        self.retrieved_trajectory_crossings.append(cross)

    def get_crossings_from_trajectory(self, elemlist):
        if elemlist[1] == "ERROR" or (int(float(elemlist[1])) is 0 and int(float(elemlist[2])) is 0):
            return "|"

        if self.is_detailed_out.get():
            n_crossings = (5, 5 + int(elemlist[1]))
            c_crossings = (5 + int(elemlist[1]) + 1, 5 + int(elemlist[1]) + int(elemlist[2]) + 1)
        else:
            n_crossings = (4, 4 + int(elemlist[1]))
            c_crossings = (4 + int(elemlist[1]) + 1,  4 + int(elemlist[1]) + int(elemlist[2]) + 1)

        crossings = []
        if int(float(elemlist[1])) is not 0:
            for i in range(n_crossings[0], n_crossings[1]):
                crossings.append(elemlist[i])
        if int(float(elemlist[2])) is not 0:
            for i in range(c_crossings[0], c_crossings[1]):
                crossings.append(elemlist[i])
        return crossings

    def set_trajectory_analysis_log(self):
        self.traj_log = Pmw.ScrolledText(self.win_trajectory_analysis_log.interior(), usehullsize=1, hull_width=625,
                                         hull_height=220)
        self.traj_log.importfile(self.file_with_trajectory)
        self.traj_log.grid(column=1, row=1)
        self.traj_log.configure(text_state='disabled')

    def calculate_lasso_in_trajectory(self):
        subprocess.Popen(self.python_compiler + [self._full_path_to_file],
                         stdout=subprocess.PIPE).communicate()[0].splitlines()

        self.update_trajectory_name("lasso")
        tmp_filename = self._filename + "_" + self.chains[0] + "_lasso.xyz"
        self.user_data = [self.program_execution + self._full_path_to_dir + os.sep + tmp_filename + " " +
                          self.trajectory_chain_loop_indexes[0] + " " + self.trajectory_chain_loop_indexes[1] +
                          " " + self.get_trajectory_advanced()]
        self.call_lasso_detection()
        self.move_files_to_polymer_directory()

    def update_trajectory_name(self, name):
        for f in os.listdir(self._full_path_to_dir):
            if f.__contains__(self._filename) and f != self._filename \
                    and not (f.endswith("_.pdb") or f.endswith("_.xyz") or f.__contains__("_lasso")):
                os.rename(self._full_path_to_dir + os.sep + f,
                          self._full_path_to_dir + os.sep + f[:-4] + "_" + name + f[-4:])

    def get_trajectory_advanced(self):
        adv = ["-f", "2"]

        if len(self.smooth_val.getvalue()) > 0:
            adv.append("-sm_nr")
            adv.append(str(self.smooth_val.getvalue()))
        if self.is_bad_caca_enabled.get():
            adv.append("-cd")
            adv.append("0")
        if not self.is_stable.get():
            adv.append("-redAC")
            adv.append(self.min_dist_crossings.getvalue())
            adv.append("-redBr")
            adv.append(self.min_dist_cross_end.getvalue())
            adv.append("-redEnd")
            adv.append(self.min_dist_cross_loop.getvalue())
        return " ".join(adv)

    def draw_error_charts(self, text, *charts):
        chart_lassos_type = mplt.figure.Figure(figsize=(5, 2), dpi=65, facecolor="lightgray")
        chart_lassos_type.subplots_adjust(top=0.92, bottom=0.26, left=0.07, right=0.97)
        ax = chart_lassos_type.add_subplot(111)
        ax.annotate(text, xy=(0.5, 0.5), xytext=(0.5, 0.5), ha="center")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        for i in charts:
            canvas = FigureCanvasTkAgg(chart_lassos_type, master=i)
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            canvas.show()

    def draw_lassos_type_chart(self):
        chart_lassos_type = mplt.figure.Figure(figsize=(6, 4), dpi=90, facecolor="lightgray")
        chart_lassos_type.subplots_adjust(top=0.96, bottom=0.27, left=0.07, right=0.98)
        ax = chart_lassos_type.add_subplot(111)
        ax.set_xlabel('Frame')
        ax.set_ylabel('Lasso type')

        xticks_divider = 140.0 if not self.is_pymol_2 else 100.0
        xticks_dist = int(ceil(len(self.retrieved_frames) / xticks_divider))
        self.lasso_info_tuple = {}
        x_min = float(self.retrieved_frames[0])
        x_max = float(self.retrieved_frames[-1])
        tmp_pos = []
        for idx, elem in enumerate((set(self.retrieved_trajectory_lassos))):
            try:
                tmp_pos.append(lassos.index(elem))
            except Exception:
                tmp_pos.append(lassos.index('Other'))
        tmp_pos.sort()
        for idx, elem in enumerate(tmp_pos):
            self.lasso_info_tuple[lassos[elem]] = [idx, colors[elem]]

        # configure x and y axis labels and rotation
        ax.tick_params(axis='both', labelsize=9)
        for tick in ax.get_xticklabels():
            tick.set_rotation(50)
        y_values = [self.lasso_info_tuple[elem][0] for elem in list(self.lasso_info_tuple.keys())]
        ax.yaxis.set_ticks(y_values)
        ax.set_ylim([-1, len(set(self.retrieved_trajectory_lassos))])
        for elem in list(self.lasso_info_tuple.keys()):
            pos = y_values.index(self.lasso_info_tuple[elem][0]) # change integer values to string equal to type of lasso
            y_values[pos] = elem
        ax.set_yticklabels(y_values)
        if not self.is_pymol_2:
            ax.set_xlim([x_min, x_max])
            ax.set_ylim([-1, len(set(self.retrieved_trajectory_lassos))])

        # draw chart, where x - frames and y - types of lasso
        for idx in range(0, len(self.retrieved_frames[:-1]), xticks_dist):
            if not self.retrieved_trajectory_lassos[idx] == "ERR" and idx < len(self.retrieved_frames[:-1]):
                try:
                    color = self.lasso_info_tuple[self.retrieved_trajectory_lassos[idx]][1]
                    ax.plot(self.retrieved_frames[idx], self.lasso_info_tuple[self.retrieved_trajectory_lassos[idx]][0],
                            marker="o", c=color, picker=3)
                except KeyError:
                    color = self.lasso_info_tuple['Other'][1]
                    ax.plot(self.retrieved_frames[idx], self.lasso_info_tuple['Other'][0], marker="o", c=color,
                            picker=3)
            else:
                color = self.lasso_info_tuple['ERR'][1]
                ax.plot(self.retrieved_frames[idx], self.lasso_info_tuple['ERR'][0], marker="o", c=color, picker=3)

        canvas = FigureCanvasTkAgg(chart_lassos_type, master=self.win_lasso_type.interior())
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvas.show()
        self.create_annotations(ax.lines)
        canvas.mpl_connect('pick_event', self.display_frame_in_pymol_on_pick)

    def view_possible_trajectory_loop(self):
        if not cmd.get_names(type="selections").__contains__("sele") and not len(self.output_data):
            self.raise_popup_menu('No atoms has been chosen.')

        for i in cmd.get_names(type="all"):
            if str(i).startswith("TMP") or str(i).startswith("DIST_") or str(i).startswith("CHAIN_") or \
                    str(i).startswith("SEQ") or str(i).startswith("PIERC") or str(i).startswith("TRIANG") or \
                    str(i).startswith("BR*") or str(i).startswith("SMOOTH_CHAIN_"):
                cmd.hide(representation="everything", selection=i)
                cmd.delete(name=i)

        atom = "ca"
        chain = self.chains[0]
        cmd.set(name="seq_view", value=1)
        cmd.set(name="mouse_selection_mode", value=6)
        cmd.hide(representation="sticks", selection="all")
        cmd.spectrum(palette="rainbow", selection="all")
        cmd.select("TMP_BR_" + self.trajectory_chain_loop_indexes[0] + "_" + self.trajectory_chain_loop_indexes[1],
                   selection="(chain " + chain + " and residue " + self.trajectory_chain_loop_indexes[0] +
                             " and name " + atom + ")+(chain " + chain + " and residue " +
                             self.trajectory_chain_loop_indexes[1] + " and name " + atom + ")")
        cmd.select("TMP_SEQ",
                   selection="chain " + chain + " and residue " + self.trajectory_chain_loop_indexes[0] + "-" +
                             self.trajectory_chain_loop_indexes[1])
        cmd.bond(atom1="chain " + chain + " and residue " + self.trajectory_chain_loop_indexes[0] + " and name " + atom,
                 atom2="chain " + chain + " and residue " + self.trajectory_chain_loop_indexes[1] + " and name " + atom)

        cmd.color(color="gray", selection="TMP_SEQ")
        cmd.show(representation='sphere', selection="TMP_BR_" + self.trajectory_chain_loop_indexes[0] + "_" +
                                                    self.trajectory_chain_loop_indexes[1])
        cmd.show(representation="sticks", selection="TMP_BR_" + self.trajectory_chain_loop_indexes[0] + "_" +
                                                    self.trajectory_chain_loop_indexes[1])
        cmd.show(representation="lines", selection="TMP_SEQ")
        cmd.set(name="sphere_color", value="orange", selection="all")
        cmd.set(name="stick_color", value="orange", selection="all")
        cmd.set("sphere_scale", value=0.5)
        cmd.deselect()

        if len(cmd.get_names(type="selections")) > 0:
            cmd.delete(name="sele")

    def display_frame_in_pymol_on_pick(self, event):
        if isinstance(event.artist, Line2D):
            if hasattr(self, "lasinf_surface_button") and self.lasinf_surface_button.winfo_exists():
                self.lasinf_surface_button.configure(state="disabled")
                self.lasinf_smooth_button.configure(state="disabled")
            artist = event.artist
            selected_frame = artist.get_xdata()[0]
            pos_frame = self.retrieved_frames.index(selected_frame)
            crossings = self.retrieved_trajectory_crossings[pos_frame]
            step = int(self.step.getvalue()) if len(self.step.getvalue()) != 0 else 1

            cmd.set(name="state", value=(pos_frame+1)*step)
            self.delete_crossings_selections()
            if crossings != "|":
                self.mark_crossings_on_trajectory(crossings)

    def create_annotations(self, artist):
        if not mplt.cbook.iterable(artist):
            artist = [artist]
        self.display_all = False
        self.axes = list(set(art.axes for art in artist))
        self.figures = set(ax.figure for ax in self.axes)
        self.annotations = {}

        for ax in self.axes:
            self.annotations[ax] = self.annotate(ax)
        for art in artist:
            art.set_picker(4)
        for fig in self.figures:
            fig.canvas.mpl_connect('pick_event', self.display_annotation)

    def annotate(self, ax):
        annotation = ax.annotate('', xy=(0, 0), xytext=(15, 15), textcoords='offset points', va='bottom', ha="left",
                                 bbox=dict(boxstyle='round,pad=0.5', fc='lightblue', alpha=0.5),
                                 arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        annotation.set_visible(False)
        return annotation

    def display_annotation(self, event):
        x, y = event.artist.get_xdata()[0], event.artist.get_ydata()[0]
        annotation = self.annotations[event.artist.axes]

        if x is not None:
            if not self.display_all:
                for ann in list(self.annotations.values()):
                    ann.set_visible(False)
            annotation.xy = x, y
            pos = list(self.axes[0].get_yticks()).index(y)
            type_lasso = list(self.axes[0].get_yticklabels())[pos].get_text()

            if type_lasso == "ERR":
                annotation.set_text("ERROR!\nFrame: %s" % (x))
            elif type_lasso == 'L0':
                annotation.set_text("Type of lasso: %s\nFrame: %s" % ('L0', x))
            else:
                pos = self.retrieved_frames.index(x)
                crossings = self.retrieved_trajectory_crossings[pos]
                cross = ""
                for i in crossings:
                    cross += i + " "
                if type_lasso == "Other":
                    type_lasso = self.retrieved_trajectory_lassos[int(float(x))]
                annotation.set_text("Type of lasso: %s\nPiercings: %s\nFrame: %s" % (type_lasso, cross, x))

            if int(float(x)) > int(float(self.retrieved_frames[len(self.retrieved_frames) / 2])):
                annotation.set_ha("right")
            else:
                annotation.set_ha("left")
            y_middle = max(list(self.axes[0].get_yticks()))
            if y > ceil(y_middle/2):
                annotation.set_y(-60)
                annotation.set_x(20) if annotation.get_ha() == "left" else annotation.set_x(-20)
            else:
                annotation.set_y(20)
                annotation.set_x(20) if annotation.get_ha() == "left" else annotation.set_x(20)
            annotation.set_visible(True)
            event.canvas.draw()

    def mark_crossings_on_trajectory(self, crossings):
        if crossings.__contains__("|"):
            return

        atom = "ca"
        pos_pierc = ""
        neg_pierc = ""

        if crossings.__contains__(",") or (type(crossings) == list and len(crossings) == 1):
            piercings = list(filter(len, crossings))
            for i in piercings:
                if i.startswith("+"):
                    pos_pierc += "(residue " + str(i[1:]) + " and name " + atom + ") "
                    pos_pierc += "(residue " + str(int(i[1:]) + 1) + " and name " + atom + ") "
                elif i.startswith("-"):
                    neg_pierc += "(residue " + str(i[1:]) + " and name " + atom + ") "
                    neg_pierc += "(residue " + str(int(i[1:]) + 1) + " and name " + atom + ") "
        else:
            for i in crossings:
                if i.startswith("+"):
                    if len(i) >= 3:
                        pos_pierc += "(residue " + str(i[1:]) + " and name " + atom + ") "
                        pos_pierc += "(residue " + str(int(i[1:]) + 1) + " and name " + atom + ") "
                    else:
                        pos_pierc += "(residue " + str(i[-1]) + " and name " + atom + ") "
                        pos_pierc += "(residue " + str(int(i[-1]) + 1) + " and name " + atom + ") "
                elif i.startswith("-"):
                    if len(i) >= 3:
                        neg_pierc += "(residue " + str(i[1:]) + " and name " + atom + ") "
                        neg_pierc += "(residue " + str(int(i[1:]) + 1) + " and name " + atom + ") "
                    else:
                        neg_pierc += "(residue " + str(i[-1]) + " and name " + atom + ") "
                        neg_pierc += "(residue " + str(int(i[-1]) + 1) + " and name " + atom + ") "

        if len(pos_pierc) > 0:
            cmd.select(name="POS_PIERC", selection=pos_pierc[:-1])
            cmd.color(color="lightblue", selection="POS_PIERC")
        if len(neg_pierc) > 0:
            cmd.select(name="NEG_PIERC", selection=neg_pierc[:-1])
            cmd.color(color="palegreen", selection="NEG_PIERC")
        cmd.deselect()

    def delete_crossings_selections(self):
        if None is not self.displayed_lasso:
            self.view_possible_trajectory_loop()

        selections = cmd.get_names(type="selections")
        for i in cmd.get_names(type="all"):
            if str(i).startswith("POS_PIERC") or str(i).startswith("NEG_PIERC") or str(i).startswith("SEQ") or \
                    str(i).startswith("TRIANG") or str(i).startswith("PIERC") or str(i).startswith("TMP_BR") or \
                    str(i).startswith("SMOOTH_CHAIN_") or str(i).startswith("CHAIN_"):
                cmd.delete(name=i)

        cmd.spectrum(palette="rainbow", selection="all")
        if selections.__contains__("TMP_SEQ"):
            cmd.color(color="gray", selection="TMP_SEQ")
        cmd.show(representation="lines", selection=cmd.get_names(type="objects")[0])
        if self.is_trajectory and len(self.chains) >= 2:
            self.display_first_chain_in_trajectory()

    def draw_atoms_piercing_lasso_chart(self):
        chart_atoms_piercing = mplt.figure.Figure(figsize=(5, 4), dpi=90, facecolor="lightgray")
        chart_atoms_piercing.subplots_adjust(top=0.96, bottom=0.27, left=0.06, right=0.98)
        ax = chart_atoms_piercing.add_subplot(111)
        ax.set_xlabel('Frame')
        ax.set_ylabel('Atom index')

        xticks_divider = 140.0 if not self.is_pymol_2 else 100.0
        xticks_dist = int(ceil(len(self.retrieved_frames) / xticks_divider))
        y_max = int(float(self.trajectory_chain_range[0]))
        y_min = int(float(self.trajectory_chain_range[-1]))
        ax.tick_params(axis='both', labelsize=9)
        for tick in ax.get_xticklabels():
            tick.set_rotation(50)

        # draw chart, where x - frames and y - atom crossing
        zipped_frames_crossings = list(zip(self.retrieved_frames, self.retrieved_trajectory_crossings))
        rect_width = 0
        for idx in range(0, len(self.retrieved_frames[:-1]), xticks_dist):
            if not zipped_frames_crossings[idx][1].__contains__("|") or zipped_frames_crossings[idx][1] != "|":
                if not zipped_frames_crossings[idx][1].__contains__("ERR"):
                    for i in range(len(zipped_frames_crossings[idx][1])):
                        color = "#008000" if zipped_frames_crossings[idx][1][i][0] == "+" else "#0000FF"
                        ax.plot(self.retrieved_frames[idx], float(zipped_frames_crossings[idx][1][i][1:]), marker="o",
                                c=color)
            elif self.is_pymol_2:
                ax.plot(self.retrieved_frames[idx], 0, marker="o", c="white", markeredgewidth=0.0)
            rect_width += 1

        if not self.is_pymol_2:
            x_min = float(self.retrieved_frames[0])
            x_max = float(self.retrieved_frames[-1])
            ax.set_xlim([x_min, x_max])
            ax.set_ylim([y_max, y_min])
            # draw orange rectangle
            rect_x = int(self.trajectory_chain_loop_indexes[0])
            rect_y = int(self.trajectory_chain_loop_indexes[-1])
            ax.add_patch(
                Rectangle((x_min - 1, rect_x), x_max - x_min + 2, rect_y - rect_x, facecolor="orange", linewidth=0))
        else:
            x_min = int(float(self.retrieved_frames[0]))
            # draw orange rectangle
            rect_x = int(self.trajectory_chain_loop_indexes[0])
            rect_y = int(self.trajectory_chain_loop_indexes[-1])
            ax.add_patch(Rectangle((x_min, rect_x), rect_width, rect_y - rect_x, facecolor="orange", linewidth=0))

        canvas = FigureCanvasTkAgg(chart_atoms_piercing, master=self.win_atoms_piercing_lasso.interior())
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canvas.show()

    ####################################################################################################################
    #                                   CALCULATE DATA FOR INSERTED FRAME
    ####################################################################################################################

    def calculate_trajectory_frame(self):
        decimal.getcontext().prec = 5

        tmp_frames_validate = []
        for i in self.frame_list:
            frame = i.getvalue()
            if len(frame) == 0:
                self.raise_popup_menu('There are empty fields in frame fields. Please fill all of them.')
            if frame.__contains__("."):
                num_zeros = len(frame) - frame.index(".")
                frame += ("0" * (6 - num_zeros))
            else:
                frame += ".00000"
            if not self.retrieved_frames.__contains__(frame):
                self.raise_popup_menu('No such frame (' + str(i.getvalue()) + ') in the above trajectory or there is an '
                                      'error in the selected frame. Please try another frame.')
            tmp_frames_validate.append(frame)

        tmp_frames_command = []
        self.frames_to_invoke = []
        for i in tmp_frames_validate:
            self.given_frame = i
            os.chdir(self._full_path_to_dir)
            self.create_file_containing_frame()
            self.create_frame_file_dir()

            detailed = (" -sframe " + ("2 " if self.is_detailed_alg.get() else "1 ")) \
                if not self.is_detailed_out_frame.get() else ""
            command = self.program_execution + self.file_frame_name + " " + self.trajectory_chain_loop_indexes[0] \
                      + " " + self.trajectory_chain_loop_indexes[1] + " " + self.get_trajectory_advanced() + detailed
            if not self.frames_to_invoke.__contains__(command):
                self.frames_to_invoke.append(command)
                tmp_frames_command.append(i)

        self.given_frames = tmp_frames_command

        self.user_data = self.frames_to_invoke

        self.call_lasso_detection()

        os.chdir(self._full_path_to_dir)
        tmp = self.current_working_dir
        self.current_working_dir = self._full_path_to_dir

        for i in self.given_frames:
            self.separate_files_to_directory(self._full_path_to_dir, os.getcwd() + os.sep +
                                             self._filename.replace(".", "_") + os.sep + "frame_" + str(i),
                                             "frame_" + str(i))
        self.current_working_dir = tmp
        os.chdir(self.current_working_dir)

        for line in self.array_of_results:
            for row in line:
                row.grid_forget()
        for i in self.btns_view_details:
            i.grid_forget()

        if hasattr(self, "lasinf_shallow_lasso_button") and self.lasinf_shallow_lasso_button.winfo_exists():
            self.lasinf_shallow_lasso_button.grid_forget()
        if hasattr(self, "lasinf_surface_button") and self.lasinf_surface_button.winfo_exists():
            self.lasinf_surface_button.grid_forget()
        if hasattr(self, "lasinf_smooth_button") and self.lasinf_smooth_button.winfo_exists():
            self.lasinf_smooth_button.grid_forget()

        self.display_lasso_information_table()
        self.create_lasso_information_buttons()
        self.win_trajectory_analysis.yview('moveto', 0)

    def create_file_containing_frame(self):
        self.file_frame_name = self._full_path_to_dir + os.sep + self._filename + "_" + self.chains[0] + "_" + \
                               "_frame_" + str(self.given_frame) + ".xyz"
        frame_file = open(self.file_frame_name, "w")
        frame_found = False
        path_to_traj_dir = self._full_path_to_dir + os.sep + self._filename.replace(".", "_") + os.sep + \
                           self._filename + "_" + self.chains[0] + ".xyz"

        with open(path_to_traj_dir) as f:
            for idx, line in enumerate(f):
                if frame_found:
                    frame_file.write(line[:-1] + " ABC\n")
                if line.__contains__(str(self.given_frame)):
                    frame_found = True
                elif line.__contains__("t"):
                    frame_found = False
        frame_file.close()

    def create_frame_file_dir(self):
        frame_dir = self._full_path_to_dir + os.sep + self._filename.replace(".", "_") + os.sep + \
                    "frame_" + str(self.given_frame)

        if not os.path.exists(frame_dir):
            os.mkdir(frame_dir)
            print("  Creating separate directory for frame data...")
        else:
            print("  Directory found. Replacing current frame files with new data...")

    ####################################################################################################################
    #                                       2. CREATE PROTEIN RESULT WINDOW
    ####################################################################################################################

    def create_protein_window(self):
        self.win_lasso_info = Pmw.Dialog(self.parent, title='Lasso information window [' + self._filename +
                                                            ' chain ' + self.chain_index.get() + ']',
                                         buttons=["Show bridges", "Exit"],
                                         command=self.invoke_protein_window_buttons)
        self.win_lasso_info.resizable(0, 0)
        self.show_bridges()

        if not all(x.__contains__("ERROR") for x in self.output_data):
            if not all((len(elem[0]) == 0 and len(elem[1]) == 0) for elem in self.lassos):
                self.create_shallow_lasso_button()
            if self.is_gln_checkbutton_selected.get():
                self.create_gln_button()
            self.create_lasso_information_buttons()
            self.create_surface_hints()

    def show_bridges(self):
        self.win_lasso_info.invoke("Show bridges")
        self.win_lasso_info.component("buttonbox").delete(0)

    def invoke_protein_window_buttons(self, clicked_button):
        if clicked_button == "Show bridges":
            self.window_parent = self.win_lasso_info.interior()
            self.display_lasso_information_table()
        else:
            if hasattr(self, "lasinf_smooth_button") and self.lasinf_smooth_button.winfo_exists():
                self.lasinf_smooth_button.grid_forget()
            if hasattr(self, "lasinf_surface_button") and self.lasinf_surface_button.winfo_exists():
                self.lasinf_surface_button.grid_forget()
            if hasattr(self, "lasinf_shallow_lasso_button") and self.lasinf_shallow_lasso_button.winfo_exists():
                self.lasinf_shallow_lasso_button.grid_forget()
            if hasattr(self, "lasinf_gln_button") and self.lasinf_gln_button.winfo_exists():
                self.lasinf_gln_button.deselect()
                self.lasinf_gln_button.grid_forget()
            self.win_lasso_info.withdraw()
            self.delete_pymol_objects()
            self.display_pymol_chain()

    def create_shallow_lasso_button(self):
        self.lasinf_shallow_display = tk.IntVar()
        self.lasinf_shallow_lasso_button = tk.Checkbutton(self.window_parent, text="Shallow lasso", anchor="center",
                                                          indicatoron=0, width=12, height=2, relief="ridge",
                                                          variable=self.lasinf_shallow_display,
                                                          command=self.pymol_display_shallow_lassos)
        hint_shallow_lasso = Pmw.Balloon(self.window_parent, relmouse="both")
        hint_shallow_lasso.bind(self.lasinf_shallow_lasso_button, textwrap.fill("Shows/hides shallow lassos. If shallow "
                                                                              "lassos are not visible press "
                                                                              "'View details' button.",
                                                                              self.hint_width))

    def create_gln_button(self):
        self.lasinf_is_gln_selected = tk.IntVar()
        self.lasinf_gln_button = tk.Checkbutton(self.window_parent, text="GLN", anchor="center", indicatoron=0,
                                                width=10, height=2, relief="ridge",
                                                variable=self.lasinf_is_gln_selected,
                                                command=self.display_gln_objects)
        hint_gln_button = Pmw.Balloon(self.window_parent, relmouse="both")
        hint_gln_button.bind(self.lasinf_gln_button, "Shows/Hides GLN matrices.")

    def create_lasso_information_buttons(self):
        self.lasinf_surface_button = tk.Checkbutton(self.window_parent, text="Surface", anchor="center",
                                                    indicatoron=0, width=10, height=2, relief="ridge",
                                                    variable=self.lasinf_surface_display,
                                                    command=self.pymol_display_surface)
        self.lasinf_surface_button.select()
        self.lasinf_smooth_button = tk.Checkbutton(self.window_parent, text="Smooth", anchor="center",
                                                   indicatoron=0, width=10, height=2, relief="ridge",
                                                   variable=self.lasinf_smooth_display,
                                                   command=self.pymol_display_smooth)

        if (hasattr(self, "lasinf_shallow_lasso_button") and self.lasinf_shallow_lasso_button.winfo_exists()) \
                and (hasattr(self, "lasinf_gln_button") and self.lasinf_gln_button.winfo_exists()):
            self.lasinf_surface_button.grid(column=2, row=len(self.output_data) + 1, padx=5, pady=8)
            self.lasinf_smooth_button.grid(column=3, row=len(self.output_data) + 1, padx=5, pady=8, )
            self.lasinf_shallow_lasso_button.grid(column=4, row=len(self.output_data) + 1, padx=5, pady=8)
            self.lasinf_gln_button.grid(column=7, row=len(self.output_data) + 1, padx=5, pady=8)
        elif (hasattr(self, "lasinf_shallow_lasso_button") and self.lasinf_shallow_lasso_button.winfo_exists()) \
                and not (hasattr(self, "lasinf_gln_button") and self.lasinf_gln_button.winfo_exists()):
            self.lasinf_surface_button.grid(column=2, row=len(self.output_data) + 1, padx=10, pady=8)
            self.lasinf_smooth_button.grid(column=3, row=len(self.output_data) + 1, padx=10, pady=8)
            self.lasinf_shallow_lasso_button.grid(column=4, row=len(self.output_data) + 1, pady=8)
        if (hasattr(self, "lasinf_gln_button") and self.lasinf_gln_button.winfo_exists()) \
                and not (hasattr(self, "lasinf_shallow_lasso_button") and self.lasinf_shallow_lasso_button.winfo_exists()):
            self.lasinf_surface_button.grid(column=2, row=len(self.output_data) + 1, padx=10, pady=8, columnspan=2)
            self.lasinf_smooth_button.grid(column=4, row=len(self.output_data) + 1, padx=10, pady=8, columnspan=2)
            self.lasinf_gln_button.grid(column=6, row=len(self.output_data) + 1, padx=10, pady=8, columnspan=2)
        else:
            self.lasinf_surface_button.grid(column=2, row=len(self.output_data) + 1, padx=10, pady=8, columnspan=2)
            self.lasinf_smooth_button.grid(column=5, row=len(self.output_data) + 1, padx=10, pady=8, columnspan=2)

    def create_surface_hints(self):
        hint_surface = Pmw.Balloon(self.window_parent, relmouse="both")

        hint_surface.bind(self.lasinf_surface_button,
                          "Shows/hides the minimal surface spanned on the closed loop.")
        hint_surface.bind(self.lasinf_smooth_button,
                          "Shows/hides the smoothed representation of the chain.")

    ####################################################################################################################
    #                                     DISPLAY LASSO INFORMATION TABLE
    ####################################################################################################################

    def display_lasso_information_table(self):
        self.lasinf_smooth_display = tk.IntVar()
        self.lasinf_surface_display = tk.IntVar()
        self.lasinf_shallow_display = None
        self.displayed_lasso = None
        self.prev_displayed_lasso = None
        self.lassos = []
        self.array_of_results = []
        cell_height = 1

        self.create_row_names()
        self.textwrap_error_list()
        self.list_img_lasso = self.get_lasso_images()
        chain_ends = self.get_chain_atoms()

        if hasattr(self, "type_loop_closing_bridge") and str(
                self.type_loop_closing_bridge.get()) != "choose two atoms to form a bridge":
            self.list_img_bridges = self.get_bridge_images()

        if len(self.output_data) == 0:
            no_lassos = tk.Label(self.window_parent, bg="white", width=14,
                                 text="No lasso found in a given protein chain .")
            no_lassos.grid(sticky="swen", column=1, columnspan=9, row=0)
            return

        for idx, elem in enumerate(self.output_data):
            elem = elem.split("SMOOTH")[0]
            row_element = []
            if elem.__contains__("ERROR"):
                row_element.append(self.create_error_loop_range(idx, cell_height))
                row_element.append(self.create_error(cell_height, elem))
            else:
                split_elements = re.split(' ', elem)
                n_end_length = self.get_n_crossings(split_elements)
                c_end_length = self.get_c_crossings(split_elements)

                l = self.itemise_types_of_lasso(split_elements, n_end_length.split(" "), c_end_length.split(" "))
                self.lassos.append(l)

                n_end_length = self.textwrap_n_crossings(n_end_length)
                c_end_length = self.textwrap_c_crossings(c_end_length)

                cell_height = max(textwrap.fill(n_end_length, 16).count("\n"),
                                  textwrap.fill(c_end_length, 16).count("\n")) + 1
                for idx2, clicked_button in enumerate(
                        (split_elements[-4], (split_elements[1] + "-" + split_elements[2]),
                         int(split_elements[2]) - int(split_elements[1]) + 1,
                         n_end_length, c_end_length, int(split_elements[1]) - int(chain_ends[0]),
                         int(chain_ends[1]) - int(split_elements[2]), str(int(float(split_elements[-2]))))):
                    if clicked_button == (split_elements[1] + "-" + split_elements[2]):
                        if hasattr(self, "warning_gaps") and self.warning_gaps:
                            not_found = True
                            i = 0
                            chain = self.chain_index.get()
                            while i < self.warning_gaps.__len__() and not_found:
                                if chain == self.warning_gaps[i][0] and int(split_elements[1]) <= \
                                        int(self.warning_gaps[i][1]) and int(split_elements[2]) >= \
                                        int(self.warning_gaps[i][2]):
                                    row_element.append(self.create_loop_range_interior(clicked_button, cell_height))
                                    not_found = False
                                i += 1
                            if not_found:
                                row_element.append(
                                    tk.Text(self.window_parent, bg="white", height=cell_height + 1, width=14,
                                            padx=0, pady=0, wrap="word", bd=0, highlightthickness=0))
                                self.configure_tkinter_text(row_element[-1], str(clicked_button), cell_height)
                        else:
                            row_element.append(
                                tk.Text(self.window_parent, bg="white", height=cell_height + 1, width=14,
                                        padx=0, pady=0, wrap="word", bd=0, highlightthickness=0))
                            self.configure_tkinter_text(row_element[-1], str(clicked_button), cell_height)
                    elif clicked_button == n_end_length or clicked_button == c_end_length:
                        row_element.append(
                            tk.Text(self.window_parent, height=cell_height + 1, width=17, padx=0, pady=0,
                                    wrap="word", bd=0,
                                    highlightthickness=0, background="white"))
                        self.color_crossings(row_element[-1], clicked_button, 0)
                    elif clicked_button == split_elements[-4]:
                        row_res_element = self.create_bridge_type_interior(split_elements, clicked_button, idx,
                                                                           cell_height)
                        row_element.append(row_res_element)
                    else:
                        row_element.append(
                            tk.Text(self.window_parent, bg="white", height=cell_height + 1, width=17,
                                    padx=0, pady=0, wrap="word", bd=0, highlightthickness=0))
                        self.configure_tkinter_text(row_element[-1], str(clicked_button), cell_height)
            self.array_of_results.append(row_element)

        self.display_results_table()
        self.create_view_details_buttons()

        if hasattr(self, "given_frames") and len(self.given_frames) != 0:
            self.display_selected_frames()

        if len(self.btns_view_details) == 1:
            self.btns_view_details[0].invoke()

    def create_row_names(self):
        rows = []

        for idx, i in enumerate(["Lasso type", "Loop range", "Loop length", "N-term piercings", "C-term piercings",
                                 "N-end length", "C-end length", "Lasso surface\narea"]):
            rows.append(tk.Label(self.window_parent, text=i, width=self.table_name_width, font=self.row_name_size,
                                 justify="center"))
            rows[-1].grid(column=1 + idx, row=0)
        rows[0].configure(width=9)
        rows[1].configure(width=9)
        rows[2].configure(width=9)

    def textwrap_error_list(self):
        for elem in errors:
            errors[elem] = textwrap.fill(errors[elem], 140)

    def get_lasso_images(self):
        img_lassos = {}
        lasso_types = {"L0": "L0", "L1": "L1", "L2": "L2", "L3": "L3", "L4": "L4", "L5": "L5", "L6": "L6",
                       "LL1,1": "LL1,1", "LL1,2": "LL1,2", "LL2,1": "LL2,1", "LL2,2": "LL2,2", "LL1,4": "LL1,4",
                       "LL4,1": "LL4,1", "LL4,2": "LL4,2", "LL4,3": "LL4,3",
                       "LS2": "LS2", "LS3": "LS3", "LS4": "LS4", "LS5": "LS5", "LS7": "LS7",
                       "ERR": "ERR"  # for trajectories only
                       }

        for k, v in lasso_types.items():
            try:
                img_lassos[k] = tk.PhotoImage(file=os.path.join(plugin_path + os.sep + "lassos", v + self._img_extension))
            except:
                # if directory does not contain appropriate lasso type picture
                img_lassos[k] = tk.PhotoImage(file=os.path.join(plugin_path + os.sep + "lassos", "error"
                                                                + self._img_extension))
        tmp_label_lasso = tk.Label(self.window_parent)  # without label, pics will disappear
        tmp_label_lasso.image = img_lassos
        return img_lassos

    def get_chain_atoms(self):
        if self.is_trajectory:
            chain = self.chains[0]
        else:
            chain = self.chain_index.get()

        atoms = []
        reg = re.compile('ATOM\s\s+\d+')

        with open(self._full_path_to_dir + os.sep + self._filename, "r") as f:
            for line in f:
                clear_line = list(filter(len, line.split(" ")))
                if reg.match(line) and self.is_trajectory:
                    try:
                        atoms.append(int(clear_line[4]))
                    except:
                        atoms.append(clear_line[5])
                elif reg.match(line) and not self.is_trajectory:
                    if clear_line[4].isalpha() and clear_line[4] == chain:
                        atoms.append(clear_line[5])
                    elif clear_line[4][0] == chain:
                        idx = clear_line[4][1:]
                        atoms.append(idx)

        return [int(atoms[0]), int(atoms[-1])]

    def get_bridge_images(self):
        img_bridges = {}

        bridges = {"SS": "SS", "AMIDE": "Amide", "AMIDE-like": "Amide", "ESTER": "Ester", "ESTER-like": "Ester",
                   "THIOESTER": "Thioester", "THIOESTER-like": "Thioester", "OTHER": "Other"}
        for k, v in bridges.items():
            img_bridges[k] = tk.PhotoImage(file=os.path.join(plugin_path + os.sep + "img", v + self._img_extension))
        tmp_label_bridge = tk.Label(self.fr_loop_closing_bridge.interior())
        tmp_label_bridge.image = img_bridges
        return img_bridges

    def create_error_loop_range(self, idx, spacing):
        loop_range = tk.Text(self.window_parent, bg="white", height=spacing + 1, width=14,
                             padx=0, pady=0, wrap="word", bd=0, highlightthickness=0)
        loop_range.insert("1.0", self.user_data[idx].split(" ")[2] + "-" + self.user_data[idx].split(" ")[3])
        loop_range.configure(state="disabled")
        loop_range.tag_add("center", 1.0, "end")
        loop_range.tag_configure("center", justify='center')
        loop_range.tag_add("spacing1", 1.0, "end")
        loop_range.tag_config("spacing1", spacing1=6 + spacing)
        return loop_range

    def create_error(self, spacing, txt):
        width = 115 if self.is_trajectory else 14
        error = tk.Text(self.window_parent, bd=0, height=spacing + 1, padx=5, pady=0, bg="white", width=width,
                        highlightthickness=0)
        error.insert("1.0", errors[int(re.split(' ', txt)[0][-3])])
        error.tag_add("center", 1.0, "end")
        error.tag_configure("center", justify='center')
        error.tag_add("spacing1", 1.0, "end")
        error.tag_config("spacing1", spacing1=2 + spacing)
        error.configure(state="disabled")
        return error

    def get_n_crossings(self, elem):
        n_end_length = ""

        if int(elem[4]) is not 0:
            for j in range(8, (8 + int(elem[4]))):
                n_end_length += elem.__getitem__(j) + ", "
        return n_end_length

    def get_c_crossings(self, elem):
        c_end_length = ""

        if int(elem[5]) is not 0:
            for j in range(5 + 2 + 2 + int(elem[4]), 5 + 2 + 2 + int(elem[4]) + int(elem[5])):
                c_end_length += elem.__getitem__(j) + ", "
        return c_end_length

    def textwrap_n_crossings(self, n):
        if None is not self.lassos[-1] and len(self.lassos[-1][0]) > 0:
            return textwrap.fill(" ".join(self.lassos[-1][0][1])[:-1], 16)
        else:
            return textwrap.fill(n[:-2], 16)

    def textwrap_c_crossings(self, c):
        if len(self.lassos[-1][1]) > 0:
            return textwrap.fill(" ".join(self.lassos[-1][1][1])[:-1], 16)
        else:
            return textwrap.fill(c[:-2], 16)

    def configure_tkinter_text(self, tk_text_elem, text, height):
        tk_text_elem.insert("1.0", str(text))
        tk_text_elem.tag_add("center", 1.0, "end")
        tk_text_elem.tag_configure("center", justify='center')
        tk_text_elem.tag_add("spacing1", 1.0, "end")
        tk_text_elem.tag_config("spacing1", spacing1=5 + height + 1)
        tk_text_elem.configure(state="disabled")

    def create_loop_range_interior(self, text, height):
        chain_img = "broken_chain"
        broken_chain_img = tk.PhotoImage(file=os.path.join(plugin_path + os.sep + "img", chain_img + self._img_extension))
        tmp_broken_chain_lasso = tk.Label(self.window_parent)  # without label, pics will disappear
        tmp_broken_chain_lasso.image = broken_chain_img

        row_res_element = tk.Text(self.window_parent, bg="white", height=height + 1, width=20,
                                  padx=4, pady=0, wrap="word", bd=0, highlightthickness=0)

        broken_chain_img = tk.Label(row_res_element, image=broken_chain_img, bg="white", height=16, width=13, bd=0)
        broken_chain_img.grid(column=0, row=0)

        hints_broken_chain = Pmw.Balloon(self.window_parent, relmouse="both")
        hints_broken_chain.bind(broken_chain_img, "The chain is broken within the loop. The broken part of the chain "
                                                  "has been\nreplaced by a straight segment, which may affect what "
                                                  "lasso types are\ndetected - be careful with interpreting the results")

        loop_range = tk.Text(row_res_element, bg="white", width=11, height=height + 1, padx=0, pady=0, wrap="word",
                             bd=0, highlightthickness=0)
        loop_range.grid_configure(column=1, row=0)

        self.configure_tkinter_text(loop_range, str(text), height - 1)
        row_res_element.configure(state="disabled")
        return row_res_element

    def create_bridge_type_interior(self, var, text, idx, height):
        hints_result_window = Pmw.Balloon(self.window_parent, relmouse="both")
        img_lasso = str(var[-4]).replace("+", "").replace("-", "").replace("N", "").replace("C", "")

        row_res_elem = tk.Text(self.window_parent, bg="white", width=15, height=height + 1, padx=0, pady=0,
                               wrap="word", bd=0, highlightthickness=0)
        try:
            lasso_img = tk.Label(row_res_elem, image=self.list_img_lasso[img_lasso], bg="white", height=20, width=20,
                                 bd=0)
            hints_result_window.bind(lasso_img, lasso_description[img_lasso])
        except KeyError: # for trajectories only if type of lasso is very complicated
            lasso_img = tk.Label(row_res_elem, image=self.list_img_lasso["ERR"], bg="white", height=20, width=20,
                                 bd=0)
            hints_result_window.bind(lasso_img, lasso_description["ERR"])
        lasso_img.grid(column=1, row=0, padx=4)

        if hasattr(self, "type_loop_closing_bridge") and not self.is_trajectory and str(
                self.type_loop_closing_bridge.get()) != "choose two atoms to form a bridge":
            bridge_img = tk.Label(row_res_elem, image=self.list_img_bridges[self.list_bridges[idx][0]],
                                  height=16, bg="white")
            type_of_bridge = self.list_bridges[idx][0]
            if type_of_bridge.__contains__("like"):
                type_of_bridge = type_of_bridge.replace("-like", "")
            hints_result_window.bind(bridge_img,
                                     (type_of_bridge.title() if len(type_of_bridge) > 2 else type_of_bridge) + " bridge")
            bridge_img.grid(column=0, row=0, padx=3)
            lasso_img.grid_configure(column=1, row=0, padx=3)

        type_lasso = tk.Text(row_res_elem, bg="white", height=height + 1, width=7, wrap="word", padx=0, pady=0, bd=0,
                             highlightthickness=0)
        type_lasso.tag_configure("subscript", offset=-5)

        lasso = self.simplify_type_of_lasso(text)

        type_lasso.insert("insert", lasso, "", str(text)[len(lasso):], "subscript")
        type_lasso.tag_add("font", 1.0 + float(len(lasso)) / 10, END)
        type_lasso.tag_configure("font", font=self.subscript_font)
        type_lasso.tag_add("left", 1.0, "end")
        type_lasso.tag_configure("left", justify='left')
        type_lasso.grid(column=2, row=0)
        type_lasso.tag_add("spacing1", 1.0, "end")
        type_lasso.tag_config("spacing1", spacing1=5 + height + 1)
        type_lasso.configure(state="disabled")
        row_res_elem.configure(state="disabled")
        return row_res_elem

    def simplify_type_of_lasso(self, text):
        if re.finditer('L|S', str(text)):
            l = [(m.start(0), m.end(0)) for m in re.finditer('L|S', str(text))]
            pos = [l[0][0], l[0][1]] if len(l) == 1 else [l[0][0], l[1][1]]
            return str(text)[int(pos[0]):int(pos[1])]

    def display_results_table(self):
        row = 0

        for i in self.array_of_results:
            row += 1
            if len(i) == 2:
                i[0].grid(column=2, row=row)
                if self.is_trajectory:
                    i[1].grid(sticky="swen", column=3, columnspan=10, row=row)
                else:
                    i[1].grid(sticky="swen", column=3, columnspan=9, row=row)
            else:
                for idx, elem in enumerate(i):
                    elem.grid(column=idx + 1, row=row)

    def display_gln_objects(self):
        if self.displayed_lasso is None:
            self.lasinf_gln_button.deselect()
            self.raise_popup_menu('No lasso chosen. Please load an appropriate lasso into the PyMOL viewer '
                                  '(e.g. press the button ''view details'') and try again.')

        if self.lasinf_is_gln_selected.get():
            print("  Matrices are loading. Please wait...")
            if hasattr(self, "win_gln_matrices") and self.win_gln_matrices.winfo_exists():
                self.win_gln_matrices.grid_forget()
            self.prev_displayed_gln_segment = None
            self.win_gln_matrices = Pmw.Group(self.window_parent, tag_text="Gaussian Linking Number Matrices "
                                                                           "(entanglement between the loop and the "
                                                                           "tails)")
            self.win_gln_matrices.grid(sticky="swen", column=0, row=len(self.output_data) + 2, columnspan=10,
                                       padx=5, pady=5)

            chain = self.chain_index.get()
            res_beg = self.output_data[self.displayed_lasso].split(" ")[1]
            res_end = self.output_data[self.displayed_lasso].split(" ")[2]

            file_path = self._filename.replace(".", "_")
            is_smoothed = "_smooth" if self.lasinf_smooth_display.get() else ""

            filename1 = self._full_path_to_dir + os.sep + file_path + os.sep + "_GLN" + os.sep + "matrixGLN_" + \
                        self._filename + "_" + chain + "_" + res_beg + "_" + res_end + "_t1" + is_smoothed + ".py"
            filename2 = self._full_path_to_dir + os.sep + file_path + os.sep + "_GLN" + os.sep + "matrixGLN_" + \
                        self._filename + "_" + chain + "_" + res_beg + "_" + res_end + "_t2" + is_smoothed + ".py"

            gln_matrix_1 = mplt.figure.Figure(figsize=(self.gln_figsize, self.gln_figsize), dpi=self.gln_dpi,
                                              facecolor='lightgray')
            ax1 = gln_matrix_1.add_subplot(111, aspect="equal", adjustable='box-forced')
            gln_matrix_1.subplots_adjust(top=0.97, bottom=0.17, left=0.11, right=0.97)
            gln_matrix_2 = mplt.figure.Figure(figsize=(self.gln_figsize, self.gln_figsize), dpi=self.gln_dpi,
                                              facecolor='lightgray')
            ax2 = gln_matrix_2.add_subplot(111, aspect="equal", adjustable='box-forced')
            gln_matrix_2.subplots_adjust(top=0.97, bottom=0.17, left=0.11, right=0.97)

            is_displayed = [False, False]
            n_end = int(self.array_of_results[self.displayed_lasso][5].get("1.0", 'end-1c'))
            c_end = int(self.array_of_results[self.displayed_lasso][6].get("1.0", 'end-1c'))
            n_terminus_error = textwrap.fill("The GLN matrix has not been generated.The length of the N-terminus "
                                             "is less than 5 amino acids.", 35)
            c_terminus_error = textwrap.fill("The GLN matrix has not been generated.The length of the C-terminus "
                                             "is less than 5 amino acids.", 35)

            if n_end < 5:
                self.draw_gln_empty(ax1, n_terminus_error)
            else:
                imp.load_source("matrix_1", filename1)
                import matrix_1
                try:
                    matrix_1.draw_gln_matrices(ax1)
                    is_displayed[0] = True
                except:
                    self.draw_gln_empty(ax1, n_terminus_error)
            if c_end < 5:
                self.draw_gln_empty(ax2, c_terminus_error)
            else:
                imp.load_source("matrix_2", filename2)
                import matrix_2
                try:
                    matrix_2.draw_gln_matrices(ax2)
                    is_displayed[1] = True
                except:
                    self.draw_gln_empty(ax2, c_terminus_error)

            canvas = FigureCanvasTkAgg(gln_matrix_1, master=self.win_gln_matrices.interior())
            canvas.get_tk_widget().grid(column=0, row=0)
            canvas._tkcanvas.grid(column=0, row=0)
            canvas.show()
            canvas2 = FigureCanvasTkAgg(gln_matrix_2, master=self.win_gln_matrices.interior())
            canvas2.get_tk_widget().grid(column=1, row=0)
            canvas2._tkcanvas.grid(column=1, row=0)

            if any(is_displayed):
                hint_view = Pmw.Balloon(self.win_gln_matrices.interior(), relmouse="both")
                hint_view.bind(self.win_gln_matrices.interior(),
                               "Click on an element of the matix to highlight (in green) the\ncorresponding segment "
                               "of the tail (on the protein structure\nin the PyMOL viewer).")

            if is_displayed[0]:
                canvas.mpl_connect('button_press_event', self.show_gln_segment)
            if is_displayed[1]:
                canvas2.mpl_connect('button_press_event', self.show_gln_segment)
            canvas2.show()
            self.pymol_color_structure_gln()
        else:
            if hasattr(self, "win_gln_matrices") and self.win_gln_matrices.winfo_exists():
                self.win_gln_matrices.grid_forget()
            self.pymol_view_details(self.displayed_lasso)

    def draw_gln_empty(self, ax, text):
        ax.text(0.15, 0.5, text, fontsize=12)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    def show_gln_segment(self, event):
        if event.xdata < event.ydata:
            if self.prev_displayed_gln_segment is not None:
                self.decolor_gln_segment()
        self.color_gln_segment(event)
        cmd.deselect()

    def decolor_gln_segment(self):
        chain = self.chain_index.get()
        atom = "ca"

        for idx in range(self.prev_displayed_gln_segment[0], self.prev_displayed_gln_segment[1] + 1):
            atom_idx = "chain " + chain + " and residue " + str(idx) + " and name " + atom
            atom_idx_1 = "chain " + chain + " and residue " + str(int(idx) + 1) + " and name " + atom
            cmd.set_color(name=atom_idx, rgb=self.gln_protein_colors[idx])
            cmd.color(color=atom_idx, selection=atom_idx)
            cmd.color(color=atom_idx, selection=atom_idx_1)

    def color_gln_segment(self, event):
        if event.xdata and event.ydata:
            gln_tuple = (int(event.xdata), int(event.ydata))
            if gln_tuple[1] >= gln_tuple[0]:
                selected_atoms_gln = ""
                chain = self.chain_index.get()
                atom = "ca"
                chain = "" if self._file_extension == "xyz" else "chain " + chain + " and "

                for i in range(gln_tuple[0], gln_tuple[1] + 1):
                    selected_atoms_gln += "(" + chain + "residue " + str(i) + " and name " + atom + ") + "
                selected_atoms_gln = selected_atoms_gln[:-3]
                cmd.select(name="GLN_SELE", selection=selected_atoms_gln)
                cmd.color(color="green", selection="GLN_SELE")
                print("  Sequence " + str(gln_tuple[0]) + "-" + str(gln_tuple[1]) + " has been colored...")
            self.prev_displayed_gln_segment = (int(event.xdata), int(event.ydata))

    def create_view_details_buttons(self):
        self.btns_view_details = []

        for idx, elem in enumerate(self.output_data):
            if len(elem) is not 0 and not elem.__contains__("ERROR"):
                self.btns_view_details.append(tk.Button(self.window_parent, text="view details",
                                                        command=lambda x=idx: self.pymol_view_details(x)))
                self.btns_view_details[-1].grid(column=0, row=1 + idx)
        self.create_view_details_hints()

    def create_view_details_hints(self):
        hint_view = Pmw.Balloon(self.window_parent, relmouse="both")

        for i in self.btns_view_details:
            hint_view.bind(i, textwrap.fill("Displays the chain with the minimal surface spanned on the closed loop.",
                                            self.hint_width))

    def display_selected_frames(self):
        chosen_frames = "Given frames:"

        for i in self.given_frames:
            chosen_frames += "\n" + i
        if hasattr(self, "selected_frames") and self.selected_frames.winfo_exists():
            self.selected_frames.grid_forget()
        self.selected_frames = tk.Label(self.window_parent, text=chosen_frames, justify="center")
        self.selected_frames.grid(row=1 + len(self.output_data), rowspan=2, column=0)

    ####################################################################################################################
    #                                           VIEW CHOSEN LASSO
    ####################################################################################################################

    def color_line_in_array_of_results(self):
        for idx, elem in enumerate(self.array_of_results):
            if idx == self.displayed_lasso:
                for i in elem:
                    if len(i.winfo_children()):
                        for j in i.winfo_children():
                            j.configure(bg="lightblue")
                    i.configure(bg="lightblue")
            else:
                for i in elem:
                    if len(i.winfo_children()):
                        for j in i.winfo_children():
                            j.configure(bg="white")
                    i.configure(bg="white")

    def decolor_line_in_array_of_results(self):
        if self.prev_displayed_lasso:
            for idx, elem in enumerate(self.array_of_results):
                if idx == self.prev_displayed_lasso:
                    for i in elem:
                        if len(i.winfo_children()):
                            for j in i.winfo_children():
                                j.configure(bg="white")
                        i.configure(bg="white")

    def color_crossings(self, text, str_crossings, show_shallow=0, res_idx=None):
        if len(str_crossings) != 0:
            text.configure(state="normal")
            for i in range(int(float(text.index(END)))):
                text.delete(i + 1.0, "end-1c")
            text.delete(1.0, "end-1c")

            if type(str_crossings) == str:
                tmp_crossings = str_crossings.split("\n")
            elif type(str_crossings) == list and show_shallow == 1:
                deep_only = str_crossings[1]
                tmp_crossings = textwrap.fill(" ".join(str_crossings[0])[:-1], 15).split("\n")
            else:
                tmp_crossings = textwrap.fill(" ".join(str_crossings[1])[:-1], 16).split("\n")

            repetitive_crossings = []
            column = 0.0
            for line in tmp_crossings:
                column += 1.0
                text.insert(column, line + "\n")
                i = 0
                while i < len(line):
                    piercing = ""
                    if line[i] == "-" or line[i] == "+":
                        beg = i
                        color = self.set_piercing_color(line[i])
                        while i < len(line) and line[i] != ",":
                            piercing += line[i]
                            i += 1
                        if show_shallow == 1:
                            piercing += ","
                            if not deep_only.__contains__(piercing) or repetitive_crossings.__contains__(piercing):
                                color = "gray"
                        text.tag_add(color,
                                     ("%.1f" % (float(beg) / 10 + column)) if len(str(beg)) == 1 else
                                     ("%.2f" % ((float(beg)) / 100 + column)),
                                     ("%.1f" % (float(i) / 10 + column)) if len(str(i)) == 1 else
                                     ("%.2f" % ((float(i)) / 100 + column)))
                        text.tag_config(color, foreground=color)
                        repetitive_crossings.append(piercing)
                    i += 1
            text.tag_add("left", 1.0, "end-1c")
            text.tag_configure("left", justify='left')
            text.tag_add("spacing1", 1.0, "end-1c")
            text.tag_config("spacing1", spacing1=9 - (3 * len(tmp_crossings)))

            if res_idx is not None:
                for elem in self.array_of_results[res_idx]:
                    if len(elem.winfo_children()) == 2:
                        elem.winfo_children()[1].configure(height=len(tmp_crossings) + 1)
                    elif len(elem.winfo_children()) == 3:
                        elem.winfo_children()[1].configure(height=len(tmp_crossings) + 1) if \
                            self._file_extension == "xyz" else \
                            elem.winfo_children()[2].configure(height=len(tmp_crossings) + 1)
                    elem.configure(height=1 + len(tmp_crossings))
                text.tag_config("spacing1", spacing1=0)
            text.configure(state="disabled")

    def itemise_types_of_lasso(self, elem, n_end, c_end):
        """
            Method reads from an element given as an argument n_end lengths and c_end lengths extended on shallow
            lassos. The values are then put into a list, which structure is equal to
                        [[[deep_n_end_lengths], [all_n_end_lengths],[[deep_c_end_lengths], [all_c_end_lengths]]
        """
        types_of_lassos = [[], []]
        if elem[4] != elem[4 + 3 + int(elem[4]) + 1 + int(elem[5]) + 2]:
            shallow_n_end = []
            if int(elem[4]) is not 0:
                if len(n_end) == 0:
                    types_of_lassos[0] = []
                else:
                    for j in range((8 + int(elem[4])) + 2 + int(elem[5]) + 4,
                                   (8 + int(elem[4])) + 2 + int(elem[5]) + 4 + int(elem[4 + 3 + int(elem[4]) + 1
                                           + int(elem[5]) + 2])):
                        shallow_n_end.append(elem.__getitem__(j) + ",")
                    types_of_lassos[0] = [n_end[:-1], shallow_n_end]
        if elem[5] != elem[5 + 2 + int(elem[4]) + 1 + int(elem[5]) + 3]:
            if len(c_end) == 0:
                types_of_lassos[1] = []
            else:
                shallow_c_end = []
                if int(elem[5]) is not 0:
                    for j in range(9 + int(elem[4]) + int(elem[5]) + 5 + int(
                            elem[9 + int(elem[4]) + int(elem[5]) + 1]) + 1,
                                   9 + int(elem[4]) + int(elem[5]) + 5 + int(
                                       elem[9 + int(elem[4]) + int(elem[5]) + 1]) + 1
                                           + int(elem[9 + int(elem[4]) + int(elem[5]) + 2])):
                        shallow_c_end.append(elem.__getitem__(j) + ",")
                    types_of_lassos[1] = [c_end[:-1], shallow_c_end]
        return types_of_lassos

    def mark_crossings_on_sequence(self):
        piercings = []
        atom = "ca"
        pos_pierc = ""
        neg_pierc = ""
        chain = self.chain_index.get()

        if self.lasinf_smooth_display.get():
            piercings = self.smooth_crossings[self.displayed_lasso]
        else:
            piercings += str(self.array_of_results[self.displayed_lasso][3].get("1.0", "end-1c")). \
                replace("\n", " ").split(" ")
            piercings += str(self.array_of_results[self.displayed_lasso][4].get("1.0", "end-1c")). \
                replace("\n", " ").split(" ")
        piercings = list(filter(len, piercings))

        if len(piercings) != 1:
            for i in piercings:
                if i.__contains__(","):
                    if i.startswith("+"):
                        pos_pierc += "(chain " + chain + " and residue " + str(i[1:-1]) + " and name " + atom + ") "
                        pos_pierc += "(chain " + chain + " and residue " + str(
                            int(i[1:-1]) + 1) + " and name " + atom + ") "
                    elif i.startswith("-"):
                        neg_pierc += "(chain " + chain + " and residue " + str(i[1:-1]) + " and name " + atom + ") "
                        neg_pierc += "(chain " + chain + " and residue " + str(
                            int(i[1:-1]) + 1) + " and name " + atom + ") "
                else:
                    if i.startswith("+"):
                        pos_pierc += "(chain " + chain + " and residue " + str(i[1:]) + " and name " + atom + ") "
                        pos_pierc += "(chain " + chain + " and residue " + str(
                            int(i[1:]) + 1) + " and name " + atom + ") "
                    elif i.startswith("-"):
                        neg_pierc += "(chain " + chain + " and residue " + str(i[1:]) + " and name " + atom + ") "
                        neg_pierc += "(chain " + chain + " and residue " + str(
                            int(i[1:]) + 1) + " and name " + atom + ") "
        else:
            piercings = piercings[0]
            if piercings.endswith(","):
                piercings = piercings[:-1]
            if piercings.startswith("+"):
                pos_pierc += "(chain " + chain + " and residue " + str(piercings[1:]) + " and name " + atom + ") "
                pos_pierc += "(chain " + chain + " and residue " + str(
                    int(piercings[1:]) + 1) + " and name " + atom + ") "
            elif piercings.startswith("-"):
                neg_pierc += "(chain " + chain + " and residue " + str(piercings[1:]) + " and name " + atom + ") "
                neg_pierc += "(chain " + chain + " and residue " + str(
                    int(piercings[1:]) + 1) + " and name " + atom + ") "

        if len(pos_pierc) > 0:
            cmd.select(name="POS_PIERC", selection=pos_pierc[:-1])
            cmd.color(color="lightblue", selection="POS_PIERC")
        if len(neg_pierc) > 0:
            cmd.select(name="NEG_PIERC", selection=neg_pierc[:-1])
            cmd.color(color="palegreen", selection="NEG_PIERC")
        cmd.deselect()

    def get_triangles_coordinates(self, path_to_file):
        if not os.path.isfile(path_to_file):
            self.raise_popup_menu('File with coordinates of vertices not found.')

        input_file = open(path_to_file, "r")
        self.crossing_coord = []
        self.triang_coord = []
        self.shallow_lassos = []
        prevline = ""

        for idx, line in enumerate(input_file.readlines()):
            triangle_coord = list((re.findall(r"(color \$polygon_int\d|-?\d{1,3}[,.]\d{1,2})", line)))
            if len(triangle_coord) == 1:
                coord = list((re.findall(r"-?\d{1,3}[,.]\d{1,2}", prevline)))
                pierc_color = line.split(" ")[-2][:-1]
                if pierc_color == "blue":
                    pierc_color = [0.0, 0.0, 1.0]
                elif pierc_color == "green":
                    pierc_color = [0.0, 1.0, 0.0]
                elif pierc_color == "gray":
                    pierc_color = [0.8, 0.8, 0.8]
                    self.shallow_lassos.append([coord, pierc_color])
                self.crossing_coord.append([coord, pierc_color])
                self.triang_coord.pop(-1)
            elif len(triangle_coord) == 9:
                self.triang_coord.append(triangle_coord)
            prevline = line

    ####################################################################################################################
    #                                    METHODS OPERATING ON OBJECTS IN PYMOL
    ####################################################################################################################

    def pymol_display_chain(self):
        chain = self.chain_index.get()

        if self._file_extension == "pdb":
            cmd.hide(representation="everything", selection="all")
            cmd.delete(name="CHAIN_*")
            cmd.select(name="CHAIN_" + chain, selection="all and chain " + chain)
            cmd.show(representation="cartoon", selection="CHAIN_" + chain)
            cmd.cartoon(type="tube", selection="all")
            cmd.spectrum(palette="rainbow", selection="all")
            cmd.deselect()

    def pymol_view_details(self, chosen_lasso):
        """
            :param chosen_lasso: index of lasso in self.output_data.
        """
        self.delete_pymol_objects()

        if hasattr(self, "lasinf_smooth_button") and self.lasinf_smooth_button.winfo_exists():
            self.lasinf_smooth_button.configure(state="normal")
            self.lasinf_smooth_button.deselect()
        if hasattr(self, "lasinf_surface_button") and self.lasinf_surface_button.winfo_exists():
            self.lasinf_surface_button.configure(state="normal")
            self.lasinf_surface_button.select()

        self.displayed_lasso = chosen_lasso

        self.color_line_in_array_of_results()
        if self.prev_displayed_lasso is not self.displayed_lasso:
            self.decolor_line_in_array_of_results()

        protein = cmd.get_names(type="all")
        if not protein.__contains__(self._filename[:-4]):
            cmd.load(self._filename)

        atom = "ca"
        res_beg = self.output_data[chosen_lasso].split(" ")[1]
        res_end = self.output_data[chosen_lasso].split(" ")[2]
        file_ = self._filename.replace(' ', '')[:-4]

        if not self.is_trajectory:
            chain = self.chain_index.get()
            br_selection = "(chain " + chain + " and residue " + res_beg + " and name " + atom + ")+(chain " + \
                           chain + " and residue " + res_end + " and name " + atom + ")"
            seq_selection = "chain " + chain + " and residue " + res_beg + "-" + res_end
            res_begin = "chain " + chain + " and residue " + res_beg + " and name " + atom
            res_endin = "chain " + chain + " and residue " + res_end + " and name " + atom
        if self._file_extension == "xyz" or not self.is_original_pdb:
            cmd.delete(name=self._filename[:-9] + "*")
            cmd.load(filename=self._full_path_to_file)
            self.connect_xyz_points()
            cmd.show(representation="lines", selection=self._filename[:-9])
            cmd.set(name="line_width", value="4")
            cmd.spectrum(palette="rainbow", selection=self._filename[:-9])
        if self.is_trajectory:
            chain = self.chains[0]
            br_selection = "(chain " + chain + " and residue " + res_beg + " and name " + atom + ")+(chain " + \
                           chain + " and residue " + res_end + " and name " + atom + ")"
            seq_selection = "chain " + chain + " and residue " + res_beg + "-" + res_end
            res_begin = "chain " + chain + " and residue " + res_beg + " and name " + atom
            res_endin = "chain " + chain + " and residue " + res_end + " and name " + atom

            cmd.delete(name=self._filename[:-4] + "*")  #
            cmd.load(filename=self._full_path_to_file)
            self.connect_xyz_points()
            cmd.show(representation="lines", selection=self._filename[:-4])
            cmd.set(name="line_width", value="4")
            cmd.spectrum(palette="rainbow", selection=self._filename[:-4])
        else:
            cmd.select("CHAIN_" + chain, selection="chain " + chain)
            cmd.spectrum(palette="rainbow", selection="CHAIN_" + chain)
            cmd.cartoon(type="tube", selection="CHAIN_" + chain)
            cmd.show(representation="lines", selection="CHAIN_" + chain)
            cmd.spectrum(palette="rainbow", selection="CHAIN_" + chain)

        cmd.select("BR_" + res_beg + "_" + res_end, selection=br_selection)
        cmd.select("SEQ", selection=seq_selection)
        cmd.bond(atom1=res_begin, atom2=res_endin)
        cmd.hide(representation='everything', selection=file_)
        cmd.color(color="gray", selection="SEQ")
        cmd.show(representation='sphere', selection="BR_" + res_beg + "_" + res_end)
        cmd.hide(representation="sticks", selection="SEQ")
        cmd.show(representation="sticks", selection="BR_" + res_beg + "_" + res_end)

        if self._file_extension == "xyz" or self.is_trajectory or not self.is_original_pdb:
            cmd.show(representation="lines", selection=self._filename[:-9])
            if self.is_trajectory and len(self.chains) >= 2:
                self.display_first_chain_in_trajectory()
        else:
            cmd.show(representation="cartoon", selection="CHAIN_" + chain)

        cmd.set(name="sphere_color", value="orange", selection="all")
        cmd.set(name="stick_color", value="orange", selection="all")
        cmd.set("sphere_scale", value=0.5)
        cmd.set(name="seq_view", value=1)
        cmd.set(name="line_width", value=4)
        cmd.deselect()
        print("  Chain, sequence, residues and bridge drawn in PyMOL...")

        file_path = self._filename.replace(".", "_")
        if self.is_trajectory:
            frame = self.output_data[chosen_lasso].split(" ")[0].split("_")[-3]
            if hasattr(self, "given_frames") and len(self.given_frames) is not 0:
                file_with_coord = file_path + os.sep + "frame_" + frame + os.sep + "surface_" + self._filename + "_" + \
                                  chain + "__frame_" + frame + "_" + res_beg + "_" + res_end + ".jms"
                step = int(self.step.getvalue()) if len(self.step.getvalue()) != 0 else 1
                pos_frame = self.retrieved_frames.index(frame)
                cmd.set(name="state", value=(pos_frame+1)*step)
            else:
                file_with_coord = file_path + os.sep + "_surfaces" + os.sep + "surface_" + self._filename + "_" + \
                                  chain + "_lasso_" + res_beg + "_" + res_end + ".jms"
                cmd.set(name="state", value=1)
        else:
            chain = self.chain_index.get()
            file_with_coord = file_path + os.sep + "_surfaces" + os.sep + "surface_" + \
                              self._filename + "_" + chain + "_" + res_beg + "_" + res_end + ".jms"

        self.get_triangles_coordinates(self._full_path_to_dir + os.sep + file_with_coord)

        if None is not self.lasinf_shallow_display and self.lasinf_shallow_display.get():
            self.pymol_draw_triangles(self.shallow_lassos, "SHALLOW_PIERC", 1)
            print("  Shallow lassos drawn in PyMOL...")
        self.pymol_draw_triangles(self.crossing_coord, "PIERC")
        print("  Piercings drawn in PyMOL...")
        self.pymol_draw_surface(self.triang_coord, "TRIANG")
        print("  Surface drawn in PyMOL...")

        cmd.hide(representation="lines", selection="all and sc.")
        cmd.move(axis="z", distance=-130)
        if self._file_extension == "xyz" or self.is_trajectory:
            cmd.orient(selection=self._filename[:-9] + "*")
            cmd.center(selection=self._filename[:-9] + "*")
        else:
            cmd.orient(selection="CHAIN_*")
            cmd.center(selection="CHAIN_*")
        cmd.clip(mode="slab", distance="1000")
        cmd.set(name="two_sided_lighting", value=1)
        if self.is_trajectory:
            if hasattr(self, "given_frames") and len(self.given_frames) is not 0:
                pos_frame = self.retrieved_frames.index(frame)
                self.mark_crossings_on_trajectory(self.retrieved_trajectory_crossings[pos_frame])
            else:
                self.mark_crossings_on_trajectory(self.retrieved_trajectory_crossings[self.displayed_lasso])
        else:
            self.mark_crossings_on_sequence()
        self.prev_displayed_lasso = chosen_lasso

        if hasattr(self, "lasinf_gln_button") and self.lasinf_is_gln_selected.get():
            self.display_gln_objects()

    def pymol_display_shallow_lassos(self):
        if self.lasinf_shallow_display.get():
            for idx, l in enumerate(self.lassos):
                if not len(l[0]) is 0:
                    txt_elem = self.array_of_results[idx][3]
                    self.color_crossings(txt_elem, l[0], 1, idx)
                if not len(l[1]) is 0:
                    txt_elem = self.array_of_results[idx][4]
                    self.color_crossings(txt_elem, l[1], 1, idx)

            if self.displayed_lasso is not None:
                self.pymol_draw_triangles(self.shallow_lassos, "SHALLOW_PIERC", 1)
                cmd.hide(representation="cgo", selection="TRIANG*")
                cmd.hide(representation="cgo", selection="PIERC*")
                cmd.show(representation="cgo", selection="PIERC*")
                cmd.show(representation="cgo", selection="TRIANG*")
            cmd.hide(representation="lines", selection="all and sc.")
            cmd.move(axis="z", distance=-130)

            if self._file_extension == "xyz":
                cmd.orient(selection=self._filename[:-9] + "*")
                cmd.center(selection=self._filename[:-9] + "*")
            else:
                if hasattr(self, "lasinf_smooth_button") and self.lasinf_smooth_button.winfo_exists() \
                        and self.lasinf_smooth_display.get():
                    cmd.orient(selection="SMOOTH_CHAIN_*")
                    cmd.center(selection="SMOOTH_CHAIN_*")
                else:
                    cmd.orient(selection="CHAIN_*")
                    cmd.center(selection="CHAIN_*")

            cmd.clip(mode="slab", distance="1000")
            cmd.set(name="two_sided_lighting", value=1)
            print("  Shallow lassos displayed in PyMOL and in array")
        else:
            for idx, l in enumerate(self.lassos):
                if not len(l[0]) is 0:
                    txt_elem = self.array_of_results[idx][3]
                    self.color_crossings(txt_elem, l[0], 0, idx)
                if not len(l[1]) is 0:
                    txt_elem = self.array_of_results[idx][4]
                    self.color_crossings(txt_elem, l[1], 0, idx)

            if self.displayed_lasso is not None:
                cmd.delete(name="SHALLOW_PIERC*")

    def pymol_display_surface(self):
        if self.displayed_lasso is None:
            self.lasinf_surface_button.select()
            self.raise_popup_menu('No lasso chosen. Please load an appropriate lasso into the PyMOL viewer '
                                  '(e.g. press the button ''view details'') and try again.')

        atom = "ca"
        res_beg = self.output_data[self.displayed_lasso].split(" ")[1]
        res_end = self.output_data[self.displayed_lasso].split(" ")[2]

        if self.is_trajectory:
            chain = self.chains[0]
            atom1 = "residue " + res_beg + " and name " + atom
            atom2 = "residue " + res_end + " and name " + atom
            frame = self.output_data[self.displayed_lasso].split(" ")[0].split("_")[-3]

            if hasattr(self, "given_frames") and len(self.given_frames) is not 0:
                step = int(self.step.getvalue()) if len(self.step.getvalue()) != 0 else 1
                pos_frame = self.retrieved_frames.index(frame)
                cmd.set(name="state", value=(pos_frame + 1) * step)
            else:
                cmd.set(name="state", value=1)
        else:
            chain = self.chain_index.get()
            atom1 = "chain " + chain + " and residue " + res_beg + " and name " + atom
            atom2 = "chain " + chain + " and residue " + res_end + " and name " + atom

        if self.lasinf_surface_display.get():
            cmd.show(representation="cgo", selection="TRIANG")
            cmd.show(representation="cgo", selection="PIERC")
            if None is not self.lasinf_shallow_display and self.lasinf_shallow_display.get():
                cmd.show(representation="cgo", selection="SHALLOW_PIERC")
            cmd.bond(atom1=atom1, atom2=atom2)
            cmd.show(representation="spheres", selection="BR_*")
            cmd.show(representation="sticks", selection="BR_*")
            cmd.set("sphere_scale", value=0.5)
            cmd.set(name="sphere_color", value="orange", selection="all")
            cmd.set(name="stick_color", value="orange", selection="all")
            if not self.is_gln_checkbutton_selected.get() or (hasattr(self, "lasinf_is_gln_selected")
                                                              and not self.lasinf_is_gln_selected.get()):
                if self._file_extension == "xyz":
                    cmd.spectrum(palette="rainbow", selection=self._filename[:-9])
                cmd.spectrum(palette="rainbow", selection="CHAIN_" + chain + "*")
            for name in cmd.get_names(type="all"):
                if name == "SEQ":
                    cmd.color(color="gray", selection="SEQ")
            if not self.is_gln_checkbutton_selected.get() or (hasattr(self, "lasinf_is_gln_selected")
                                                              and not self.lasinf_is_gln_selected.get()):
                if self.is_trajectory:
                    if hasattr(self, "given_frames") and len(self.given_frames) is not 0:
                        pos_frame = self.retrieved_frames.index(frame)
                        self.mark_crossings_on_trajectory(self.retrieved_trajectory_crossings[pos_frame])
                    else:
                        self.mark_crossings_on_trajectory(self.retrieved_trajectory_crossings[self.displayed_lasso])
                else:
                    self.mark_crossings_on_sequence()
        else:
            cmd.hide(representation="cgo", selection="TRIANG")
            cmd.hide(representation="cgo", selection="PIERC")
            if None is not self.lasinf_shallow_display and self.lasinf_shallow_display.get():
                cmd.hide(representation="cgo", selection="SHALLOW_PIERC")
            cmd.unbond(atom1=atom1, atom2=atom2)
            cmd.hide(representation='sphere', selection="BR_*")
            cmd.hide(representation="sticks", selection="BR_*")
            for name in cmd.get_names(type="all"):
                if name == "SEQ":
                    cmd.hide(representation="sticks", selection="SEQ and sc.")
                    cmd.hide(representation="sticks", selection="SEQ")

            if not self.is_gln_checkbutton_selected.get() or (hasattr(self, "lasinf_is_gln_selected")
                                                              and not self.lasinf_is_gln_selected.get()):
                if self.is_trajectory:
                    cmd.spectrum(palette="rainbow", selection=self._filename[:-4])
                elif self._file_extension == "xyz":
                    cmd.spectrum(palette="rainbow", selection=self._filename[:-9])
                if self.lasinf_smooth_display.get():
                    cmd.spectrum(palette="rainbow", selection="SMOOTH_CHAIN_" + chain)
                else:
                    cmd.spectrum(palette="rainbow", selection="CHAIN_" + chain + "*")
            print("  Surface without area of triangulation showed...")

    def pymol_color_structure_gln(self):
        if self.lasinf_is_gln_selected.get():
            if self.displayed_lasso is None:
                self.lasinf_gln_button.deselect()
                self.raise_popup_menu('No lasso chosen. Please load an appropriate lasso into the PyMOL viewer '
                                      '(e.g. press the button ''view details'') and try again.')

            atom = "ca"
            chain = self.chain_index.get()
            self.get_gln_colors()

            for idx, elem in enumerate(self.gln_protein_colors):
                atom_idx = "chain " + chain + " and residue " + str(elem) + " and name " + atom
                atom_idx_1 = "chain " + chain + " and residue " + str(int(elem) + 1) + " and name " + atom
                cmd.set_color(name=atom_idx, rgb=self.gln_protein_colors[elem])
                cmd.color(color=atom_idx, selection=atom_idx)
                cmd.color(color=atom_idx, selection=atom_idx_1)
            print("  Crossings has been colored...")
        else:
            cmd.spectrum(palette="rainbow", selection="CHAIN_*")
            cmd.color(color="gray", selection="SEQ")
            self.mark_crossings_on_sequence()

    def get_gln_colors(self):
        self.gln_protein_colors = {}
        chain = self.chain_index.get()

        file_path = self._filename.replace(".", "_")
        is_smoothed = "_smooth" if self.lasinf_smooth_display.get() else ""

        res_beg = self.output_data[self.displayed_lasso].split(" ")[1]
        res_end = self.output_data[self.displayed_lasso].split(" ")[2]
        filename = self._full_path_to_dir + os.sep + file_path + os.sep + "_surfaces" + os.sep + "surface_" + \
                   self._filename + "_" + chain + "_" + res_beg + "_" + res_end + "_GLN1" + is_smoothed + ".txt"

        with open(filename) as f:
            for elem in f:
                words = list(filter(len, elem.split(" ")))
                self.gln_protein_colors[int(words[3])] = [round(float(words[5]), 3), round(float(words[6]), 3),
                                                          round(float(words[7]), 3)]

    def pymol_display_smooth(self):
        if self.displayed_lasso is None:
            self.lasinf_smooth_button.deselect()
            self.raise_popup_menu('No lasso chosen. Please load an appropriate lasso into the PyMOL viewer '
                                  '(e.g. press the button ''view details'') and try again.')

        if self.lasinf_smooth_display.get():
            self.lasinf_surface_button.select()

            atom = "ca"
            res_beg = self.output_data[self.displayed_lasso].split(" ")[1]
            res_end = self.output_data[self.displayed_lasso].split(" ")[2]

            file_path = self._filename.replace(".", "_")
            if self.is_trajectory:
                chain = self.chains[0]
                seq_selection = "residue " + res_beg + "-" + res_end
                atom1 = "residue " + res_beg + " and name " + atom
                atom2 = "residue " + res_end + " and name " + atom
                br_selection = "(SMOOTH_CHAIN_" + chain + " and residue " + res_beg + " and name " + atom \
                               + ")+(SMOOTH_CHAIN_" + chain + " and residue " + res_end + " and name " + atom + ")"
                frame = self.output_data[self.displayed_lasso].split(" ")[0].split("_")[-3]

                if hasattr(self, "given_frames") and len(self.given_frames) is not 0:
                    file_with_smooth_vert = file_path + os.sep + "frame_" + frame + os.sep + self._filename + "_" + \
                                            chain + "__frame_" + frame + "_" + res_beg + "_" + res_end + "_smooth.pdb"
                    surface_triang_coord = file_path + os.sep + "frame_" + frame + \
                                           os.sep + "surface_" + self._filename + "_" + chain + "__frame_" + frame + \
                                           "_" + res_beg + "_" + res_end + "_smooth.jms"
                    step = int(self.step.getvalue()) if len(self.step.getvalue()) != 0 else 1
                    pos_frame = self.retrieved_frames.index(frame)
                    cmd.set(name="state", value=(pos_frame + 1) * step)
                else:
                    file_with_smooth_vert = file_path + os.sep + "_smooth" + os.sep + self._filename + "_" + chain + \
                                            "_lasso_" + res_beg + "_" + res_end + "_smooth.pdb"
                    surface_triang_coord = file_path + os.sep + "_surfaces" + os.sep + "surface_" + self._filename + \
                                           "_" + self.chains[0] + "_lasso_" + res_beg + "_" + res_end + "_smooth.jms"
                    cmd.set(name="state", value=1)
            else:
                chain = self.chain_index.get()
                seq_selection = "chain " + chain + " and residue " + res_beg + "-" + res_end
                atom1 = "chain " + chain + " and residue " + res_beg + " and name " + atom
                atom2 = "chain " + chain + " and residue " + res_end + " and name " + atom
                br_selection = "(SMOOTH_CHAIN_" + chain + " and chain " + chain + " and residue " + res_beg + \
                               " and name " + atom + ")+(SMOOTH_CHAIN_" + chain + " and chain " + chain + \
                               " and residue " + res_end + " and name " + atom + ")"
                file_with_smooth_vert = file_path + os.sep + "_smooth" + os.sep + self._filename + "_" + chain + "_" \
                                        + res_beg + "_" + res_end + "_smooth.pdb"
                surface_triang_coord = file_path + os.sep + "_surfaces" + os.sep + "surface_" + self._filename + "_" \
                                       + chain + "_" + res_beg + "_" + res_end + "_smooth.jms"

            if not os.path.isfile(self._full_path_to_dir + os.sep + file_with_smooth_vert) and not self.is_artifact:
                self.raise_popup_menu('The smoothed configuration has not been generated. This chain cannot be '
                                      'smoothed while preserving its topology.')
            if not os.path.isfile(self._full_path_to_dir + os.sep + file_with_smooth_vert) and self.is_artifact:
                self.raise_popup_menu('No smooth chain is available when the detected lasso may be artificial.')

            self.delete_pymol_objects()
            cmd.hide(representation="everything", selection="all")
            cmd.load(filename=self._full_path_to_dir + os.sep + file_with_smooth_vert, object="SMOOTH_CHAIN_" + chain)
            cmd.spectrum(palette="rainbow", selection="SMOOTH_CHAIN_" + chain)
            cmd.select("BR_" + res_beg + "_" + res_end, selection=br_selection)
            cmd.bond(atom1=atom1, atom2=atom2)
            cmd.set(name="line_width", value=4)
            cmd.set("sphere_scale", value=0.5)
            cmd.select("SEQ", selection=seq_selection)
            cmd.color(color="gray", selection="SEQ")

            atoms = {'atoms': []}
            atom = "ca"
            cmd.iterate_state(state=1, selection="SMOOTH_CHAIN_" + chain, expression='atoms.append(resi)', space=atoms)
            atoms['atoms'].sort(key=int)
            smooth_atoms = atoms['atoms']

            for idx in range(len(smooth_atoms) - 1):
                cmd.bond(
                    atom1="SMOOTH_CHAIN_" + chain + " and id " + str(smooth_atoms[idx]) + " and name " + atom,
                    atom2="SMOOTH_CHAIN_" + chain + " and id " + str(smooth_atoms[idx + 1]) + " and name " + atom)
            print("  Smoothed chain and bridge drawn in PyMOL...")

            self.get_triangles_coordinates(self._full_path_to_dir + os.sep + surface_triang_coord)

            if hasattr(self, "lasinf_shallow_lasso_button") and self.lasinf_shallow_lasso_button.winfo_exists() \
                    and self.lasinf_shallow_display.get():
                self.pymol_draw_triangles(self.shallow_lassos, "SHALLOW_PIERC", 1)
            self.pymol_draw_triangles(self.crossing_coord, "PIERC")
            self.pymol_draw_surface(self.triang_coord, "TRIANG")
            cmd.show(representation='sphere', selection="BR_" + res_beg + "_" + res_end)
            cmd.show(representation="sticks", selection="BR_" + res_beg + "_" + res_end)
            cmd.set(name="stick_color", value="orange", selection="all")
            cmd.set(name="sphere_color", value="orange", selection="all")
            cmd.move(axis="z", distance=-130)
            cmd.orient(selection="SMOOTH_CHAIN_" + chain)
            cmd.center(selection="SMOOTH_CHAIN_" + chain)
            cmd.clip(mode="slab", distance="1000")
            cmd.set(name="two_sided_lighting", value=1)

            if hasattr(self, "lasinf_gln_button") and self.lasinf_is_gln_selected.get():
                self.display_gln_objects()
            else:
                if self.is_trajectory:
                    if hasattr(self, "given_frames") and len(self.given_frames) is not 0:
                        pos_frame = self.retrieved_frames.index(frame)
                        self.mark_crossings_on_trajectory(self.retrieved_trajectory_crossings[pos_frame])
                    else:
                        self.mark_crossings_on_trajectory(self.retrieved_trajectory_crossings[self.displayed_lasso])
                else:
                    self.mark_crossings_on_sequence()
        else:
            cmd.delete(name="BR_*")
            cmd.set(name="connect_mode", value=0)
            self.pymol_view_details(self.displayed_lasso)
            self.lasinf_surface_button.select()

    def pymol_draw_surface(self, obj_with_coord, pymol_cgo_name):
        triang = [BEGIN, TRIANGLES]

        for vertex in obj_with_coord:
            n = self.compute_normal(float(vertex[0]), float(vertex[1]), float(vertex[2]),
                                    float(vertex[3]), float(vertex[4]), float(vertex[5]),
                                    float(vertex[6]), float(vertex[7]), float(vertex[8]))
            triang.extend([
                COLOR, 0.8, 0.8, 0.8,
                NORMAL, n[0], n[1], n[2],
                VERTEX, float(vertex[0]), float(vertex[1]), float(vertex[2]),
                VERTEX, float(vertex[3]), float(vertex[4]), float(vertex[5]),
                VERTEX, float(vertex[6]), float(vertex[7]), float(vertex[8]),
            ])
        triang.append(END)
        cmd.load_cgo(triang, pymol_cgo_name)

    def compute_normal(self, x1, y1, z1, x2, y2, z2, x3, y3, z3):
        """
            Using coordinates this method computes normal perpendicular to the surface. Used to change the direction
            of light in PyMOL thus the surface of triangles is more visible.
        :return: vector containing normal perpendicular to the surface.
        """
        nx = (y2 - y1) * (z3 - z2) - (z2 - z1) * (y3 - y2)
        ny = (z2 - z1) * (x3 - x2) - (x2 - x1) * (z3 - z2)
        nz = (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2)
        return nx, ny, nz

    def pymol_draw_triangles(self, obj_with_coord, pymol_cgo_name, show_shallow=0):
        """
            Method displays in PyMOL triangles in defined color(blue, green or pink). It creates a cgo object in PyMOL
            containing a finite number of colored triangles.
        :param show_shallow: flag to tell if shallow lassos should be displayed
        :param obj_with_coord: object in class with a list of vertexes
        :param pymol_cgo_name: name under which surface spanned on the loops will be seen
        """
        triang = [BEGIN, TRIANGLES]

        for vertex in obj_with_coord:
            n = self.compute_normal(float(vertex[0][0]), float(vertex[0][1]), float(vertex[0][2]),
                                    float(vertex[0][3]), float(vertex[0][4]), float(vertex[0][5]),
                                    float(vertex[0][6]), float(vertex[0][7]), float(vertex[0][8]), )
            val = [0.2, -0.8, 0.2] if show_shallow == 1 else [0.0, 0.0, 0.0]
            triang.extend([
                COLOR, vertex[1][0] + val[0], vertex[1][1] + val[1], vertex[1][2] + val[2],
                NORMAL, n[0], n[1], n[2],
                VERTEX, float(vertex[0][0]), float(vertex[0][1]), float(vertex[0][2]),
                VERTEX, float(vertex[0][3]), float(vertex[0][4]), float(vertex[0][5]),
                VERTEX, float(vertex[0][6]), float(vertex[0][7]), float(vertex[0][8]),
                       ])
        triang.append(END)
        cmd.load_cgo(triang, pymol_cgo_name)

    @staticmethod
    def set_piercing_color(c):
        if c == "+":
            return "blue"
        elif c == "-":
            return "green"
