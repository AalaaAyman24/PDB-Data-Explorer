import os
import tkinter as tk
from tkinter import *
from prody import *
from Bio.PDB import *
from tkinter import ttk
from tkinter import messagebox
from Bio.PDB import PDBParser, PDBList, Selection
from Bio.PDB.Polypeptide import PPBuilder
import biotite.sequence as biotite_seq
import biotite.sequence.align as align
import biotite.database.rcsb as rcsb
from Bio.Blast import NCBIWWW, NCBIXML


def load_pdb_data():
    pdb_id = (number1Entry.get())
    if len(pdb_id) == 0:
        messagebox.showwarning("Warning", "Invalid entry")
        return None
    pdb.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    file = 'pdb' + pdb_id + '.ent'
    return parser.get_structure(pdb_id, file)

################################ Show Details ########################################


def load_and_show():
    try:
        if len(number1Entry.get()) == 0:
            messagebox.showwarning("Warning", "InValid entry")
            return

        pdb_id = (number1Entry.get())
        pdb.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        file = os.path.join(os.getcwd(), 'pdb' + pdb_id + '.ent')
        data = parser.get_structure(pdb_id, file)

        for i in data.header.keys():
            mylist1.insert(tk.END, f"{i}: {data.header[i]}")

    except Exception as e:
        messagebox.showerror(
            "Error", f"An error occurred while loading or displaying the data: {e}")


################################ Chain ########################################


def chainn():
    try:
        data = load_pdb_data()
        if data is None:
            return

        for chain in data[0]:
            mylist1.insert(tk.END, str("Chain ID: " + str(chain.id)))

    except Exception as e:
        messagebox.showerror(
            "Error", f"An error occurred while fetching chain data: {e}")


def seq():
    try:
        data = load_pdb_data()
        if data is None:
            return

        mylist1.insert(tk.END, "Chain Sequences:")  # Insert header

        for pp in ppb.build_peptides(data):
            mylist1.insert(tk.END, str(pp.get_sequence()), "\n")

    except Exception as e:
        messagebox.showerror(
            "Error", f"An error occurred while fetching sequence data: {e}")


def chain_count():
    try:
        data = load_pdb_data()
        if data is None:
            return

        chain_count = len(list(data[0].get_chains()))
        mylist1.insert(tk.END, f"Total number of chains: {chain_count}")

    except Exception as e:
        messagebox.showerror(
            "Error", f"An error occurred while counting chains: {e}")


def chain_length():
    try:
        data = load_pdb_data()
        if data is None:
            return

        for chain in data[0]:
            length = sum(1 for _ in chain.get_residues())
            mylist1.insert(
                tk.END, f"Chain ID: {chain.id}, Length: {length}")

    except Exception as e:
        messagebox.showerror(
            "Error", f"An error occurred while fetching chain lengths: {e}")


################################ Residue ########################################


def show_residue_info(attribute):
    try:
        data = load_pdb_data()
        if data is None:
            return

        residues = list(data.get_residues())
        for residue in residues:
            if attribute == "name":
                mylist1.insert(tk.END, f"Name: {residue.resname}")
            elif attribute == "parent":
                mylist1.insert(
                    tk.END, f"{residue.resname} its parent is: {residue.get_parent()}")
            elif attribute == "location":
                mylist1.insert(
                    tk.END, f"{residue.resname} its location is: {residue.full_id}")
            elif attribute == "mass":
                mylist1.insert(
                    tk.END, f"{residue.resname} its mass is: {residue.center_of_mass()}")
            elif attribute == "atom":
                mylist1.insert(
                    tk.END, f"{residue.resname} its atoms are: {residue.child_list}")
            elif attribute == "count":
                mylist1.insert(
                    tk.END, f"Total number of residues: {len(residues)}")
                return
            elif attribute == "unique":
                unique_residues = set(res.resname for res in residues)
                mylist1.insert(
                    tk.END, f"Unique residues: {', '.join(unique_residues)}")
                return

            else:
                mylist1.insert(tk.END, f"Unknown attribute: {attribute}")

    except Exception as e:
        messagebox.showerror(
            "Error", f"An error occurred while fetching residue information: {e}")

################################ Atoms ########################################


def show_atom_info(attribute):
    try:
        data = load_pdb_data()
        if data is None:
            messagebox.showwarning("Warning", "No data found")
            return

        atoms = list(data.get_atoms())
        if not atoms:
            messagebox.showwarning(
                "Warning", "No atoms found in the structure")
            return

        if attribute == "name":
            for atom in atoms:
                mylist1.insert(tk.END, f"Name: {atom.get_name()}")

        elif attribute == "mass":
            for atom in atoms:
                mylist1.insert(
                    tk.END, f"{atom.get_name()} its mass is: {atom.mass}")

        elif attribute == "parent":
            for atom in atoms:
                mylist1.insert(
                    tk.END, f"{atom.get_name()} its parent is: {atom.get_parent()}")

        elif attribute == "coord":
            for atom in atoms:
                mylist1.insert(
                    tk.END, f"{atom.get_name()} its coordinate: {atom.get_coord()}")

        elif attribute == "element":
            for atom in atoms:
                mylist1.insert(
                    tk.END, f"{atom.get_name()} element is: {atom.element}")

        elif attribute == "occupancy":
            for atom in atoms:
                mylist1.insert(
                    tk.END, f"{atom.get_name()} its occupancy is: {atom.get_occupancy()}")

        else:
            mylist1.insert(tk.END, f"Unknown attribute: {attribute}")

    except Exception as e:
        messagebox.showerror(
            "Error", f"An error occurred while fetching atom information: {e}")

################################ Search ########################################


def search_in_pdb_file():
    pdb_id = (number1Entry.get())
    if not pdb_id:
        messagebox.showwarning("Input Error", "Please enter a valid PDB ID.")
        return
    try:
        pdb.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        file = 'pdb'+pdb_id+'.ent'
        data = parser.get_structure(pdb_id, file)
    except Exception as e:
        messagebox.showerror(
            "Error", f"Error retrieving or parsing PDB file: {e}")
        return

    amino_acid_name = (number2Entry.get())
    if not amino_acid_name:
        messagebox.showwarning("Input Error", "Please enter a residue name.")
        return

    found_residue = None
    for r in data[0].get_residues():
        if r.resname == str(amino_acid_name):
            found_residue = r
            break

    if found_residue:
        mylist1.insert(400, str(found_residue))
    else:
        mylist1.insert(400, "Residue not found.")


def search_in_pdb_database():
    amino_acid_name = number2Entry.get()
    if not amino_acid_name:
        messagebox.showwarning(
            "Input Error", "Please enter a residue name.")
        return
    try:
        query = rcsb.BasicQuery(amino_acid_name)
        pdb_ids = rcsb.search(query)
        mylist1.insert(400, str(pdb_ids))
        mylist1.insert(400, f"Count: {rcsb.count(query)}")

    except Exception as e:
        messagebox.showerror("Error", f"Error during advanced search: {e}")


# Blast Function ########################################

def blastPDB(sequence):
    try:
        print(f"Submitting sequence to BLAST: {sequence}")
        result_handle = NCBIWWW.qblast("blastp", "pdb", sequence)

        blast_records = NCBIXML.parse(result_handle)
        best_hit = None
        best_alignment = None
        for blast_record in blast_records:
            if blast_record.alignments:
                best_hit = blast_record.alignments[0]
                best_alignment = best_hit.hsps[0]
                break

        if best_hit:
            return {
                "accession": best_hit.accession,
                "e_value": best_alignment.expect,
                "score": best_alignment.score,
                "alignment_length": best_alignment.align_length,
                "identity": best_alignment.identities,
                "query_start": best_alignment.query_start,
                "query_end": best_alignment.query_end
            }
        else:
            return None

    except Exception as e:
        print(f"BLAST error: {e}")
        return None


def Blast():
    try:
        pdb_id = number1Entry.get().strip()
        blast_input = e1.get().strip()

        if not pdb_id or not blast_input:
            messagebox.showerror(
                "Error", "Please enter a valid PDB ID and sequence.")
            return

        print(f"Retrieving PDB file for ID: {pdb_id}")
        pdb = PDBList()
        file_name = pdb.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")

        if not os.path.exists(file_name):
            messagebox.showerror("Error", f"PDB file '{file_name}' not found.")
            return

        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        data = parser.get_structure(pdb_id, file_name)

        print(f"Running BLAST for input sequence: {blast_input}")
        blast_result = blastPDB(blast_input)

        if blast_result is None:
            messagebox.showerror(
                "Error", "BLAST search failed or returned no results.")
            return

        print("Best hit information:", blast_result)

        mylist1.insert(40, "BLAST search completed. Best hit details:")
        mylist1.insert(41, f"Accession: {blast_result['accession']}")
        mylist1.insert(42, f"E-value: {blast_result['e_value']}")
        mylist1.insert(43, f"Score: {blast_result['score']}")
        mylist1.insert(
            44, f"Alignment Length: {blast_result['alignment_length']}")
        mylist1.insert(45, f"Identity: {blast_result['identity']}")
        mylist1.insert(46, f"Query start: {blast_result['query_start']}")
        mylist1.insert(47, f"Query end: {blast_result['query_end']}")

    except Exception as e:
        print(f"An exception occurred: {e}")
        messagebox.showerror("Error", f"An error occurred: {e}")


################################ Align ########################################

def alignn():
    try:
        pdb_id = (number1Entry.get())
    except Exception as e:
        mylist1.insert(END, f"Error getting PDB ID: {e}")
        return

    try:
        pdb.retrieve_pdb_file(pdb_id, pdir=".", file_format="pdb")
    except Exception as e:
        mylist1.insert(END, f"Error retrieving PDB file: {e}")
        return

    try:
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        file = 'pdb' + pdb_id + '.ent'
        data = parser.get_structure(pdb_id, file)
    except Exception as e:
        mylist1.insert(END, f"Error parsing PDB structure: {e}")
        return

    try:
        seq1_input = e2.get().strip()
        seq2_input = e3.get().strip()

        if not seq1_input or not seq2_input:
            mylist1.insert(
                END, "Both sequences must be provided for alignment.")
            return

        seq1 = biotite_seq.ProteinSequence(seq1_input)
        seq2 = biotite_seq.ProteinSequence(seq2_input)
    except Exception as e:
        mylist1.insert(END, f"Error creating ProteinSequence: {e}")
        return

    try:
        matrix = align.SubstitutionMatrix.std_protein_matrix()
    except Exception as e:
        mylist1.insert(END, f"Error creating substitution matrix: {e}")
        return

    try:
        mylist1.insert(END, "Local Alignment")
        alignments = align.align_optimal(seq1, seq2, matrix, local=True)
        for ali in alignments:
            ali_str = str(ali)
            for line in ali_str.splitlines():
                mylist1.insert(END, line)
    except Exception as e:
        mylist1.insert(END, f"Error with local alignment: {e}")
        return

    mylist1.insert(END, "\n")

    try:
        mylist1.insert(END, "Global Alignment")
        alignments = align.align_optimal(seq1, seq2, matrix, local=False)
        for ali in alignments:
            ali_str = str(ali)
            for line in ali_str.splitlines():
                mylist1.insert(END, line)
            mylist1.insert(END, "\n")
    except Exception as e:
        mylist1.insert(END, f"Error with global alignment: {e}")


################################ Clear ########################################


def delete():
    mylist1.delete(0, END)


################################ Download ########################################


def download_pdb_file():
    num1 = (number1Entry.get())
    pdb.retrieve_pdb_file(num1, pdir=".", file_format="pdb")
    messagebox.showinfo(" Information", "File downloaded successfully")

################################ GUI ########################################


ppb = PPBuilder()
pdb = PDBList()
root = tk.Tk()
root.title("PDB Database")
root.geometry("1000x700")
root.minsize(800, 600)


Font_tuple = ("Times New Roman", 10, "bold")
Font_tuple1 = ("Times New Roman", 12, "bold")
Font_tuple2 = ("Arial", 12)


############################ sidebar #############################


sidebar = tk.Frame(root, bg="#133E87", width=200, height=700)
sidebar.pack_propagate(False)
sidebar.pack(side=LEFT, fill=Y)


btn1 = Button(sidebar, text="Protein Details", bg="#608BC1",
              fg="white", font=Font_tuple1, command=load_and_show)
btn1.pack(fill=X, padx=15, pady=10)


### Chain ###
btn2 = Menubutton(sidebar, text="Chain", bg="#608BC1",
                  fg="white", font=Font_tuple1)
btn2.pack(fill=X, padx=15, pady=10)

btn2.menu = Menu(btn2, tearoff=0)
btn2["menu"] = btn2.menu
btn2.menu.add_separator()
btn2.menu.add_command(label="Display chains", command=chainn)
btn2.menu.add_separator()
btn2.menu.add_command(label="Sequence of chain", command=seq)
btn2.menu.add_separator()
btn2.menu.add_command(label="Display chains", command=chain_count)
btn2.menu.add_separator()
btn2.menu.add_command(label="Sequence of chain", command=chain_length)
btn2.menu.add_separator()

### Residues ###
btn3 = Menubutton(sidebar, text="Residues", bg="#608BC1",
                  fg="white", font=Font_tuple1)
btn3.pack(fill=X, padx=15, pady=10)

btn3.menu = Menu(btn3, tearoff=0)
btn3["menu"] = btn3.menu
btn3.menu.add_separator()
btn3.menu.add_command(label="Residue Name",
                      command=lambda: show_residue_info("name"))
btn3.menu.add_separator()
btn3.menu.add_command(label="Residue Parent",
                      command=lambda: show_residue_info("parent"))
btn3.menu.add_separator()
btn3.menu.add_command(label="Residue Mass",
                      command=lambda: show_residue_info("mass"))
btn3.menu.add_separator()
btn3.menu.add_command(label="Residue Atoms",
                      command=lambda: show_residue_info("atom"))
btn3.menu.add_separator()
btn3.menu.add_command(label="Residue Location",
                      command=lambda: show_residue_info("location"))
btn3.menu.add_separator()
btn3.menu.add_command(label="Residue Count",
                      command=lambda: show_residue_info("count"))
btn3.menu.add_separator()
btn3.menu.add_command(label="Unique Residue",
                      command=lambda: show_residue_info("unique"))
btn3.menu.add_separator()

### Atoms ###
btn4 = Menubutton(sidebar, text="Atoms", bg="#608BC1",
                  fg="white", font=Font_tuple1)
btn4.pack(fill=X, padx=15, pady=10)

btn4.menu = Menu(btn4, tearoff=0)
btn4["menu"] = btn4.menu
btn4.menu.add_separator()
btn4.menu.add_command(
    label="Atom Name", command=lambda: show_atom_info("name"))
btn4.menu.add_separator()
btn4.menu.add_command(label="Atoms Mass",
                      command=lambda: show_atom_info("mass"))
btn4.menu.add_separator()
btn4.menu.add_command(label="Atoms Parent",
                      command=lambda: show_atom_info("parent"))
btn4.menu.add_separator()
btn4.menu.add_command(label="Atoms Coordinate",
                      command=lambda: show_atom_info("coord"))
btn4.menu.add_separator()
btn4.menu.add_command(label="Atoms Element",
                      command=lambda: show_atom_info("element"))
btn4.menu.add_separator()
btn4.menu.add_command(label="Atoms Occubancy",
                      command=lambda: show_atom_info("occupancy"))
btn4.menu.add_separator()

### Search File ###
btn5 = Button(sidebar, text="Search in File", bg="#608BC1",
              fg="white", font=Font_tuple1, command=search_in_pdb_file)
btn5.pack(fill=X, padx=15, pady=10)


### Search Database ###
btn6 = Button(sidebar, text="Search in Database", bg="#608BC1",
              fg="white", font=Font_tuple1, command=search_in_pdb_database)
btn6.pack(fill=X, padx=15, pady=10)

### Blast ###
btn7 = Button(sidebar, text="Blast Sequence", bg="#608BC1",
              fg="white", font=Font_tuple1, command=Blast)
btn7.pack(fill=X, padx=15, pady=10)


### Align ###
btn8 = Button(sidebar, text="Align Sequence", bg="#608BC1",
              fg="white", font=Font_tuple1, command=alignn)
btn8.pack(fill=X, padx=15, pady=10)


### Clear ###
btn9 = Button(sidebar, text="Clear All", bg="#608BC1",
              fg="white", font=Font_tuple1, command=delete)
btn9.pack(fill=X, padx=15, pady=10)


### Download ###
btn10 = Button(sidebar, text="Download File", bg="#608BC1",
               fg="white", font=Font_tuple1, command=download_pdb_file)
btn10.pack(fill=X, padx=15, pady=10)


frame_main = tk.Frame(root, bg="#f0f0f0", padx=20, pady=20)
frame_main.pack(side=RIGHT, fill=BOTH, expand=True)


### Exit ###
btn11 = Button(sidebar, text="Exit", bg="#608BC1",
               fg="white", font=Font_tuple1, command=frame_main.quit)
btn11.pack(fill=X, padx=15, pady=10)


############################ Info #############################

number1label = Label(
    frame_main, text="Enter The Code of Protein", fg="darkblue", font=Font_tuple2)
number1label.pack(anchor="w", pady=5)

number1Entry = Entry(frame_main, font=Font_tuple2)
number1Entry.pack(fill=X, pady=10)

############################ Search  #############################

number2label = Label(
    frame_main, text="Search For Specific Molecule", fg="darkblue", font=Font_tuple2)
number2label.pack(anchor="w", pady=5)
number2Entry = Entry(frame_main, font=Font_tuple2)
number2Entry.pack(fill=X, pady=10)


############################ blast #############################

e1_label = Label(frame_main, text="Blast Sequence In Protein",
                 fg="darkblue", font=Font_tuple2)

e1_label.pack(anchor="w", pady=5)
e1 = Entry(frame_main, font=('calibre', 15, 'normal'))
e1.pack(side=TOP, fill=X, padx=5, pady=5)


############################ align #############################

e2_label = Label(frame_main, text="Enter The First Sequence",
                 fg="darkblue", font=Font_tuple2)

e2_label.pack(anchor="w", pady=5)
e2 = Entry(frame_main, font=('calibre', 15, 'normal'))
e2.pack(side=TOP, fill=X, padx=5, pady=5)


e3_label = Label(frame_main, text="Enter The Second Sequence",
                 fg="darkblue", font=Font_tuple2)
e3_label.pack(anchor="w", pady=5)
e3 = Entry(frame_main, font=('calibre', 15, 'normal'))
e3.pack(side=TOP, fill=X, padx=5, pady=5)


############################ Results #############################


results_label = Label(frame_main, text="Results Output",
                      fg="darkblue", font=Font_tuple2)
results_label.pack(anchor="w", pady=5)


################################ Results Output Box #############################


mylist1 = Listbox(frame_main, font='Arial 10 bold',
                  fg='navy blue', width=50, height=20)
mylist1.pack(fill=BOTH, expand=True, padx=10, pady=10)

scroll_bar = Scrollbar(frame_main, command=mylist1.yview)
scroll_bar.pack(side=RIGHT, fill=Y)
mylist1.config(yscrollcommand=scroll_bar.set)


############################## Start GUI ##########################################

root.mainloop()
