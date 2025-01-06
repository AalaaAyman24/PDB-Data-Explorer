---

# PDB Data Explorer 🔬🌐

This project provides a graphical user interface (GUI) for analyzing Protein Data Bank (PDB) files. The tool is built using Python and leverages various bioinformatics libraries to retrieve, parse, and display structural information about proteins, chains, residues, and atoms. Additionally, the tool supports sequence alignment and BLAST search functionalities.

---

## Features

1. **Load PDB Data**: Retrieve and parse PDB files using their unique IDs.
2. **View Details**: Display metadata and header information of the loaded PDB structure.
3. **Chain Analysis**:
   - List chain IDs.
   - Display chain sequences.
   - Count chains.
   - Show chain lengths.
4. **Residue Information**:
   - Display residue names, parent chains, locations, and unique residues.
   - Calculate residue mass and list associated atoms.
5. **Atom Information**:
   - Display atom names, coordinates, element types, and occupancy.
6. **Search Functionality**:
   - Search for specific residues in a PDB file.
   - Perform advanced searches in the PDB database.
7. **BLAST Integration**:
   - Perform BLAST searches using protein sequences.
   - Display best hit information including accession number, e-value, and alignment details.

---

## Installation

### Requirements
- Python 3.8+
- Libraries: tkinter, prody, biopython, numpy, requests

### Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/AalaaAyman24/PDB-Data-Explorer
   ```
   
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
   
3. Run the tool:
   ```bash
   python main.py
   ```

---

## Usage

1. Launch the app.
2. Enter a PDB ID to load data.
3. Use intuitive buttons to explore chains, residues, and atoms or perform BLAST searches.

---

## Acknowledgments

- *PDB*: For structural data.
- *NCBI*: For BLAST services.
- *Biopython*: For simplifying PDB analysis.

---

## License

This project is licensed under the MIT License. Feel free to use, modify, and distribute it as per the license terms.

---

## Contributing

Contributions are welcome! If you'd like to contribute, please follow these steps:
1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Commit your changes and push to your fork.
4. Submit a pull request with a detailed description of your changes.

---

## Contact

For any inquiries or feedback, please contact:
**Aalaa Ayman**
- Email: [aalaasalah389@gmail.com]
- GitHub: [https://github.com/AalaaAyman24]

---
