# PAD-Genomic_Sequence_Clustering

This repository contains the code for the **PAD 2025 Genomic Sequence Clustering** project.  
The project is implemented in **Python** and structured into four modules (`P1.py`–`P4.py`), each corresponding to one step of the assignment workflow.

---

## 📁 Project Structure

```
PAD-Genomic_Sequence_Clustering/
├─ P1.py          # Sequence file parsing  
├─ P2.py          # Pairwise alignment using Dynamic Programming  
├─ P3.py          # Distance matrix computation  
├─ P4.py          # WPGMA clustering  
├─ environment.yml 
├─ data/ 
└─ README.md  
```


---

## ⚙️ Setup Instructions

This project uses **Anaconda** to manage dependencies.

### 1️⃣ Create the environment from the `.yml` file

If this is your first time setting up the project, run:

conda env create -f environment.yml

### 2️⃣ Activate the environment

conda activate pad2025

To verify the active environment:

conda info --envs


---

## 🧑‍💻 Authors

Project developed as part of **PAD 2025 – Genomic Sequence Clustering** coursework.  
Contributors: *[Your Name(s)]*

---

## 🧾 Notes

- All modules raise `ValueError("malformed input")` for invalid inputs, following the assignment rules.  
- The project is compatible with **Python 3.11**.  
- The file structure and naming must be preserved for submission.  

---