# Step 1: Import necessary libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Step 2: Load the datasets from GitHub
sift_url = "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/sift.tsv"
foldx_url = "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/datasets/foldX.tsv"
sift = pd.read_csv(sift_url, sep=r'\s+')
foldx = pd.read_csv(foldx_url, sep=r'\s+')

# Step 3: Create 'specific_protein_aa' column
sift["specific_protein_aa"] = sift["Protein"] + "_" + sift["Amino_Acid"]
foldx["specific_protein_aa"] = foldx["Protein"] + "_" + foldx["Amino_Acid"]

# Step 4: Merge both datasets using the 'specific_protein_aa' column
merged = pd.merge(sift, foldx, on="specific_protein_aa", how="inner")
merged = merged.drop(columns=["Protein_y", "Amino_Acid_y"])
merged = merged.rename(columns={"Protein_x": "Protein", "Amino_Acid_x": "Amino_Acid"})

# Print preview of merged data
print("Merged Data:")
print(merged.head())

# Step 5: Extract the first amino acid from the 'Amino_Acid' column
merged['First_Amino_Acid'] = merged['Amino_Acid'].apply(lambda x: x[0])

# Print the updated DataFrame with the first amino acid
print("\nFirst Amino Acid Extracted:")
print(merged[['Amino_Acid', 'First_Amino_Acid']].head())

# Step 6: Create a frequency table of all the first amino acids
first_aa_freq = merged['First_Amino_Acid'].value_counts()

# Step 7: Generate Bar Plot and Pie Chart from the frequency table
# Bar plot
first_aa_freq.plot(kind="bar", color="blue")
plt.xlabel("Amino Acids")
plt.ylabel("Frequency")
plt.title("Frequency of First Amino Acids")
plt.show()

# Pie chart
data = first_aa_freq.values
labels = first_aa_freq.index
colors = sns.color_palette("pastel")[0:len(labels)]
plt.pie(data, labels=labels, colors=colors, autopct="%.1f%%")
plt.title("Distribution of First Amino Acids")
plt.show()

# Step 8: Identify the amino acid with the highest impact on structure and function
# Filtering for the most impactful mutations
most_impact = merged[(merged["sift_Score"] < 0.05) & (merged["foldX_Score"] > 2)]
print("\nMost Impactful Mutations (SIFT < 0.05 and FoldX > 2):")
print(most_impact.head())

# Identifying the most impactful amino acid
impact_aa = most_impact['First_Amino_Acid'].value_counts().idxmax()  # Amino acid with highest impact
print("\nAmino Acid with Highest Impact on Structure and Function:", impact_aa)

# Step 9: Analyze amino acids with more than 100 occurrences in terms of their structural and functional properties
# Filtering for amino acids with more than 100 occurrences in the most impactful mutations
freq_above_100 = most_impact['First_Amino_Acid'].value_counts()
aa_above_100 = freq_above_100[freq_above_100 > 100].index.tolist()

# Analyzing structural and functional properties for these amino acids
for aa in aa_above_100:
    aa_data = most_impact[most_impact['First_Amino_Acid'] == aa]
    avg_sift = aa_data['sift_Score'].mean()
    avg_foldx = aa_data['foldX_Score'].mean()
    print(f"\nAmino Acid: {aa}")
    print(f"Average SIFT Score (Functional Impact): {avg_sift}")
    print(f"Average FoldX Score (Structural Impact): {avg_foldx}")
