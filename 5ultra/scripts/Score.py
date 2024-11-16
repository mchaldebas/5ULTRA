#!/usr/bin/env python3

# Import necessary libraries
from sklearn.impute import SimpleImputer
import sys
import joblib
import pandas as pd

def usage():
    print(f"Usage: {sys.argv[0]} input_file output_file")
    sys.exit(1)

if len(sys.argv) != 3:
    usage()

input_file = sys.argv[1]
output_file = sys.argv[2]

# Load the saved Random Forest model and encoders
rf = joblib.load('./scripts/random_forest_model.pkl')
encoder = joblib.load('./scripts/onehot_encoder.pkl')
label_encoders = joblib.load('./scripts/label_encoders.pkl')

# Read the input data
input_df = pd.read_csv(input_file, sep='\t', low_memory=False)

# Define the function to filter and transform the data
def filter_and_transform(df):
    # Exclude rows where 'CSQ' is 'mKozak'
    df = df[df['CSQ'] != 'mKozak'].copy()
    # Update 'uORF_TYPE' based on 'CSQ' values
    df['uORF_TYPE'] = df['CSQ'].map({
        "uStop_loss to N-terminal extension": "N-terminal extension",
        "uStop_gain to Non-Overlapping": "Non-Overlapping",
        "uStop_loss to Overlapping": "Overlapping"
    }).fillna(df['uORF_TYPE'])
    return df

# Backup the original dataframe before modifications
original_df = input_df.copy()

# Add an identifier column for merging later
original_df['__original_index'] = original_df.index
input_df['__original_index'] = input_df.index

# Continue with transformations
input_df = filter_and_transform(input_df)

# Adding LOEUF and pLI gene annotation
pLI_data = pd.read_csv("./data/pli&LOEUFByGene.tsv", sep="\t").drop_duplicates(subset="GENE")
input_df = input_df.merge(pLI_data, on="GENE", how="left")
input_df['pLI'] = pd.to_numeric(input_df['pLI'], errors='coerce')

# Imputing Missing Values
impute_columns = ['pLI', 'LOEUF', 'uSTART_PHYLOP', 'uSTART_PHASTCONS']
median_imputer = SimpleImputer(strategy='median')
input_df[impute_columns] = median_imputer.fit_transform(input_df[impute_columns])

# Drop unnecessary columns
columns_to_drop = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
                   '5UTR_START', '5UTR_END', 'STRAND', 'GENE', 'TRANSCRIPT', 'mSTART', 
                   'mSTART_CODON', 'START_EXON', 'mKOZAK', 'MANE', 'uORF_START',
                   'uORF_END', 'uSTART_CODON', 'uKOZAK', 'uORF_AA_LENGTH', 'uORF_SEQ', 
                   'Phastcons', 'PhyloP', 'uORF_mSTART_DIST', "Nterminal_count", 
                   "NonOverlapping_count", "Overlapping_count", "5ULTRA_Score"]
input_df = input_df.drop(columns=columns_to_drop, errors='ignore')

# One-Hot Encode 'CSQ'
categorical_columns = ['CSQ']
if 'CSQ' in input_df.columns:
    encoded_features_input = encoder.transform(input_df[categorical_columns])
    encoded_df_input = pd.DataFrame(encoded_features_input, columns=encoder.get_feature_names_out(categorical_columns))
    input_df = input_df.drop(columns=categorical_columns).reset_index(drop=True)
    input_df = pd.concat([input_df, encoded_df_input], axis=1)
else:
    csq_columns = encoder.get_feature_names_out(['CSQ'])
    encoded_df_input = pd.DataFrame(0, index=input_df.index, columns=csq_columns)
    input_df = input_df.reset_index(drop=True)
    input_df = pd.concat([input_df, encoded_df_input], axis=1)

# Label Encode other categorical features using the same mappings as in training
label_encoded_columns = ['translation', 'mKOZAK_STRENGTH', 'ribo_sorfs_uORFdb', 'uSTOP_CODON', 'uORF_TYPE', 'uKOZAK_STRENGTH']
mapping = {
    'translation': {'increased': 0, 'N-terminal extension': 1, 'decreased': 2},
    'mKOZAK_STRENGTH': {'Weak': 0, 'Adequate': 1, 'Strong': 2},
    'ribo_sorfs_uORFdb': {1: 0, 101: 1, 100: 1, 111: 2, 11: 1, 10: 1, 110: 2, 0: 1},
    'uSTOP_CODON': {'TAA': 3, 'TAG': 2, 'TGA': 1, 'TGA > TGA':0, 'TAG > TAG':0, 'TAA > TGA':3, 'TAG > TAA':1, 
                    'TAA > TAA':0, 'TAG > TGA':2, 'TGA > TAA':3, 'TGA > TAG':2, 'TAA > TAG':1},
    'uORF_TYPE': {'N-terminal extension': 1, 'Non-Overlapping': 0, 'Overlapping': 2},
    'uKOZAK_STRENGTH': {'Weak': 0, 'Adequate': 1, 'Strong': 2}
}

for col in label_encoded_columns:
    if col in input_df.columns:
        input_df[col] = input_df[col].map(mapping[col])

# Handling Missing Values and Setting Defaults
input_df['uORF_rank'] = input_df['uORF_rank'].apply(lambda x: 1 if pd.isna(x) else int(str(x).split('_')[0]))
input_df['uSTOP_CODON'] = input_df['uSTOP_CODON'].apply(lambda x: 4 if pd.isna(x) else x)
input_df['uKOZAK_STRENGTH'] = input_df['uKOZAK_STRENGTH'].apply(lambda x: 1 if pd.isna(x) else x)

# Prepare data for prediction
feature_names = rf.feature_names_in_
# Include '__original_index' in columns to retain
input_df = input_df.reindex(columns=list(feature_names) + ['__original_index'], fill_value=0)

# Prepare data for prediction
X_input = input_df[feature_names]  # Use only feature columns
# Predict probabilities
y_pred_proba = rf.predict_proba(X_input)[:, 1]

# Add the predicted probabilities to the processed dataframe
input_df['5ULTRA_Score'] = y_pred_proba

# Merge scores back into the original dataframe
original_df = original_df.merge(
    input_df[['__original_index', '5ULTRA_Score']],
    on='__original_index',
    how='left'
).drop(columns='__original_index')

# Save the results
original_df.to_csv(output_file, sep='\t', index=False)
