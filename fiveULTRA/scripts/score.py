import os
import joblib
import pandas as pd
from sklearn.impute import SimpleImputer
import logging

def filter_and_transform(df):
    """
    Filters and transforms the input DataFrame.
    """
    df = df[df['CSQ'] != 'mKozak'].copy()
    df['uORF_TYPE'] = df['CSQ'].map({
        "uStop_loss to N-terminal extension": "N-terminal extension",
        "uStop_gain to Non-Overlapping": "Non-Overlapping",
        "uStop_loss to Overlapping": "Overlapping"
    }).fillna(df['uORF_TYPE'])
    return df

def score_variants(input_file, output_file, data_dir=os.path.join(os.path.expanduser("~"), ".5ULTRA", "data"), full_anno=False, mane=False):
    """
    Scores variants and handles pLI/LOEUF merging correctly.
    """
    rf_model_path = os.path.join(os.path.expanduser(data_dir), 'random_forest_model.pkl')
    encoder_path = os.path.join(os.path.expanduser(data_dir), 'onehot_encoder.pkl')
    rf = joblib.load(rf_model_path)
    encoder = joblib.load(encoder_path)

    input_df = pd.read_csv(input_file, sep='\t', low_memory=False, comment=None)
    if not isinstance(input_df.index, pd.RangeIndex):
        input_df = input_df.reset_index()

    # Rename the first 5 columns by position to guarantee they are #CHROM, POS, ID, REF, ALT
    vcf_cols = {input_df.columns[0]: '#CHROM', input_df.columns[1]: 'POS', 
                input_df.columns[2]: 'ID', input_df.columns[3]: 'REF', 
                input_df.columns[4]: 'ALT'}
    input_df.rename(columns=vcf_cols, inplace=True)

    # Clean the rescued column names
    input_df.columns = [str(c).replace('level_0', '#CHROM').replace('level_1', 'POS').replace('level_2', 'ID').replace('level_3', 'REF').replace('level_4', 'ALT') for c in input_df.columns]

    if input_df.empty:
        logging.warning("No variants found. Exiting.")
        return False

    for col in ['5UTR_LENGTH', 'uSTART_mSTART_DIST', 'uSTART_PHYLOP', 'uSTART_PHASTCONS']:
        if col in input_df.columns:
            input_df[col] = pd.to_numeric(input_df[col], errors='coerce')

    input_df['uSTART_CAP_DIST'] = input_df['5UTR_LENGTH'] - input_df['uSTART_mSTART_DIST']

    # --- pLI/LOEUF Handling ---
    pLI_file = os.path.join(os.path.expanduser(data_dir), "pli_LOEUFByGene.tsv")
    if not os.path.exists(pLI_file):
        logging.error(f"pLI file not found: {pLI_file}")
        return False

    pLI_data = pd.read_csv(pLI_file, sep="\t")
    pLI_data['pLI'] = pd.to_numeric(pLI_data['pLI'], errors='coerce')
    pLI_data.dropna(subset=["pLI"], inplace=True)
    pLI_data.drop_duplicates(subset="GENE", keep='first', inplace=True)
    pLI_data['GENE'] = pLI_data['GENE'].str.strip().str.upper()
    pLI_data = pLI_data.dropna(subset=['GENE'])

    unique_genes_df = input_df[['GENE']].drop_duplicates()
    unique_genes_df['GENE'] = unique_genes_df['GENE'].str.strip().str.upper()
    gene_pli_loeuF_df = unique_genes_df.merge(pLI_data, on="GENE", how="left")

    if 'pLI' not in gene_pli_loeuF_df.columns or 'LOEUF' not in gene_pli_loeuF_df.columns:
        gene_pli_loeuF_df['pLI'] = float('nan')
        gene_pli_loeuF_df['LOEUF'] = float('nan')

    input_df = input_df.merge(gene_pli_loeuF_df, on="GENE", how="left")

    rename_mapping = {'ribo_sorfs_uORFdb': 'Ribo_seq', 'translation': 'Translation', 'type': 'Splicing_CSQ'}
    input_df.rename(columns={k: v for k, v in rename_mapping.items() if k in input_df.columns}, inplace=True)
    original_df = input_df.copy()

    columns_to_keep = ['Translation', '5UTR_LENGTH', 'mKOZAK_STRENGTH', 'uORF_count',
        'Ribo_seq', 'uSTART_mSTART_DIST', 'uSTOP_CODON', 'uORF_TYPE',
        'uKOZAK_STRENGTH', 'uORF_LENGTH', 'uORF_rank', 'uSTART_PHYLOP',
        'uSTART_PHASTCONS', 'uSTART_CAP_DIST', 'CSQ', 'GENE', 'pLI', 'LOEUF']

    input_df = input_df[columns_to_keep]
    input_df = filter_and_transform(input_df)

    if input_df.empty:
        # Final output selection (ensure correct names are used)
        cols_final = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'CSQ', 'Translation', '5ULTRA_Score', 'GENE', 'TRANSCRIPT']
        columns_to_keep = [col for col in cols_final if col in original_df.columns]
        original_df[columns_to_keep].to_csv(output_file, sep='\t', index=False)
        return True

    impute_columns = ['pLI', 'LOEUF', 'uSTART_PHYLOP', 'uSTART_PHASTCONS']
    specific_impute_values = {'pLI': 0.5, 'LOEUF': 1, 'uSTART_PHYLOP': 0, 'uSTART_PHASTCONS': 0.5}
    
    for col in impute_columns:
        input_df[col] = input_df[col].fillna(specific_impute_values[col])

    # Model Encoding and Prediction
    categorical_columns = ['CSQ']
    if 'CSQ' in input_df.columns:
        encoded_features_input = encoder.transform(input_df[categorical_columns])
        encoded_df_input = pd.DataFrame(encoded_features_input, columns=encoder.get_feature_names_out(categorical_columns), index=input_df.index)
        input_df = pd.concat([input_df.drop(columns=['CSQ', 'GENE'], errors='ignore'), encoded_df_input], axis=1)

    mapping = {'Translation': {'increased': 0, 'N-terminal extension': 1, 'decreased': 2}, 'mKOZAK_STRENGTH': {'Weak': 0, 'Adequate': 1, 'Strong': 2}, 'Ribo_seq': {'False':0, 'New uORF':1, 'True':2}, 'uSTOP_CODON': {'TAA': 3, 'TAG': 2, 'TGA': 1}, 'uORF_TYPE': {'N-terminal extension': 1, 'Non-Overlapping': 0, 'Overlapping': 2}, 'uKOZAK_STRENGTH': {'Weak': 0, 'Adequate': 1, 'Strong': 2}}
    for col in mapping:
        if col in input_df.columns: input_df[col] = input_df[col].map(mapping[col]).fillna(0)

    input_df['uORF_rank'] = input_df['uORF_rank'].apply(lambda x: 1 if pd.isna(x) else int(str(x).split('_')[0]))
    input_df['5ULTRA_Score'] = rf.predict_proba(input_df[rf.feature_names_in_])[:, 1]
    
    # Merge Score back
    original_df = original_df.merge(input_df[['5ULTRA_Score']], left_index=True, right_index=True, how="left")

    # Output selection
    cols_order = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'CSQ', 'Translation', '5ULTRA_Score', 'GENE', 'TRANSCRIPT']
    columns_to_keep = [col for col in cols_order if col in original_df.columns]
    
    original_df[columns_to_keep].to_csv(output_file, sep='\t', index=False)
    return True