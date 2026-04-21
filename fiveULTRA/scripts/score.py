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
    }).fillna(df.get('uORF_TYPE', ''))
    return df

def score_variants(input_file, output_file, data_dir=os.path.join(os.path.expanduser("~"), ".5ULTRA", "data"), full_anno=False, mane=False):
    rf_model_path = os.path.join(data_dir, 'random_forest_model.pkl')
    encoder_path = os.path.join(data_dir, 'onehot_encoder.pkl')
    rf = joblib.load(rf_model_path)
    encoder = joblib.load(encoder_path)

    # 1. READ AND FIX HEADER (Restores #CHROM, POS, ID, REF, ALT)
    input_df = pd.read_csv(input_file, sep='\t', low_memory=False, comment=None)
    if not isinstance(input_df.index, pd.RangeIndex):
        input_df = input_df.reset_index()

    # Direct positional rename ensures the first 5 columns are always correct
    input_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT'] + list(input_df.columns[5:])
    input_df.columns = [c.strip().lstrip('\ufeff') for c in input_df.columns]

    if input_df.empty:
        logging.warning("No variants found.")
        return False

    # 2. FORCE NUMERIC (Kills the ' NA' crash)
    # We must convert ALL columns that the model uses to avoid string-to-float errors
    for col in ['5UTR_LENGTH', 'uSTART_mSTART_DIST', 'uSTART_PHYLOP', 'uSTART_PHASTCONS', 'pLI', 'LOEUF', 'uORF_count', 'uORF_LENGTH']:
        if col in input_df.columns:
            input_df[col] = pd.to_numeric(input_df[col], errors='coerce')

    input_df['uSTART_CAP_DIST'] = input_df['5UTR_LENGTH'] - input_df['uSTART_mSTART_DIST']

    # --- pLI/LOEUF Handling (Your working logic) ---
    pLI_file = os.path.join(data_dir, "pli_LOEUFByGene.tsv")
    if os.path.exists(pLI_file):
        pLI_data = pd.read_csv(pLI_file, sep="\t")
        pLI_data['pLI'] = pd.to_numeric(pLI_data['pLI'], errors='coerce')
        pLI_data.dropna(subset=["pLI"], inplace=True)
        pLI_data.drop_duplicates(subset="GENE", keep='first', inplace=True)
        pLI_data['GENE'] = pLI_data['GENE'].str.strip().str.upper()
        unique_genes_df = input_df[['GENE']].drop_duplicates()
        unique_genes_df['GENE'] = unique_genes_df['GENE'].str.strip().str.upper()
        gene_pli_loeuF_df = unique_genes_df.merge(pLI_data, on="GENE", how="left")
        input_df = input_df.merge(gene_pli_loeuF_df, on="GENE", how="left", suffixes=('', '_dup'))

    # Standardize names for logic
    rename_mapping = {'ribo_sorfs_uORFdb': 'Ribo_seq', 'translation': 'Translation', 'type': 'Splicing_CSQ'}
    input_df.rename(columns={k: v for k, v in rename_mapping.items() if k in input_df.columns}, inplace=True)

    # 3. SCORING WORKSPACE (Fixes --full)
    # We use proc_df for the model so we don't delete columns from input_df
    columns_model = ['Translation', '5UTR_LENGTH', 'mKOZAK_STRENGTH', 'uORF_count',
        'Ribo_seq', 'uSTART_mSTART_DIST', 'uSTOP_CODON', 'uORF_TYPE',
        'uKOZAK_STRENGTH', 'uORF_LENGTH', 'uORF_rank', 'uSTART_PHYLOP',
        'uSTART_PHASTCONS', 'uSTART_CAP_DIST', 'CSQ', 'GENE', 'pLI', 'LOEUF']

    proc_df = input_df[[c for c in columns_model if c in input_df.columns]].copy()
    proc_df = filter_and_transform(proc_df)

    if not proc_df.empty:
        impute_vals = {'pLI': 0.5, 'LOEUF': 1, 'uSTART_PHYLOP': 0, 'uSTART_PHASTCONS': 0.5}
        for col, val in impute_vals.items():
            if col in proc_df.columns: proc_df[col] = proc_df[col].fillna(val)

        # Encoding
        encoded = pd.DataFrame(encoder.transform(proc_df[['CSQ']]), columns=encoder.get_feature_names_out(['CSQ']), index=proc_df.index)
        proc_df = pd.concat([proc_df.drop(columns=['CSQ', 'GENE'], errors='ignore'), encoded], axis=1)

        mapping = {'Translation': {'increased': 0, 'N-terminal extension': 1, 'decreased': 2}, 'mKOZAK_STRENGTH': {'Weak': 0, 'Adequate': 1, 'Strong': 2}, 'Ribo_seq': {'False':0, 'New uORF':1, 'True':2}, 'uSTOP_CODON': {'TAA': 3, 'TAG': 2, 'TGA': 1}, 'uORF_TYPE': {'N-terminal extension': 1, 'Non-Overlapping': 0, 'Overlapping': 2}, 'uKOZAK_STRENGTH': {'Weak': 0, 'Adequate': 1, 'Strong': 2}}
        for col in mapping:
            if col in proc_df.columns: proc_df[col] = proc_df[col].map(mapping[col]).fillna(0)

        proc_df['uORF_rank'] = proc_df['uORF_rank'].apply(lambda x: 1 if pd.isna(x) else int(str(x).split('_')[0]))
        
        # Predict and map scores back to the original intact input_df
        input_df.loc[proc_df.index, '5ULTRA_Score'] = rf.predict_proba(proc_df[rf.feature_names_in_])[:, 1]

    # 4. FINAL SELECTION (Handles standard and --full)
    if mane and 'MANE' in input_df.columns:
        input_df = input_df[input_df['MANE'].astype(str) != '[]']

    core_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'CSQ', 'Translation', '5ULTRA_Score', 'GENE', 'TRANSCRIPT']
    
    if full_anno:
        # Core columns first, then everything else
        final_list = core_cols + [c for c in input_df.columns if c not in core_cols]
    else:
        final_list = core_cols

    # Filter for columns that exist and save
    existing_cols = [col for col in final_list if col in input_df.columns]
    input_df[existing_cols].to_csv(output_file, sep='\t', index=False)
    
    return True