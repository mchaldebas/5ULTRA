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

    input_df = pd.read_csv(input_file, sep='\t', low_memory=False)

    if input_df.empty:
        logging.warning("No variants found. Exiting.")
        return False

    input_df['5UTR_LENGTH'] = pd.to_numeric(input_df['5UTR_LENGTH'], errors='coerce')
    input_df['uSTART_mSTART_DIST'] = pd.to_numeric(input_df['uSTART_mSTART_DIST'], errors='coerce')
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


    # 1. Create a DataFrame with UNIQUE genes from input_df
    unique_genes_df = input_df[['GENE']].drop_duplicates()
    unique_genes_df['GENE'] = unique_genes_df['GENE'].str.strip().str.upper()


    # 2. Merge the unique genes with pLI_data
    gene_pli_loeuF_df = unique_genes_df.merge(pLI_data, on="GENE", how="left")

    # 3. Handle potential missing pLI/LOEUF (if no match was found)
    if 'pLI' not in gene_pli_loeuF_df.columns or 'LOEUF' not in gene_pli_loeuF_df.columns:
        logging.warning("No matching genes found in pLI data. pLI/LOEUF will be NaN.")
        gene_pli_loeuF_df['pLI'] = float('nan')
        gene_pli_loeuF_df['LOEUF'] = float('nan')
    # Handle cases where all values are NaN.
    elif gene_pli_loeuF_df['pLI'].isnull().all():
        logging.warning("All Pli Values are Nan")

    # 4. Now, merge the complete input_df with the gene_pli_loeuF_df
    input_df = input_df.merge(gene_pli_loeuF_df, on="GENE", how="left")

    # --- Rest of your scoring logic ---
    rename_mapping = {
        'ribo_sorfs_uORFdb': 'Ribo_seq',
        'translation': 'Translation',
        'type': 'Splicing_CSQ'
    }
    input_df.rename(columns={k: v for k, v in rename_mapping.items() if k in input_df.columns}, inplace=True)
    original_df = input_df.copy()

    columns_to_keep = ['Translation', '5UTR_LENGTH', 'mKOZAK_STRENGTH', 'uORF_count',
        'Ribo_seq', 'uSTART_mSTART_DIST', 'uSTOP_CODON', 'uORF_TYPE',
        'uKOZAK_STRENGTH', 'uORF_LENGTH', 'uORF_rank', 'uSTART_PHYLOP',
        'uSTART_PHASTCONS', 'uSTART_CAP_DIST', 'CSQ', 'GENE', 'pLI', 'LOEUF']

    input_df = input_df[columns_to_keep]
    input_df = filter_and_transform(input_df)

    if input_df.empty:
        logging.warning("No variants to score after filtering. Exiting.")
        if mane:
            original_df = original_df[original_df['MANE'] != '[]']
        # ... (rest of your column selection logic) ...
        if full_anno:
            columns_to_remove = ["mSTART", "mSTART_CODON", "minimum_uORF_mSTART_DIST", "uSTART_CODON", "uORF_SEQ"]
            columns_order_to_keep = [
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'CSQ',
            'Translation', '5ULTRA_Score', 'SpliceAI', 'Splicing_CSQ', 'GENE',
            'TRANSCRIPT', 'MANE', '5UTR_START', '5UTR_END', 'STRAND', '5UTR_LENGTH',
            'START_EXON', 'mKOZAK', 'mKOZAK_STRENGTH', 'uORF_count', 'Overlapping_count',
            'Nterminal_count', 'NonOverlapping_count', 'uORF_START', 'uORF_END',
            'Ribo_seq', 'uSTART_mSTART_DIST', 'uSTART_CAP_DIST', 'uSTOP_CODON',
            'uORF_TYPE', 'uKOZAK', 'uKOZAK_STRENGTH', 'uORF_LENGTH', 'uORF_AA_LENGTH',
            'uORF_rank', 'uSTART_PHYLOP', 'uSTART_PHASTCONS', 'pLI', 'LOEUF'
            ]
            columns_to_keep = [col for col in columns_order_to_keep if col in original_df.columns]
            remaining_columns = [col for col in original_df.columns if col not in columns_to_keep]
            columns_to_keep.extend(remaining_columns)
            original_df = original_df[columns_to_keep]
            original_df = original_df.drop(columns=columns_to_remove, errors='ignore')

        else:
            columns_order_to_keep = [
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'CSQ',
            'Translation', '5ULTRA_Score', 'SpliceAI', 'Splicing_CSQ',
            'GENE', 'TRANSCRIPT'
            ]
            columns_to_keep = [col for col in columns_order_to_keep if col in original_df.columns]
            original_df = original_df[columns_to_keep]
        original_df.to_csv(output_file, sep='\t', index=False)
        return True

    impute_columns = ['pLI', 'LOEUF', 'uSTART_PHYLOP', 'uSTART_PHASTCONS']
    median_imputer = SimpleImputer(strategy='median')
    input_df[impute_columns] = median_imputer.fit_transform(input_df[impute_columns])

    columns_to_drop = ['GENE']
    input_df = input_df.drop(columns=columns_to_drop, errors='ignore')

    categorical_columns = ['CSQ']
    if 'CSQ' in input_df.columns:
        encoded_features_input = encoder.transform(input_df[categorical_columns])
        encoded_df_input = pd.DataFrame(
            encoded_features_input,
            columns=encoder.get_feature_names_out(categorical_columns),
            index=input_df.index
        )
        input_df = input_df.drop(columns=categorical_columns)
        input_df = pd.concat([input_df, encoded_df_input], axis=1)

    label_encoded_columns = ['Translation', 'mKOZAK_STRENGTH', 'Ribo_seq', 'uSTOP_CODON', 'uORF_TYPE', 'uKOZAK_STRENGTH']
    mapping = {
        'Translation': {'increased': 0, 'N-terminal extension': 1, 'decreased': 2},
        'mKOZAK_STRENGTH': {'Weak': 0, 'Adequate': 1, 'Strong': 2},
        'Ribo_seq': {'False':0, 'New uORF':1, 'True':2},
        'uSTOP_CODON': {'TAA': 3, 'TAG': 2, 'TGA': 1, 'TGA > TGA':0, 'TAG > TAG':0, 'TAA > TGA':3, 'TAG > TAA':1,
                        'TAA > TAA':0, 'TAG > TGA':2, 'TGA > TAA':3, 'TGA > TAG':2, 'TAA > TAG':1},
        'uORF_TYPE': {'N-terminal extension': 1, 'Non-Overlapping': 0, 'Overlapping': 2},
        'uKOZAK_STRENGTH': {'Weak': 0, 'Adequate': 1, 'Strong': 2}
    }
    for col in label_encoded_columns:
        if col in input_df.columns:
            input_df[col] = input_df[col].map(mapping[col])

    input_df['uORF_rank'] = input_df['uORF_rank'].apply(lambda x: 1 if pd.isna(x) else int(str(x).split('_')[0]))
    input_df['uSTOP_CODON'] = input_df['uSTOP_CODON'].apply(lambda x: 4 if pd.isna(x) else x)
    input_df['uKOZAK_STRENGTH'] = input_df['uKOZAK_STRENGTH'].apply(lambda x: 1 if pd.isna(x) else x)
    input_df['Ribo_seq'] = input_df['Ribo_seq'].apply(lambda x: 0 if pd.isna(x) else x)

    feature_names = rf.feature_names_in_
    X_input = input_df[feature_names]
    y_pred_proba = rf.predict_proba(X_input)[:, 1]
    input_df['5ULTRA_Score'] = y_pred_proba
    ultra_score = input_df[['5ULTRA_Score']]

    original_df = original_df.merge(ultra_score, left_index=True, right_index=True, how="left")
    original_df['Ribo_seq'] = original_df['Ribo_seq'].map({
        1: False, 101: True, 100: True, 111: True, 11: True,
        10: True, 110: True, 0: 'New uORF'
    })

    if mane:
        original_df = original_df[original_df['MANE'] != '[]']

    # ... (rest of your column selection and saving logic) ...
    if full_anno:
            columns_to_remove = ["mSTART", "mSTART_CODON", "minimum_uORF_mSTART_DIST", "uSTART_CODON", "uORF_SEQ"]
            columns_order_to_keep = [
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'CSQ',
            'Translation', '5ULTRA_Score', 'SpliceAI', 'Splicing_CSQ', 'GENE',
            'TRANSCRIPT', 'MANE', '5UTR_START', '5UTR_END', 'STRAND', '5UTR_LENGTH',
            'START_EXON', 'mKOZAK', 'mKOZAK_STRENGTH', 'uORF_count', 'Overlapping_count',
            'Nterminal_count', 'NonOverlapping_count', 'uORF_START', 'uORF_END',
            'Ribo_seq', 'uSTART_mSTART_DIST', 'uSTART_CAP_DIST', 'uSTOP_CODON',
            'uORF_TYPE', 'uKOZAK', 'uKOZAK_STRENGTH', 'uORF_LENGTH', 'uORF_AA_LENGTH',
            'uORF_rank', 'uSTART_PHYLOP', 'uSTART_PHASTCONS', 'pLI', 'LOEUF'
            ]
            columns_to_keep = [col for col in columns_order_to_keep if col in original_df.columns]
            remaining_columns = [col for col in original_df.columns if col not in columns_to_keep]
            columns_to_keep.extend(remaining_columns)
            original_df = original_df[columns_to_keep]
            original_df = original_df.drop(columns=columns_to_remove, errors='ignore')

    else:
        columns_order_to_keep = [
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'CSQ',
            'Translation', '5ULTRA_Score', 'SpliceAI', 'Splicing_CSQ',
            'GENE', 'TRANSCRIPT'
            ]
        columns_to_keep = [col for col in columns_order_to_keep if col in original_df.columns]
        original_df = original_df[columns_to_keep]
    original_df.to_csv(output_file, sep='\t', index=False)
    return True
