import joblib
import pandas as pd
import numpy as np

def preprocess_new_data(df, encoder, label_encoders, feature_names):
    # Start with the columns required for preprocessing
    required_columns = ['CSQ', 'translation', 'mKOZAK_STRENGTH', 'ribo_sorfs_uORFdb', 
                        'uSTOP_CODON', 'uORF_TYPE', 'uKOZAK_STRENGTH', 'uORF_rank', 
                        'uSTART_PHASTCONS', 'uSTART_PHYLOP', '5UTR_LENGTH', 'uORF_count', 
                        'uSTART_mSTART_DIST', 'uORF_LENGTH']
    # Ensure all required columns are in the DataFrame
    for col in required_columns:
        if col not in df.columns:
            df[col] = np.nan  # Assign NaN if the column is missing
    # One-Hot Encode 'CSQ'
    categorical_columns = ['CSQ']
    encoded_features = encoder.transform(df[categorical_columns])
    encoded_df = pd.DataFrame(encoded_features, columns=encoder.get_feature_names_out(categorical_columns))
    df = df.drop(columns=categorical_columns).reset_index(drop=True)
    df = pd.concat([df, encoded_df], axis=1)
    # Label Encode other categorical features
    label_encoded_columns = ['translation', 'mKOZAK_STRENGTH', 'ribo_sorfs_uORFdb', 
                             'uSTOP_CODON', 'uORF_TYPE', 'uKOZAK_STRENGTH']
    for col in label_encoded_columns:
        le = label_encoders[col]
        # Handle unseen labels by assigning a default value or using 'transform' with 'errors' parameter
        df[col] = df[col].map(lambda s: le.transform([s])[0] if s in le.classes_ else -1)
    # Apply the same transformations
    df['uORF_rank'] = df['uORF_rank'].apply(lambda x: 1 if pd.isna(x) else int(str(x).split('_')[0]))
    df['uSTART_PHASTCONS'] = df['uSTART_PHASTCONS'].apply(lambda x: 0.1 if pd.isna(x) else x)
    df['uSTART_PHYLOP'] = df['uSTART_PHYLOP'].apply(lambda x: 0 if pd.isna(x) else x)
    df['uSTOP_CODON'] = df['uSTOP_CODON'].apply(lambda x: 4 if pd.isna(x) else x)
    df['uKOZAK_STRENGTH'] = df['uKOZAK_STRENGTH'].apply(lambda x: 1 if pd.isna(x) else x)
    # Reindex to match training features
    df = df.reindex(columns=feature_names, fill_value=0)
    return df

# Load the model and encoders
rf = joblib.load('random_forest_model.pkl')
encoder = joblib.load('onehot_encoder.pkl')
label_encoders = joblib.load('label_encoders.pkl')
feature_names = joblib.load('feature_names.pkl')

# Load your new data
new_data = pd.read_csv('path_to_your_new_data.csv')
# Preprocess the data
X_new = preprocess_new_data(new_data, label_encoders, feature_names)
# Predict probabilities
y_pred_proba = rf.predict_proba(X_new)[:, 1]
# Attach probabilities to your original data if needed
new_data['5ULTRA_Score'] = y_pred_proba
# Save the results
new_data.to_csv('output.csv', index=False)
