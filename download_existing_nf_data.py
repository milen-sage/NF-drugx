import synapseclient
import pandas as pd

syn = synapseclient.login()

model_map_df = syn.tableQuery("SELECT * FROM syn16979992").asDataFrame().to_csv("model_map.csv")
drugs_preclinical_outcomes_df = syn.tableQuery("SELECT * FROM syn16943049").asDataFrame().to_csv("drug_preclinical_outcomes.csv")
clinical_trials = syn.tableQuery("SELECT * FROM syn16942974").asDataFrame().to_csv("clinical_trials.csv")
model_metadata_df = syn.tableQuery("SELECT * FROM syn16942046").asDataFrame().to_csv("model_metadata.csv")

