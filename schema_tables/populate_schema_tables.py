import pandas as pd
import numpy as np
import difflib
import sqlalchemy as sa
import os

def fuzzy_match(term, target_list):  

    match = difflib.get_close_matches(term, target_list, n = 1)
    
    if not match:
        return True

    return False


def print_table(df, name):

    print("TABLE " + name)
    #print(df.to_string())
    #print("END")

    print("=========================================================================")
    print("=========================================================================")

    return


def get_model_type(x):

    model_type = np.NaN
    if pd.isnull(x["isInVivoModel"]) and pd.isnull(x["isPDX"]) and pd.isnull(x["isCellLine"]):
        model_type = np.NaN
    elif not pd.isnull(x["isPrimary"]) and x["isPrimary"]:
        model_type = "primary cell line"
    elif not pd.isnull(x["isCellLine"]) and x["isCellLine"]:
        model_type = "cell line" 
    elif (not pd.isnull(x["isInVivoModel"]) and x["isInVivoModel"]) and (pd.isnull(x["isPDX"]) or not x["isPDX"]):
        model_type = "in vivo"
    elif not pd.isnull(x["isPDX"]) and x["isPDX"]:
        model_type = "PDX" 

    return model_type



synapse_data_path = "./synapse_data"

if __name__ == "__main__":

    # read gene targets data from drug-target explorer DB
    gene_targets = pd.read_csv(os.path.join(synapse_data_path, "drug_target_associations.csv"), sep = " ")
    
    # read model metadata
    model_meta = pd.read_csv(os.path.join(synapse_data_path, "model_metadata_v1.csv"))
    # read drug preclinical outcomes data
    preclinical_outcomes_data = pd.read_csv(os.path.join(synapse_data_path, "./drug_preclinical_outcomes.csv"))

    # get relevant model data columns from preclinical data
    preclinical_models = preclinical_outcomes_data[["model_name", "organism", "model_type"]].drop_duplicates(subset = "model_name")
    # change column names for consistency
    preclinical_models = preclinical_models.rename(columns = {"model_name":"modelSystemName", "model_type":"type"})

    
    # merge model data from preclinical and meta model data; on model system name; use only models defined in meta models data (i.e. left join)
    model_meta = model_meta.merge(preclinical_models, how = "left", on = "modelSystemName")

    # get organism value for each model either from preclinical models table or from meta model data table
    model_meta["organism"] = model_meta.apply(lambda x: x["organism_y"] if not x["organism_x"] else x["organism_x"], axis = 1)
    model_meta.drop(["organism_x", "organism_y"], axis = 1, inplace = True)

    # determine if a preclinical data model has already been included in models metadata; allow for fuzzy match so that capitalization and spaces don't affect results 
    # (note that the fuzzy match may not be the best way to do this, especially if different model names differ only by + or - sign (e.g. indicating genotype)
    # we can fine tune the match later using the cutoff argument (not using pandas join on model name here because only exact matches would be considered))
    preclinical_models["new"] = preclinical_models["modelSystemName"].apply(fuzzy_match, args = (model_meta["modelSystemName"].values,))

    # prune unnamed models and only include new models
    preclinical_models = preclinical_models[(preclinical_models["modelSystemName"] != "unnamed") & (preclinical_models["new"] == True)]

    '''
    merge preclinical models with meta models data
    '''

    # don't have model type information in preclinical models, so add NaN values
    model_nans = np.empty(len(preclinical_models.index))
    model_nans[:]  = np.nan
    preclinical_models["type"] = model_nans

    # concatenate preclinical and meta model data, at relevant columns
    # remove Synapse provenance information for now
    model_meta = pd.concat([model_meta, preclinical_models], axis = 0, ignore_index = True, sort = True).drop(["ROW_ID", "ROW_VERSION"], axis = 1)

 
    
    #####################
    #process CELL LINE
    #####################

    cell_line_table = model_meta[["modelSystemName", "cellosaurusId", "isPrimary", "isCellLine", "ATCC"]]
    cell_line_table = cell_line_table.loc[cell_line_table["isCellLine"] != False].reset_index()
    cell_line_table = cell_line_table.drop(["isCellLine"], axis = 1)

    # rename columns to match DB schema
    cell_line_table.rename(columns = {"cellosaurusId":"cellosaurus_id", "modelSystemName":"name", "index":"model_id"}, inplace = True)

    # set cell line type (e.g.non-primary vs primary)
    cell_line_table["type"] = cell_line_table["isPrimary"].apply(lambda x: np.NaN if pd.isnull(x) else ("primary" if x else "secondary"))

    # ensure a unique set of cell line models
    cell_line_table = cell_line_table.drop_duplicates(subset = "name")

    print("Processed cell line")


    
    #####################
    #process ORGANISM
    #####################
    
    # get data relevant for organism table
    organism_table = model_meta[["organism"]].dropna(axis = "index").drop_duplicates().reset_index()

    # rename columns to match DB schema
    organism_table = organism_table[["organism"]].rename(columns = {"organism":"name"})

    print("Processed organism")


    #####################
    #process  DISEASE
    #####################

    # get data relevant for disease table; namely, unique disease names
    disease_table = model_meta[["disease"]].dropna(axis = "index").drop_duplicates().reset_index()
    disease_table.drop("index", axis = 1, inplace = True)

    # rename columns to match DB schema
    disease_table = disease_table.rename(columns = {"disease":"name"})

    # add EFO_id for future compatibility with OpenTargets
    disease_table_nans = np.empty(len(disease_table.index)) 
    disease_table_nans[:] = np.nan
    disease_table["efo_id"] = disease_table_nans 
    
    # add normal as a disease free state to the disease table
    normal_row = pd.DataFrame([["normal", 'nan']], columns = ["name", "efo_id"])
    disease_table = disease_table.append(normal_row, ignore_index = True)

    print("Processed disease")
   
   
   
    #####################
    #process SYMPTOM
    #####################

    # get data relevant for symptom table; namely, unique symptom names
    symptom_table = model_meta[["modelOf", "nfSymptomModel"]].dropna(axis = "index")

    # remove control tussue models and leave only tumor symptoms (unique)
    symptom_table = symptom_table[symptom_table["nfSymptomModel"] != "control"]
    symptom_table = symptom_table["modelOf"].drop_duplicates().reset_index()
    symptom_table.drop("index", axis = 1, inplace = True)

    # rename columns to match DB schema
    symptom_table = symptom_table.rename(columns = {"modelOf":"name"})

    # add EFO_id for future compatibility with OpenTargets
    symptom_table_nans = np.empty(len(symptom_table.index)) 
    symptom_table_nans[:] = np.nan
    symptom_table["efo_id"] = symptom_table_nans

    # add normal as a disease free symptom state to the symptom table
    normal_row = pd.DataFrame([["normal", "nan"]], columns = ["name", "efo_id"])
    symptom_table = symptom_table.append(normal_row, ignore_index = True)

    print("Processed symptom")
    


    #####################
    #process MODEL
    #####################
    
    # get relevant data for model table
    model_table = model_meta[["modelSystemName", "isInVivoModel", "isPDX", "organism", "modelOf", "disease"]].reset_index()
    
    # rename columns to match DB schema
    model_table.rename(columns = {"modelSystemName":"name", "index":"model_id"}, inplace = True)

    model_table["type"] = model_meta.apply(get_model_type, axis = 1) 
    
    
    # get cell line ids from cell line table that have matching meta model data id
    model_table = pd.merge(model_table, cell_line_table[["model_id"]].reset_index(), how = "left", on = "model_id")
 
    # pandas reset_index returns an index column (in this case the cell line id) w/o an option for rename; so do that manually to match DB 
    model_table.rename(columns = {"index":"cell_line_id"}, inplace = True)

    # get organism ids by joining on organism names
    model_table = pd.merge(model_table, organism_table.reset_index()[["name", "index"]], how = "left", left_on = "organism", right_on = "name")

    # pandas reset_index returns an index column (in this case the organism id) w/o an option for rename; so do that manually to match DB schema 
    # remove extraneous columns brought in by join as well as correct automatic pandas join name changes (e.g. name to name_x, etc.)
    model_table = model_table.rename(columns = {"index":"organism_id", "name_x":"name"}).drop(["name_y"], axis = 1)


    # get symptom ids by joining on symptom names
    model_table = pd.merge(model_table, symptom_table.reset_index()[["name", "index"]], how = "left", left_on = "modelOf", right_on = "name")

    # pandas reset_index returns an index column (in this case the symptom id) w/o an option for rename; so do that manually to match DB 
    # remove extraneous columns brought in by join as well as correct automatic pandas join name changes (e.g. name to name_x, etc.)
    model_table = model_table.rename(columns = {"index":"symptom_id", "name_x":"name"}).drop(["name_y", "modelOf"], axis = 1)

    # get disease ids by joining on disease names
    model_table = pd.merge(model_table, disease_table.reset_index()[["name", "index"]], how = "left", left_on = "disease", right_on = "name")

    # pandas reset_index returns an index column (in this case the disease id) w/o an option for rename; so do that manually to match DB 
    # remove extraneous columns brought in by join as well as correct automatic pandas join name changes (e.g. name to name_x, etc.)
    model_table = model_table.rename(columns = {"index":"disease_id", "name_x":"name"}).drop(["name_y", "disease", "isPDX", "isInVivoModel"], axis = 1)


    # though current data doesn't have that instance, multiple model systems can be
    # used for the same symptom or disease; and vice versa, multiple symptoms and
    # diseases can be studied in the same model system; hence we create many-to-many
    # relation tables

    # populate disease model table
    disease_model_table = model_table.reset_index()[["disease_id", "index"]]
    disease_model_table.rename(columns = {"index":"model_id"}, inplace = True)
    disease_model_table.dropna(inplace = True)
    
    # populate symptom model table
    symptom_model_table = model_table.reset_index()[["symptom_id", "index"]]
    symptom_model_table.rename(columns = {"index":"model_id"}, inplace = True)
    symptom_model_table.dropna(inplace = True)

    # drop symptom and disease id from model table
    model_table.drop(["disease_id", "symptom_id"], axis = 1, inplace = True)
    
    print("Processed model")


    
    #####################
    #process GENE
    #####################
    
    # get data relevant for gene table from model meta data
    gene_table = model_meta[["genePerturbed"]].dropna(axis = "index").reset_index()

    # get genes from drug target explorer
    gene_names = gene_targets[["hugo_gene"]].drop_duplicates(subset = ["hugo_gene"])

    # rename columns to match DB schema
    gene_table = gene_table.rename(columns = {"genePerturbed":"name", "index":"model_id"})
    gene_names = gene_names.rename(columns = {"hugo_gene":"name"})
    # don't have model information for genes in the drug target explorer
    gene_nans = np.empty(len(gene_names.index))
    gene_nans[:] = np.nan
    gene_names["model_id"] = gene_nans

    # aggregate all genes
    gene_table = pd.concat([gene_table, gene_names], axis = 0, ignore_index = True, sort = True)

    # add ensembl_id for future compatibility with OpenTargetsi; for now set to nan - look up ensembl_ids based on HUGO later
    gene_table_nans = np.empty(len(gene_table.index))
    gene_table_nans[:] = np.nan
    gene_table["ensembl_id"] = gene_table_nans


    print("Processed gene")
 

    #####################
    #process  DYSREGULATED GENE
    #####################

    # could be done more efficiently, together with processing gene table, but for now let's keep them separated

    # get data relevant for dysregulated genes, selecting only records with associated gene IDs (cannot have dysregulated gene w/o a gene ID)
    dysregulated_gene_table = model_meta[["genePerturbationType", "genePerturbationMethod"]].loc[gene_table["model_id"]].reset_index()

    # rename columns to match DB schema
    dysregulated_gene_table = dysregulated_gene_table.rename(columns = {"genePerturbationType":"gene_perturbation_type", "genePerturbationMethod":"gene_perturbation_method", "index":"model_id"}) 
 
    # join gene id information
    dysregulated_gene_table = pd.merge(gene_table[["model_id"]].dropna().reset_index(), dysregulated_gene_table, how = "left", on = "model_id")

    # rename gene table index column to match DB schema
    dysregulated_gene_table.rename(columns = {"index":"gene_id"}, inplace = True)

    # add transcript perturbation and allele frequency columns to populate as we get more data
    dysregulated_gene_nans = np.empty(len(dysregulated_gene_table.index))
    dysregulated_gene_nans[:] = np.nan
    dysregulated_gene_table["transcript_perturbation"] = dysregulated_gene_nans 
    dysregulated_gene_table["allele_frequency"] = dysregulated_gene_nans 

 
    print("Processed dysregulated gene")


    #####################
    #process  MODEL:DISREGULATED GENE
    #####################
    
    # get model_id and dysregulated_gene_id from the dysregulated_gene_table 
    # and store them in a new junction table; however, note that in general this relationship need not hold 
    # after the DB is instantiated, i.e. the relationship between model and dysregulated gene is many-to-many
    # there could be multiple models involving the same dysregulated gene and multiple dysregulated genes 
    # captured in a model

    # select relevant data from dysregulated_gene_table
    model_dysregulated_gene_table = dysregulated_gene_table.reset_index()
    model_dysregulated_gene_table = model_dysregulated_gene_table[["index", "model_id"]]

    # rename columns to match DB schema
    model_dysregulated_gene_table.rename(columns = {"index":"dysregulated_gene_id"}, inplace = True)

 


    '''
    # drop table publication for now

    #####################
    #process PUBLICATION
    #####################

    # get relevant publication data from preclinical outcomes data; dedup publication id's
    publication_table = preclinical_outcomes_data[["pub"]].drop_duplicates(subset = ["pub"]).dropna()
    
    # publications are indexed by synapse_ids; rename columns to follow DB Schema
    publication_table = publication_table.rename(columns = {"pub":"synapse_id"}).reset_index(drop = True)


    # TODO: don't have publication names for now; can query synapse later
    # for now add empty name column
    publication_nans = np.empty(len(publication_table.index))
    publication_nans[:] = np.nan
    publication_table["name"] = publication_nans


    print("Processed publication")
    '''
    
    #####################
    #process FUNDER
    #####################

    # get funder information from meta model data; for now this only includes NTAP...
    model_meta["funder"] = model_meta["ntapFunded"].apply(lambda x: "NTAP" if x else np.NaN)
    
    funder_table = model_meta[["funder"]].drop_duplicates(subset = ["funder"]).dropna().reset_index(drop = True)
    funder_table.rename(columns = {"funder":"name"}, inplace = True)


    print("Processed funder")
    
    #####################
    #process STUDY
    #####################

    # get columns relevant to study from preclinical data
    study_table = preclinical_outcomes_data[["pmid", "pub"]].drop_duplicates(subset = ["pmid"]) 

    # rename columns to match DB schema
    study_table.rename(columns= {"pub":"synapse_id"}, inplace = True)


    '''
    # drop publication table for now
    # get publication ids from publication table
    study_table = study_table.merge(publication_table.reset_index()[["synapse_id", "index"]], how = "right", on = "synapse_id")
    # rename columns to match DB schema
    study_table.rename(columns = {"index":"publication_id"}, inplace = True)
    ''' 
    # get funder name from model meta data using model pmid
    # convert modelPmid to string (since values in study table can be non numeric
    model_meta["modelPmid"] = model_meta["modelPmid"].astype(str)
    study_table = study_table.merge(model_meta[["modelPmid", "funder"]], how = "left", left_on = "pmid", right_on = "modelPmid")

    # get funder ids from funder table
    study_table = study_table.merge(funder_table.reset_index()[["name", "index"]], how = "left", left_on = "funder", right_on = "name")

    # rename columns to match DB schema
    study_table.rename(columns = {"index":"funder_id", "name":"description"}, inplace = True)

    # clean unnecessary fields to match schema (keep pmid for now)
    study_table = study_table[["synapse_id", "funder_id", "pmid"]]


    print("Processed study")

    #####################
    #process EXPERIMENTAL INTERVENTION
    #####################

    # get different intervention types from preclinical experiments
    intervention_table = preclinical_outcomes_data[["intervention_type", "std_intervention_name"]].reset_index(drop = True)

    # get other intervention types from trials
    intervention_trials = pd.read_csv(os.path.join(synapse_data_path, "nf_trials_standardized.csv"))[["std_intervention_name", "std_intervention_type"]].reset_index(drop = True)
    intervention_trials.rename(columns = {"std_intervention_type":"intervention_type"}, inplace = True)

    # aggregate all different interventions
    intervention_table = pd.concat([intervention_trials, intervention_table], axis = 0, ignore_index = True, sort = True)

    # remove all duplicate interventions
    intervention_table.drop_duplicates(subset = ["std_intervention_name"], inplace = True)

    # rename columns to match DB
    intervention_table.rename(columns = {"std_intervention_name":"name", "intervention_type":"type"}, inplace = True)
    intervention_nans = np.empty(len(intervention_table.index))
    intervention_nans[:] = np.nan
    intervention_table["description"] = intervention_nans


    print("Processed experimental intervention")

    #####################
    #process EFFICACY
    #####################

    # get different efficacy levels
    efficacy_table = preclinical_outcomes_data[["efficacy"]].dropna()

    # rename columns to match DB schema
    efficacy_table.rename(columns = {"efficacy":"code"}, inplace = True)
    # remove duplicate efficacy types
    efficacy_table = efficacy_table.drop_duplicates(subset = ["code"]).reset_index(drop = True)

    # we don't have efficacy description at this point; so populate with nan
    efficacy_nans = np.empty(len(efficacy_table.index))
    efficacy_nans[:] = np.nan
    efficacy_table["description"] = efficacy_nans
    efficacy_table["type"] = efficacy_nans


    print("Processed efficacy")
    
    #####################
    #process DRUG
    #####################
    
    print("Processing drugs")

    # read drug target explorer table
    drug_table = pd.read_csv(os.path.join(synapse_data_path, "drugs.csv"))#.head(1000)

    # rename columns to match DB schema
    drug_table.rename(columns = {"internal_id":"DT_explorer_internal_id", "std_name":"name"}, inplace = True)


    print("Processed drug")


    #####################
    #process DRUG SYNONYMS
    #####################

    print("Processing drug synonyms")

    # read drug target explorer table
    drug_synonyms_table = pd.read_csv(os.path.join(synapse_data_path, "drug_synonyms_dte.csv"))#.head(1000)

    # rename columns to match DB schema
    drug_synonyms_table.rename(columns = {"internal_id":"DT_explorer_internal_id", "common_name":"name"}, inplace = True)

    # get drug ids from drug table
    drug_synonyms_table = drug_synonyms_table.merge(drug_table.reset_index()[["DT_explorer_internal_id", "index"]], how = "left", on = "DT_explorer_internal_id")
    # rename drug id to match DB
    drug_synonyms_table.rename(columns = {"index":"drug_id"}, inplace = True)


    print("Processed drug synonyms")
    
     
    
    #####################
    #process DRUG EXPERIMENT
    #####################
    
    # first, create a drug screen table containing individual drug-assay responses or summary drug assay responses
    # second, create a drug assay table containing individual drug dosages within a given drug assay part of a given drug screen
    # or containing a range of dosage over a set of drug assays with a given drug screen
    
    print("Processing drug screen and drug screen assays")

    # get drug screens from all preclinical experiments; all but the oncolytic virus based interventions are drugs
    drug_experiment_table = preclinical_outcomes_data[["outcome", "outcome_unit", "outcome_type", "dosing", "std_intervention_name", "intervention_type"]]
    drug_experiment_table = preclinical_outcomes_data[preclinical_outcomes_data["intervention_type"] != "virus"].reset_index()

    # rename columns to match DB
    drug_experiment_table.rename(columns = {"index":"experiment_id", "outcome":"response", "outcome_unit":"response_unit", "outcome_type":"response_type", "dosing":"dosage_note"}, inplace = True)

    # standardize drug names a bit more to account for difference in capitalization and white spaces
    name_std = drug_table["name"].str.lower()
    name_std = name_std.str.strip()
    drug_table["name"] = name_std

    name_std = drug_experiment_table["std_intervention_name"].str.lower()
    name_std = name_std.str.strip()
    drug_experiment_table["std_intervention_name"] = name_std    
    
    # get drug ids based on drug names
    drug_experiment_table = drug_experiment_table.merge(drug_table[["name"]].reset_index(), how = "left", left_on = "std_intervention_name", right_on = "name")

    # rename drug id to match DB schema
    drug_experiment_table.rename(columns = {"index":"drug_id"}, inplace = True)

    # drop unnecessary columns for the drug_screen table
    drug_screen_table = drug_experiment_table[["experiment_id", "response", "response_unit", "response_type"]]    

    # drop unnecessary columns for the drug_screen_assay table
    drug_assay_table = drug_experiment_table[["dosage_note", "drug_id"]]    
    
    # add columns to match schema though we do not populate with data we don't have at this point
    # (we will populate these columns as new screening data is available for ingestion following the drug screen ingestion SOP) 
    drug_assay_fill = np.empty(len(drug_assay_table.index))
    drug_assay_fill[:] = np.NaN
    drug_assay_table["dosage"] = drug_assay_fill
    drug_assay_table["dosage_unit"] = drug_assay_fill
    drug_assay_table["type"] = "summary"
    drug_assay_table["drug_assay_id"] = np.NaN
    drug_assay_fill[:] = np.arange(0, len(drug_assay_fill)) 
    drug_assay_table["drug_screen_id"] = drug_assay_fill
    drug_assay_table["drug_screen_assay_id"] = drug_assay_table.apply(lambda x: str(int(x["drug_assay_id"])) + "_" + str(int(x["drug_screen_id"])) + "_" + x["type"] if not pd.isnull(x["drug_assay_id"]) else str(int(x["drug_screen_id"])) + "_"  + x["type"], axis = 1)
    drug_assay_table.drop(["drug_screen_id"], axis = 1, inplace = True)
    
    # copy over drug_screen_assay_id to drug_screen table (i.e. a foreign key)
    drug_screen_table["drug_screen_assay_id"] = drug_assay_table["drug_screen_assay_id"]

    print("Processed drug experiment")

    
    
    #####################
    #process EXPERIMENT
    #####################

    # this is the main juncture table relating information about diseases, symptoms, models, study, study properties (e.g. cohort size, efficacy, intervention type, etc.)
    intervention_experiment_table = model_meta.merge(preclinical_outcomes_data, how = "right", left_on = "modelSystemName", right_on = "model_name")

    
    # get study id 
    intervention_experiment_table = intervention_experiment_table.merge(study_table[["pmid"]].reset_index(), how = "left", on = "pmid")
    # rename study index to match DB schema
    intervention_experiment_table.rename(columns = {"index":"study_id"}, inplace = True)
    # remove pmid from study table since it's no longer needed
    study_table.drop(["pmid"], axis = 1, inplace = True)


    # get model id
    intervention_experiment_table = intervention_experiment_table.merge(model_table[["name"]].reset_index(), how = "left", left_on = "model_name", right_on = "name") 
    # rename model index to match DB schema
    intervention_experiment_table.rename(columns = {"index":"model_id"}, inplace = True)


    # get intervention id
    intervention_experiment_table = intervention_experiment_table.merge(intervention_table[["name"]].reset_index(), how = "left", left_on = "std_intervention_name", right_on = "name") 
    # rename model index to match DB schema
    intervention_experiment_table.rename(columns = {"index":"intervention_id"}, inplace = True)
    

    # get efficacy type id
    intervention_experiment_table = intervention_experiment_table.merge(efficacy_table[["code"]].reset_index(), how = "left", left_on = "efficacy", right_on = "code") 
    # rename model index to match DB schema
    intervention_experiment_table.rename(columns = {"index":"efficacy_id"}, inplace = True)


    # leave only columns specified in DB schema
    intervention_experiment_table = intervention_experiment_table[["intervention_id", "efficacy_id", "model_id", "study_id", "cohort_size"]]

    print("Processed experiment")



    #####################
    #process DRUG:TARGET_GENE
    #####################
 
    print("Processing drug target")
    # get relevant information from drug-target explorer associations
    drug_target_gene_table = gene_targets[["internal_id", "hugo_gene", "n_quantitative", "cv", "sd", "IC50_nM", "AC50_nM", "EC50_nM", "Potency_nM", "Ki_nM", "Kd_nM", "n_qualitative", "total_n", "confidence", "pchembl_d", "pchembl_t", "known_selectivity_index"]]

    # get gene ids from gene table based on hugo symbol
    drug_target_gene_table = drug_target_gene_table.merge(gene_table[["name"]].reset_index(), how = "left", left_on = "hugo_gene", right_on = "name")
    # rename gene index to match DB schema
    drug_target_gene_table.rename(columns = {"index":"gene_id"}, inplace = True)
  
    # get drug ids from drug table based on drug target explorer internal id
    drug_target_gene_table = drug_target_gene_table.merge(drug_table[["DT_explorer_internal_id"]].reset_index(), how = "left", left_on = "internal_id", right_on = "DT_explorer_internal_id")
    # rename drug index to match DB schema
    drug_target_gene_table.rename(columns = {"index":"drug_id"}, inplace = True)
    drug_target_gene_table = drug_target_gene_table[["drug_id", "gene_id", "n_quantitative", "cv", "IC50_nM", "AC50_nM", "EC50_nM", "Potency_nM", "Ki_nM", "Kd_nM", "n_qualitative", "total_n", "confidence", "pchembl_d", "pchembl_t", "known_selectivity_index"]]


    print("Processed drug target")
    
    # clean tables from extraneous fields to match DB schema and store as csv's matching DB schema


    print("######################################")
    print("Storing tables...")
    print("######################################")

    cell_line_table = cell_line_table.drop(["model_id"], axis = 1).reset_index().rename(columns = {"index":"id"})
    print_table(cell_line_table, "cell_line")
    #df_to_db("cell_line", cell_line_table)
    cell_line_table.to_csv("cell_line.csv", index = False)

    model_table = model_table.drop(["model_id"], axis = 1).reset_index().rename(columns = {"index":"id"})
    print_table(model_table, "model")
    #df_to_db("model", model_table)
    model_table.to_csv("model.csv", index = False)
    
    gene_table = gene_table.drop(["model_id"], axis = 1).reset_index().rename(columns = {"index":"id"})
    print_table(gene_table, "gene")
    #df_to_db("gene", gene_table)
    gene_table.to_csv("gene.csv", index = False)
    
    dysregulated_gene_table = dysregulated_gene_table.drop(["model_id"], axis = 1).reset_index().rename(columns = {"index":"id"})
    print_table(dysregulated_gene_table, "dysregulated_gene")
    #df_to_db("dysregulated_gene", dysregulated_gene_table)
    dysregulated_gene_table.to_csv("dysregulated_gene.csv", index = False)
    
    model_dysregulated_gene_table = model_dysregulated_gene_table.reset_index().rename(columns = {"index":"id"})
    print_table(model_dysregulated_gene_table, "model_dysregulated_gene")
    #df_to_db("model_dysregulated_gene", model_dysregulated_gene_table)
    model_dysregulated_gene_table.to_csv("model_dysregulated_gene.csv", index = False)
   
    organism_table = organism_table.reset_index().rename(columns = {"index":"id"}) 
    print_table(organism_table, "organism")
    #df_to_db("organism", organism_table)
    organism_table.to_csv("organism.csv", index = False)

    disease_table = disease_table.reset_index().rename(columns = {"index":"id"})
    print_table(disease_table, "disease")
    #df_to_db("disease", disease_table)
    disease_table.to_csv("disease.csv", index = False)

    disease_model_table = disease_model_table.reset_index().rename(columns = {"index":"id"})
    print_table(disease_model_table, "disease_model")
    #df_to_db("disease", disease_table)
    disease_model_table.to_csv("disease_model.csv", index = False)

    symptom_table = symptom_table.reset_index().rename(columns = {"index":"id"})
    print_table(symptom_table, "symptom")
    #df_to_db("symptom", symptom_table)
    symptom_table.to_csv("symptom.csv", index = False)
    
    symptom_model_table = symptom_model_table.reset_index().rename(columns = {"index":"id"})
    print_table(symptom_model_table, "symptom_model")
    #df_to_db("disease", disease_table)
    symptom_model_table.to_csv("symptom_model.csv", index = False)

    '''
    # drop the publication table for now
    publication_table = publication_table.reset_index().rename(columns = {"index":"id"})
    print_table(publication_table, "publication")
    #df_to_db("publication", publication_table)
    publication_table.to_csv("publication.csv", index = False)
    ''' 

    funder_table = funder_table.reset_index().rename(columns = {"index":"id"})
    print_table(funder_table, "funder")
    #df_to_db("funder", funder_table)
    funder_table.to_csv("funder.csv", index = False)

    study_table = study_table.reset_index().rename(columns = {"index":"id"})
    print_table(study_table, "study")
    #df_to_db("study", study_table)
    study_table.to_csv("study.csv", index = False)
    
    intervention_table = intervention_table.reset_index().rename(columns = {"index":"id"})
    print_table(intervention_table, "intervention")
    #df_to_db("experimental_intervention", intervention_table)
    intervention_table.to_csv("experimental_intervention.csv", index = False)
    
    efficacy_table = efficacy_table.reset_index().rename(columns = {"index":"id"}) 
    print_table(efficacy_table, "efficacy")
    #df_to_db("efficacy", efficacy_table)
    efficacy_table.to_csv("efficacy.csv", index = False)
    
    drug_table = drug_table.reset_index().rename(columns = {"index":"id"})
    print_table(drug_table.head(200), "drug")
    print("Drug records: " + str(len(drug_table.index)))
    #df_to_db("drug", drug_table)
    drug_table.to_csv("drug.csv", index = False)
    
    drug_synonyms_table = drug_synonyms_table.reset_index().rename(columns = {"index":"id"})
    print_table(drug_synonyms_table.head(200), "drug synonyms")
    #df_to_db("drug_synonyms", drug_synonyms_table)
    drug_synonyms_table.to_csv("drug_synonyms.csv", index = False)
    
    drug_screen_table = drug_screen_table.reset_index().rename(columns = {"index":"id"}) 
    print_table(drug_screen_table, "drug screen")
    #df_to_db("drug_screen", drug_screen_table)
    drug_screen_table.to_csv("drug_screen.csv", index = False)
    
    drug_assay_table = drug_assay_table.reset_index().rename(columns = {"index":"id"}) 
    print_table(drug_assay_table, "drug screen assay")
    #df_to_db("drug_screen_assay", drug_assay_table)
    drug_assay_table.to_csv("drug_screen_assay.csv", index = False)
    
    drug_target_gene_table = drug_target_gene_table.reset_index().rename(columns = {"index":"id"})
    print_table(drug_target_gene_table.head(200), "drug target-gene")
    #df_to_db("drug_target_gene", drug_target_gene_table)
    drug_target_gene_table.to_csv("drug_target_gene.csv", index = False)

    intervention_experiment_table = intervention_experiment_table.reset_index().rename(columns = {"index":"id"})
    print_table(intervention_experiment_table, "experiment")
    #df_to_db("experiment", intervention_experiment_table)
    intervention_experiment_table.to_csv("experiment.csv", index = False)
