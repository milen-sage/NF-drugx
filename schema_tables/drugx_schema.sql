CREATE TABLE cell_line
(
id INTEGER NOT NULL UNIQUE ,
name TEXT NOT NULL,
cellosaurus_id TEXT,
is_primary BOOLEAN,
PRIMARY KEY (id)
);

CREATE TABLE disease
(
id INTEGER NOT NULL UNIQUE ,
name TEXT NOT NULL,
efo_id INT,
PRIMARY KEY (id)
);

CREATE TABLE drug
(
id INTEGER NOT NULL UNIQUE ,
name TEXT NOT NULL,
DT_explorer_internal_id INTEGER,
inchikey TEXT,
inchi TEXT,
std_smiles TEXT,
PRIMARY KEY (id)
);

CREATE TABLE drug_synonyms
(
id INTEGER NOT NULL UNIQUE ,
drug_id INTEGER,
name TEXT NOT NULL,
DT_explorer_internal_id INT,
PRIMARY KEY (id)
);

CREATE TABLE efficacy
(
id INTEGER NOT NULL UNIQUE ,
description TEXT,
type TEXT,
PRIMARY KEY (id)
);

CREATE TABLE experimental_intervention
(
id INTEGER NOT NULL UNIQUE ,
description TEXT,
type TEXT,
PRIMARY KEY (id)
);

CREATE TABLE drug_experiment
(
id INTEGER NOT NULL UNIQUE ,
experiment_id INTEGER NOT NULL,
drug_id INTEGER,
dosage TEXT,
response TEXT,
response_unit TEXT,
response_type TEXT,
PRIMARY KEY (id)
);

CREATE TABLE funder
(
id INTEGER NOT NULL UNIQUE ,
name TEXT NOT NULL,
PRIMARY KEY (id)
);

CREATE TABLE gene
(
id INTEGER NOT NULL UNIQUE ,
name TEXT,
ensembl_id INT,
PRIMARY KEY (id)
);

CREATE TABLE dysregulated_gene
(
id INTEGER NOT NULL UNIQUE ,
gene_id INTEGER,
gene_perturbation_method TEXT,
transcript_perturbation TEXT,
allele_frequency TEXT,
gene_perturbation_type TEXT,
PRIMARY KEY (id)
);

CREATE TABLE drug_target_gene
(
id INTEGER NOT NULL UNIQUE ,
drug_id INTEGER,
gene_id INTEGER,
n_quantitative INTEGER,
cv FLOAT,
IC50_nM FLOAT,
AC50_nM FLOAT,
EC50_nM FLOAT,
potency_nM FLOAT,
Ki_nM FLOAT,
Kd_nM FLOAT,
n_qualitative FLOAT,
total_n INT,
confidence FLOAT,
pchembl_d FLOAT,
known_selectivity_index FLOAT,
pchembl_t FLOAT,
PRIMARY KEY (id)
);

CREATE TABLE model_disregulated_gene
(
id INTEGER NOT NULL UNIQUE ,
model_id INTEGER,
dysregulated_gene_id INTEGER,
PRIMARY KEY (id)
);

CREATE TABLE organism
(
id INTEGER NOT NULL UNIQUE ,
name TEXT NOT NULL,
PRIMARY KEY (id)
);

CREATE TABLE model
(
id INTEGER NOT NULL UNIQUE ,
name TEXT NOT NULL,
cell_line_id INTEGER,
organism_id INTEGER,
is_in_vivo_model BOOLEAN,
is_pdx BOOLEAN,
PRIMARY KEY (id)
);

CREATE TABLE publication
(
id INTEGER NOT NULL UNIQUE ,
name TEXT NOT NULL,
synapse_id TEXT,
PRIMARY KEY (id)
);

CREATE TABLE study
(
id INTEGER NOT NULL UNIQUE ,
description TEXT NOT NULL,
publication_id INTEGER,
pmid TEXT,
funder_id INTEGER,
PRIMARY KEY (id)
);

CREATE TABLE symptom
(
id INTEGER NOT NULL UNIQUE ,
name TEXT NOT NULL,
efo_id INT,
PRIMARY KEY (id)
);

CREATE TABLE experiment
(
id INTEGER NOT NULL UNIQUE ,
intervention_id INTEGER NOT NULL,
efficacy_id INTEGER,
disease_id INTEGER,
symptom_id INTEGER,
model_id INTEGER,
cohort_size INTEGER,
study_id INTEGER,
PRIMARY KEY (id)
);

/* READ CSV TABLES */
LOAD DATA LOCAL INFILE './cell_line.csv' INTO TABLE cell_line FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @name, @cellosaurus_id, @is_primary) 
SET name = nullif(@name,''), 
cellosaurus_id = nullif(@cellosaurus_id,''), 
is_primary = @is_primary = 'True';

LOAD DATA LOCAL INFILE './organism.csv' INTO TABLE organism FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @name)
SET name = nullif(@name,''); 

LOAD DATA LOCAL INFILE './model.csv' INTO TABLE model FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @name, @is_in_vivo_model, @is_pdx, @cell_line_id, @organism_id)
SET name = nullif(@name,''), 
is_in_vivo_model = @is_in_vivo_model = 'True',
is_pdx = @is_pdx = 'True',
cell_line_id = nullif(@cell_line_id,''),
organism_id = nullif(@organism_id,'');

LOAD DATA LOCAL INFILE './disease.csv' INTO TABLE disease FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id,	@name, @efo_id)
SET name = nullif(@name,''), 
efo_id = nullif(@efo_id,'');

LOAD DATA LOCAL INFILE './symptom.csv' INTO TABLE symptom FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id,	@name, @efo_id)
SET name = nullif(@name,''), 
efo_id = nullif(@efo_id,'');

LOAD DATA LOCAL INFILE './study.csv' INTO TABLE study FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @pmid, @publication_id, @description, @funder_id)
SET pmid = nullif(@pmid,''),
publication_id = nullif(@publication_id,''),
description = nullif(@description,''),
funder_id = nullif(@funder_id,'');

LOAD DATA LOCAL INFILE './model_dysregulated_gene.csv' INTO TABLE model_disregulated_gene FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @dysregulated_gene_id, @model_id)
SET dysregulated_gene_id = nullif(@dysregulated_gene_id,''),
model_id = nullif(@model_id,'');

LOAD DATA LOCAL INFILE './experiment.csv' INTO TABLE experiment FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @intervention_id, @efficacy_id, @disease_id, @symptom_id, @model_id, @study_id, @cohort_size)
SET intervention_id = nullif(@intervention_id,''),
efficacy_id = nullif(@efficacy_id,''),
disease_id = nullif(@disease_id,''),
symptom_id = nullif(@symptom_id,''),
model_id = nullif(@model_id,''),
study_id = nullif(@study_id,''),
cohort_size = nullif(@cohort_size,'');

LOAD DATA LOCAL INFILE './publication.csv' INTO TABLE publication FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @synapse_id, @name)
SET synapse_id = nullif(@synapse_id,''),
name = nullif(@name,'');

LOAD DATA LOCAL INFILE './funder.csv' INTO TABLE funder FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @name)
SET name = nullif(@name,'');

LOAD DATA LOCAL INFILE './dysregulated_gene.csv' INTO TABLE dysregulated_gene FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @gene_id, @gene_perturbation_type, @gene_perturbation_method, @transcript_perturbation, @allele_frequency)
SET gene_id = nullif(@gene_id,''),
gene_perturbation_type = nullif(@gene_perturbation_type,''),
gene_perturbation_method = nullif(@gene_perturbation_method,''),
transcript_perturbation = nullif(@transcript_perturbation,''),
allele_frequency = nullif(@allele_frequency,'');

LOAD DATA LOCAL INFILE './experimental_intervention.csv' INTO TABLE experimental_intervention FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @type, @description)
SET type = nullif(@type,''),
description = nullif(@description,'');

LOAD DATA LOCAL INFILE './efficacy.csv' INTO TABLE efficacy FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @type, @description)
SET type = nullif(@type,''),
description = nullif(@description,'');


LOAD DATA LOCAL INFILE './gene.csv' INTO TABLE gene FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @name, @ensembl_id)
SET name = nullif(@name,''),
ensembl_id = nullif(@ensembl_id,'');

LOAD DATA LOCAL INFILE './drug_experiment.csv' INTO TABLE drug_experiment FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @experiment_id, @dosage, @response, @response_unit, @response_type, @drug_id)
SET drug_id = nullif(@drug_id,''),
response_type = nullif(@response_type,''),
response_unit = nullif(@response_unit,''),
response = nullif(@response,''),
dosage = nullif(@dosage,''),
experiment_id = nullif(@experiment_id,'');

LOAD DATA LOCAL INFILE './drug_target_gene.csv' INTO TABLE drug_target_gene FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @n_quantitative, @cv, @IC50_nM, @AC50_nM, @EC50_nM, @potency_nM, @Ki_nM, @Kd_nM, @n_qualitative, @total_n, @confidence, @pchembl_d, @pchembl_t, @known_selectivity_index, @gene_id, @drug_id)
SET drug_id = nullif(@drug_id,''),
gene_id = nullif(@gene_id,''),
known_selectivity_index = nullif(@known_selectivity_index,''),
pchembl_t = nullif(@pchembl_t,''),
pchembl_d = nullif(@pchembl_d,''),
confidence = nullif(@confidence,''),
total_n = nullif(@total_n,''),
n_qualitative = nullif(@n_qualitative,''),
Kd_nM = nullif(@Kd_nM,''),
Ki_nM = nullif(@Ki_nM,''),
potency_nM = nullif(@potency_nM,''),
EC50_nM = nullif(@EC50_nM,''),
AC50_nM = nullif(@AC50_nM,''),
IC50_nM = nullif(@IC50_nM,''),
cv = nullif(@cv,''),
n_quantitative = nullif(@n_quantitative,'');

LOAD DATA LOCAL INFILE './drug.csv' INTO TABLE drug FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @DT_explorer_internal_id, @name, @inchikey, @inchi, @std_smiles)
SET std_smiles = nullif(@std_smiles,''),
inchi = nullif(@inchi,''),
inchikey = nullif(@inchikey,''),
name = nullif(@name,''),
DT_explorer_internal_id = nullif(@DT_explorer_internal_id,'');

LOAD DATA LOCAL INFILE './drug_synonyms.csv' INTO TABLE drug_synonyms FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS (id, @DT_explorer_internal_id, @name, @drug_id)
SET drug_id = nullif(@drug_id,''),
name = nullif(@name,''),
DT_explorer_internal_id = nullif(@DT_explorer_internal_id,'');
