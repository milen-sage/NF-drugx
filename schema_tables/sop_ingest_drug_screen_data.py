import pandas as pd
import os

import glob
from ingest_utilities import update_tables

# TODO: add "strong" typing


# -*- coding: utf-8 -*-
"""###########################################################################
This module provides a collection of utilities to ingest model and drug 
experiment data according to the SOP specification here: TODO

Ingesting model and drug experiment data following the SOP spec requires 
as well the specification of the following (human-readable) sub-SOPs ingesting data into the 
normalized and juncture entities below:
    * disease: SOP - TODO
    * symptom: SOP - TODO
    * organism: SOP - TODO
    * cell_line: SOP - TODO
    * model: SOP - TODO
    * experiment: SOP - TODO
    * study: SOP - TODO
    * funder: SOP - TODO
    * experimental_intervention: SOP - TODO
    * efficacy: SOP - TODO

For more information on the entities above, please refer to the specific 
SOPs and the drugx DB schema here: TODO 
###########################################################################"""


# drug-gene DB json-ld config
drugx_schema_path = "./drugx_db.json"

# drug-gene DB update tables following SOP
update_tables_path = "./sop_data/"

# drug-gene DB drug screen data SOP tables
table1 = os.path.join(update_tables_path, "combined_NF_screening_data_released_only_summary_data_table_1_Jan30.csv")
table2 = os.path.join(update_tables_path, "combined_NF_screening_data_released_only_summary_data_table_2_Jan30.csv")

table1 = pd.read_csv(table1)
table2 = pd.read_csv(table2)

update_tables(table2, drugx_schema_path)
