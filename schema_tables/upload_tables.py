import os
import pandas as pd
import synapseclient
from synapseclient import Table, Column
from synapseclient.table import build_table
from synapseclient.exceptions import SynapseHTTPError

import json

# pandas column datatypes map to Synapse table column types
DTYPE_2_TABLETYPE = {'?': 'BOOLEAN',
                     'd': 'DOUBLE', 'g': 'DOUBLE', 'e': 'DOUBLE', 'f': 'DOUBLE',
                     'b': 'INTEGER', 'B': 'INTEGER', 'h': 'INTEGER', 'H': 'INTEGER',
                     'i': 'INTEGER', 'I': 'INTEGER', 'l': 'INTEGER', 'L': 'INTEGER',
                     'm': 'INTEGER', 'q': 'INTEGER', 'Q': 'INTEGER',
                     'S': 'STRING', 'U': 'STRING', 'O': 'STRING',
                     'a': 'STRING', 'p': 'INTEGER', 'M': 'DATE'}

# login Synapse
syn = synapseclient.login()

# project location
project_root_path = "/Users/milen-sage/workspace/drug-gene-db/schema_tables/"

# db metadata location
drugx_meta_file = os.path.join(project_root_path, "drugx_db.json")

with open(drugx_meta_file, "r") as m_f:
    drugx_meta = json.load(m_f)

# project workspace
db_project_syn_id = drugx_meta["synapse_id"]


# iterate over DB tables get corresponding schema table csv and update the corresponding Synapse table
# currently schema table csv is obtained assuming file name is table_name + ".csv"; in the future we will have 
# direct uris

for table_name in drugx_meta["tables"]:

    print("=====================================")
    print("Processing table " + table_name)
    
    # get location of update table
    update_table_file_name = os.path.join(project_root_path, table_name + ".csv")
    # read csv table as a dataframe
    update_table = pd.read_csv(update_table_file_name)
    
    # get location of update target table on Synapse
    table_syn_id = drugx_meta["tables"][table_name]["synapse_id"]
    
    # remove all columns and recreate schema of existing table
    schema = syn.get(table_syn_id)
    cols = syn.getTableColumns(schema)
    for col in cols:
        schema.removeColumn(col)

    for col in update_table:
        column_type = DTYPE_2_TABLETYPE[update_table[col].dtype.char]
        new_col = syn.store(Column(name = col, columnType=column_type))
        schema.addColumn(new_col)

    # update table name if different
    if table_name != schema.name:
        schema.name = table_name
    
        # update table schema on synapse
        schema = syn.store(schema)
    
    # the schema of the target table now matches the schema of the update table
    # store update table on Synapse to replace target table
    table = syn.store(Table(table_syn_id, update_table))
    
    print("Done.")
    print("=====================================")
    

with open(drugx_meta_file, "w") as m_f:
    json.dump(drugx_meta, m_f, indent = 3)
