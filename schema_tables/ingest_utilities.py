import pandas as pd
import numpy as np
import difflib
import os
import json
import networkx as nx

import synapseclient


# TODO: add "strong" typing

# -*- coding: utf-8 -*-
"""###########################################################################
This module provides a collection of utilities to ingest data following  
a prescribed schema (e.g. by ingestion SOPs) into a relational database and
related Synapse tables following the schema.
###########################################################################"""



# handle synapse login
# TODO: think what's a good way to handle that as part of a workflow
def synlogin():
    syn = synapseclient.Synapse()
    syn.login()

    return syn



# after outer join, keep columns from right table if not NaN; otherwise keep columns from left table
def squash_row(x, columns):

    for col in columns:
        if not pd.isnull(x[col + "_y"]):
            x[col] = x[col + "_y"]

    return x



# update a normalized table
def update_table(df, df_update, primary_key):

    if not df.columns == df_update.columns:
        raise ValueError("Different column sets. Both old and updated table should have the same set of columns.")

    if primary_key != "index" and (not primary_key in df.columns):
        raise ValueError("Primary key must be in update table column set.")

    df_updated = df.reset_index().merge(df_update, how = "outer", on = primary_key)
    df_updated = df_updated.apply(squash_row, args = (df.columns,), axis = 1)
    
    # ensure to get only schema columns
    df_updated = df_updated[df.columns]

    return df_updated



# function to normalize a table (e.g. dedup, drop nan's)
def normalize_table(df, primary_key):

    try:
        # if valid primary key has been provided normalize df
        df = df.reset_index()
        df_norm = df.dropna().drop_duplicates(subset = [primary_key])
        return df_norm
    except KeyError:
        # if the primary key is not in the df; then return the same df w/o changes
        print("Specified primary key is not in table schema. Proceeding without table changes.") 
        return df



# download synapse table as a pd dataframe; return table schema and etags as results too
def get_synapse_table(synapse_id, syn):

    results = syn.tableQuery("SELECT * FROM %s" % synapse_id)
    df = results.asDataFrame()

    return df, results



# function to ingest normalized table data into existing table on Synapse with a particular primary key following
# respective table schema/data ingest sub-SOP, synapse_id and synapse connection (i.e. syn)
def ingest_entity_table_data(df_update, synapse_id, primary_key, syn):
    
    # ensure update table is normalized
    df_update = normalize_table(df_update, primary_key)

    # get existing DB data and table schema
    df, results = get_synapse_table(synapse_id, syn)

    # update table with update dataframe 
    df = update_table(df_update, primary_key)

    # store table back to Synapse
    table = Table(synapse_id, df, etag=results.etag)
    syn.store(table)
   
    # return updated table
    return df



# get a table update order given a SOP-derived update table and a db schema, 
 
# if table x contains a foreign key (i.e. relation to another table y)
# add the related table y in the list of tables to update before table x.
# This determines update order, so that normalized entity tables are 
# updated before junction updates, for instance. Essentially, 
# perform a topological sort on the db_schema subgraph of table
# relationships (assuming it is a DAG, which should be the case for a 
# well-temperered DB); this returns the update order of tables
# 
# the DB relations subgraph is determined by subset of columns (and respectively their container tables)
# that need to be updated based on the given SOP-derived update table
def get_update_order(df_update, db_schema):


    # construct db_schema relationship graph;
    # this graph is a DAG hence nx ignores adding node duplicates, so no need to check if a node exists
    DB_relations_graph = nx.DiGraph()
    for table_name, table_schema in db_schema["tables"].items():
        for attribute, attribute_map in table_schema.items(): 

            # only add tables as DB relations graph nodes, if they have attributes for update in the update data frame
            if not attribute in db_schema["ignore_attributes"] and attribute in df_update.columns:
                # check if the attribute originates from a different table (i.e. it is a foreign key)
                if isinstance(attribute_map, dict):
                   # desideratum: update first the source table from which this attribute originates (i.e. it is attribute_of)
                   # then update the current table (i.e. table_name)
                   # therefore add the following directed link indicating this relationships in the topo sort
                    
                   DB_relations_graph.add_edge(attribute_map["attribute_of"], table_name)
                else: # an attribute of this entity/table has data in the provided update df; so add this table as a node
                   DB_relations_graph.add_node(table_name)

                
    # perform topological sort on the DB relations subgraph to get update order
    update_order = list(nx.topological_sort(DB_relations_graph))
    print(update_order)
    return update_order, DB_relations_graph



# update tables
def update_tables(df_update, db_schema_path):

    # login to Synapse
    syn = synlogin()
    
    with open(db_schema_path, "r") as db_f:
        db_schema = json.load(db_f)

    # list of tables to update; tables ranked by update order
    tables_to_update, DB_relations_graph = get_update_order(df_update, db_schema)


    for table_name in tables_to_update:

        table_schema = db_schema["tables"][table_name]
        sub_sop_columns = list(set(table_schema.keys()).difference(set(db_schema["ignore_attributes"])))
         
        # get sub-SOP columns from update data frame if they exist
        df.loc[:, df.columns.isin(list('BCD'))]
        df_sub_sop = df_update[sub_sop_columns]

        # SOP and table schema columns may have different names; keep track of the mapping
        table_schema_columns = {}

        # keep track of related entitie's attributes (i.e. foreign keys; useful for joins)
        foreign_keys = {}

        for sop_column in sub_sop_columns:

            attribute = table_schema[sop_column]

            if isinstance(attribute, dict):
                # this attribute links to a different entity table (it is a foreign key)
                #TODO: might want to be more explicit about that - instead of assuming a dictionary indicates links to a different table can have an explicit term/key
                # get its schema column name based on the connecting foreign key to the other table (i.e. connected_by attribute)  
                schema_column = attribute["connected_by"]   
                
                # get table/entitiy where this attribute/foreign key originates
                # store foreign key name and its corresponding SOP column name
                foreign_keys[attribute["attribute_of"]] = attribute["connected_by"]
            else:
                schema_column = attribute
                
            table_schema_columns[sop_column] = schema_column 

        # update column names of sub-SOP data frame to match table schema
        df_sub_sop.rename(columns = table_schema_columns, inplace = True)

        # not all columns in the sub-SOP column-set may be provided as part of the update data frame;
        # add the missing columns to data frame filled w/ NaNs
        for missing_column in set(table_schema_columns).difference(df_sub_sop.columns):
            df_sub_sop[missing_column] = np.NaN
       
        # perform joins if needed; 
        if foreign_keys:

            for foreign_table_name, foreign_key in foreign_keys.items():
                
                foreign_table_schema = db_schema["tables"][foreign_table_name]

                # get foreign table from synapse
                foreign_table, foreign_table_schema_syn = get_synapse_table(foreign_table_schema["synapse_id"], syn)
                
                # join to get foreign key values in the update table data frame;
                # note that the update order of tables ensures the foreign table has already been updated as needed
                # hence the left join below updates foreign key values in this table as well.
                df_sub_sop = df_sub_sop.merge(foreign_table[foreign_table_schema["primary_key"]].reset_index(), how = "left", left_on = foreign_key, right_on = foreign_table_schema["primary_key"])

                # rename and drop columns to match DB schema
                df_sub_sop.drop([foreign_key], axis = 1, inplace = True)
                df_sub_sop.rename(columns = {"index":foreign_key}, inplace = True)
        
        # ensure the sub-SOP table contains exactly the set of columns required by the DB schema
        df_sub_sop = df_sub_sop[table_schema["columns"]]
            
        # ingest update data-frame in corresponding Synapse table
        # TODO: check if ingest_entity_table_data accounts for many-to-many relationships and joins on index correctly
        # ingest_entity_table_data(df_sub_sop, table_schema["synapse_id"], table_schema["primary_key"], syn)
