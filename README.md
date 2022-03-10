# PCR_Experiment_Setup
A Shiny app to integrate dynamic experimental setups, sample metadata, and freezer inventory

## General workflow

### Purpose
The purpose of this app is to integrate several separate but related sets of information to make experimental design, finding samples, and updating sample information in a laboratory information management system (LIMS).

This app takes two inputs:
1. A sample list that includes unique experimental conditions (dilution and volume) for each sample
2. An sample inventory including sample volume and freezer location

These two inputs are integrated to give the user information about the samples in addition to the sample freezer location. Finally, this information is re-formatted for use in an procedure master index sheet of all assayed samples. This process produces three files which are downloaded by the user:

1. Sample Inventory: Although a small version is presented here, actual sample inventories can contain millions of cells. This file contains sample metadata and location information specifically for the sample to be used. 
2. Inventory LIMS Update Sheet: Often, information on how the experiment was performed must be recorded in a LIMS. This can be a cumbersome process if samples are not easily organized as they are in the experiment. This file utilizes the unique ID number associated with a sample by the LIMS to automatically update some information and provide fast access to samples so that other information may be updated. 
3. Master Sample Sheet Output: This file contains formatted sample and experiment metadata to be added to a master index of all samples used for a given procedure
