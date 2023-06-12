# How to Add a New Dataset

## Adding the Data Files for a New Dataset

1. Download the data files onto your machine.
2. gzip the files (unless they're Excel files--then just leave them uncompressed).
3. Rename the files to include the prefix in the following format: `{source}-{cancer type}-{datatype}-{original filename}`.
4. Create a new version on the Zenodo repository, upload the data, and publish.
<br><br><br>

## Adding the Code for a New dataset

For a new dataset to function, each dataframe must contain its own load function. This is stored in the class given by `~/cancers/{source}/{source}{cancer}.py`, where **{source}** is the lowercase abbreviation for your source and **{cancer}** is the lowercase abbreviation for your cancer. This structure is given in [00_why_we_did_what_we_done](00_why_we_did_what_we_done#structure-of-cancers). To add this code, do the following:<br>
1.  If the above file does not exist, copy [child_source_template.py](child_source_template.py) into the appropriate folder and give it the name `{source}{cancer}`. Then, open the file and rename the class `NameOrAcronym` to `{Source}{Cancer}`, where `{Source}` is the acronym for the source with the first letter capitalized, and `{Cancer}` is the acronym for the cancer with the first letter capitalized.
2.  Add the entry `{datatype}: {name-of-data-file}` to the private `self.data_files`.
3.  Add the entry `{datatype}: load_{datatype}` to the private `self.load_functions`.
4.  Inside the class, define your loading function. Use the template as a guide.
5.  If you had to follow step 1, open the file containing the cancer object, which is `~/cancers/{cancer_type}.py`. You will need to add two lines of code:
    *   Underneath the other imports, add the line `from cptac.cancers.{source}.{source}{cancer} import {Source}{Cancer}`
    *   In the `__init__` function, add the line `self._sources["{source}"] = {Source}{Cancer}(no_internet=no_internet)`

**Important:** Whenever you're testing changes to your code, make sure to locally install the package using `pip`, using the following instructions. These instructions will take the local copy of the package that you've been editing and install it in your Anaconda environment's package installation directory. This will make it so that when you've opened a Python prompt or a Jupyter Notebook from that Anaconda environment and then import the package, you'll be importing your edited version of the package. This allows you to test the edits you've made, without having to push them to PyPI. So, to install your locally edited version of the package:
1.  Open your Anaconda prompt or terminal
2.  Activate your development environment (`conda activate MyEnvironment`, subbing in the name of your environment)
3.  Navigate to the cptac directory that contains the `setup.py` file (which is the upper cptac directory, not the lower one). `pip` reads this file to know how to install the package.
4.  In that directory, run this command: `pip install .` (don't forget the dot--it's a reference to your current directory, telling pip to build the package based on the `setup.py` file it finds in the current directory)
5.  Alternatively, if you're in a different directory, you could run `pip install /path/to/cptac/directory/with/setup/py/file`, subbing in the proper path to the cptac `setup.py` file, and replacing / with \ if you're on Windows. `pip` will follow that path, find the `setup.py` file, and then install the package based off of it.

<br><br><br>
## General Dataframe Formatting Requirements

First, I'll give you two pictures showing the general format that you want to parse all the data tables to conform to:

Here is an example of a table with a single level column index:


![single_indexed_df](imgs/single_indexed_df.png)


Here is an example of a table with a multi level column index:


![multiindexed_df](imgs/multiindexed_df.png)


These tables conform to these requirements:

*   Each column represents a particular variable. For omics tables, each column represents a gene/protein/etc. For metadata tables, each column represents a clinical variable (age, weight, etc.).
*   Each row represents a specific sample.
    *   Samples are indexed by their patient ID.
    *   All normal samples have '.N' appended to the end of the patient ID.
    *   Some patients donated both a tumor sample and a normal sample. Both samples should have the same patient ID; the only difference between the two will be that the normal sample has '.N' appended to the end of it.
        *   In some datasets, such as HNSCC, there are also cored normal samples. These are marked with a '.C' at the end of the patient ID, and in the clinical dataframe have the value "Normal" in the Sample_Tumor_Normal column, and the value true in a column you create called Cored_Sample. See the "Formatting requirements for specific dataframes" section of this document for details.
*   No duplicate column headers or index values.
    *   Check by calling the 'df.index.duplicated().any()' and 'df.columns.duplicated().any()' functions; both should return False.
    *   Exception: The somatic_mutation dataframe will have multiple rows for each sample, so it will have duplicate index values.
        *   Some metadata dataframes, such as the treatment dataframe in Ovarian or medical_history dataframe in Ccrcc, may also have multiple rows for each sample. If so, they will need to be excluded from join functions. Check the _valid_omics_dfs and _valid_metadata_dfs lists and make sure that tables with duplicate index values aren't included in either of those lists.
    *   There may be duplicate column headers. If that's the case, you'll probably need to use additional identifiers to uniquely identify each column. We use a multi-level column index to accomplish this (see the 00_why_we_did_what_we_done document for more details).
    *   If you have duplicated samples, either there are multiple rows for each patient (as in the treatment or somatic_mutation dataframes), or you need to append a '.N' to the Patient_IDs of the normal samples.
*   Index and column names are standardized:
    *   The value of df.index.name is "Patient_ID" for all tables.
    *   For dataframes with single-level column indices, the value of df.columns.name is "Name".
    *   For dataframes with multi-level column indices, the name of the highest index level is "Name".
    *   This is all ensured by calling the standardize_axes_names function at the end of the dataset loader; see child_dataset_template.py
*   If it has a multi level column index:
    *   Include the minimum number of levels needed to uniquely identify each column. In addition to the gene name, you may need to include a level for a database identifier, a phosphorylation/acetylation site (if it's a phosphoproteomics or acetylproteomics file), and/or a peptide. If it's not a phosphoproteomics or acetylproteomics file, the gene name and database ID levels will probably be sufficient, since database IDs will be unique for different isoforms.
    *   Whatever of the possible multiindex levels you have must have these names, and be in this order: Name, Site, Peptide, Database_ID
*   The key for each dataframe in the self._data dictionary is the standard, all lowercase name for that dataset, as determined by the getter function for it. The getter functions for each table type are named "get_" followed by the exact name of the table type.
*   If the new dataset has a table type not included in any other datasets, you must write a getter for it in the parent DataSet class, found in `cptac/dataset.py`, using the private method DataSet._get_dataframe and passing the tissue_type parameter if it's an omics dataframe.
    *   You'd also need to add the new dataframe's name to self._valid_omics_dfs if it's a valid omics df for the DataSet merge functions, or self._valid_metadata_dfs if it's a valid metadata df for DataSet.append_metadata_to_omics
        *   Note that a dataframe with multiple rows for each sample, like the treatment dataframe in the Ovarian dataset, should **NOT** be a valid dataset for joining
*   Column names are consistent--e.g., all Sample_Tumor_Normal columns should be labeled as such, not as Sample_Status or something else. Rename columns as necessary to match this.

**Formatting requirements for specific dataframes**


* Mass-spectrometry derived tables (proteomics, phosphoproteomics, acetylproteomics)
    * The expression values should be centered around zero. If they are not already centered around zero, then you probably need to find the "Reference Intensity" column and subtract it from all the other columns. Do that automatically when parsing the table. If that isn't the solution or you can't find the column, talk to someone with more information about it.
*   somatic_mutation
    *   The somatic_mutation table should have 3 columns: Gene, Mutation, and Location.
    *   If you're parsing from a MAF file, then the columns you want will probably be named as follows:
        *   Patient_ID (for the index) will be named "Tumor_Sample_Barcode"
        *   Gene will be named "Hugo_Symbol"
        *   Mutation will be named "Variant_Classification"
        *   Location will be named "HGVSp_Short"
    *   Make sure the range of possible values in the Mutation column matches the lists for sorting between truncations and missense mutations in the function DataSet._filter_multiple_mutations. If the possible values don't match the default list of possible values, add an elif statement to handle the special case, as has been done for the Colon, GBM, and HNSCC datasets.
        *   This also affects the DataSet.get_genotype_all_vars function and the utils.get_frequently_mutated functions.
        *   The standard possible values are:
            *   'Frame_Shift_Del'
            *   'Frame_Shift_Ins'
            *   'Nonsense_Mutation'
            *   'Nonstop_Mutation'
            *   'Splice_Site'
            *   'In_Frame_Del'
            *   'In_Frame_Ins'
            *   'Missense_Mutation'
            *   'Silent'
*   clinical 
    *   Must contain a Sample_Tumor_Normal column, which contains either "Tumor" or "Normal" for each sample, according to its status.
        *   Make sure that the "Tumor" and "Normal" values are capitalized properly--can use Series.str.title() to fix if necessary
    *   Don't worry about marking "normal" samples in the clinical table. We treat all samples in the clinical table as tumor samples; I'll explain why. First note that each normal sample comes from a person who also donated a tumor sample. The data in the clinical data pertain to the actual person, not to the individual samples that came from them. Therefore we classify all the rows in the clinical table as "tumor" samples, because if we had data for each normal sample as well, we would just be duplicating the data for whatever tumor sample came from the same person. We do create rows in the clinical table for the normal samples, but every column except for the Sample_Tumor_Normal column is NaN (null) for those rows. Note that adding these NaN rows to the clinical dataframe is already built into the dataset template file--it occurs when we reindex the clinical dataframe with the master index (line 127 in the cptac/devdocs/child_dataset_template.py file). The master index is a list of every sample ID that occurs in any table in the dataset, which we have to create by explicitly going through all the tables, because not all samples are in all tables due to data quality issues--if a sample had some problem for a particular data type, it was dropped from the table for that data type. Creating the master index is done by the "unionize_indices" function on line 123 of the template file. This gets the normal samples from the tables where they appear; then when we reindex the clinical table with the master index, it inserts rows of NaNs for any samples IDs that weren't originally in the clinical table, including the normal sample IDs.
    *   In some datasets, such as HNSCC, there are also cored normal samples. You should mark the patient IDs for these samples with a '.C' at the end in all tables, instead of a '.N'. You also need to give them the value "Normal" in the Sample_Tumor_Normal column, and the value True in a column called Cored_Sample that you need to create in the clinical dataframe. All tumor and non-cored normal samples will have the value False in the Cored_Sample column.

### Tips for Writing the Parser/Loader:

*   If any dataframes are split between two files--such as one file for the tumor sample proteomics, and one file for the normal sample proteomics--they'll have been read into separate dataframes, and you need to merge those into one dataframe.
    *   Make sure that samples coming from the normal file have '.N' appended to their Patient_ID numbers, to keep a record of which ones are normal samples. You may need to do that manually.
*   Multiple metadata dataframes may be contained in one file, e.g. clinical and derived_molecular data might both be in clinical.txt, as in Endometrial. If that's the case, load the whole file and then assign different columns to the proper tables. Ask Dr. Payne if you have questions about which columns should go in what table.
*   If there's a column in the clinical or other metadata dataframe indicating that some samples were excluded from analysis, like in Endometrial, we'll probably just drop those samples from the data. Check with Dr. Payne.
