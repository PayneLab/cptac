```
"WELCOME TO CPTAC!

cptac is a package of adventure, danger, and low cunning. In it you will
explore some of the most amazing data ever seen by mortals. No computer should
be without one!"
```
(a paraphrase of the opening to the beloved 1977 text-based adventure game [Zork](https://en.wikipedia.org/wiki/Zork) created at MIT by Tim Anderson, Marc Blank, Bruce Daniels, and Dave Lebling)

This document will explain what cptac is, how it is structured, and why we chose to structure it that way.<br><br><br>


# What the Package Does

The cptac package distributes cancer datasets generated by the National Institute of Health's [Clinical Proteomic Tumor Analysis Consortium](https://proteomics.cancer.gov/programs/cptac) (CPTAC). Each dataset consists of several tables of proteogenomic data about a particular type of cancer, gathered from approximately 100 tumor tissue samples and 30 normal tissue samples. Types of data in the tables include protein expression levels, RNA expression levels, protein phosphorylation levels, somatic mutations, and other molecular data. Each dataset also contains one or more tables of clinical data (age, gender, medical history, clinical attributes, etc.) for the patients the samples came from.

Our package delivers these data tables in a Python environment, as pandas DataFrames. These function similarly to R dataframes, and are ready to be directly piped into whatever statistical, machine learning, graphing, or other Python packages a researcher may want to use to analyze the data. Although the tables come to us in a variety of formats, we carefully parse them to be in the same format across all datasets, so you don't have to learn an entirely new format when you switch from working with one dataset to another. This means that a researcher can easily take an analysis they've run on one dataset, and run it almost seamlessly on any other dataset.

cptac is freely distributed through the Python Package Index (PyPI). Anyone can install the package using the pip program:

`pip install cptac`
<br><br><br>

# Basics of How the Package Works

Currently (August 2022), the package contains 10 datasets:

*   Breast cancer (BRCA)
*   Clear cell renal cell carcinoma (CCRCC)
*   Colon cancer (COAD)
*   Endometrial cancer (UCEC)
*   Glioblastoma (GBM)
*   Head and neck squamous cell carcinoma (HNSCC)
*   Lung adenocarcinoma (LUAD)
*   Lung squamous cell carcinoma (LSCC)
*   Ovarian cancer (OV)
*   Pancreatic ductal adenocarcinoma (PDAC)

This list will continue to grow as the consortium generates more datasets. All of these data files would be too much for a user to download all at once when they install the package; additionally, PyPI limits the size of our package to 60 MB. So, instead of storing the data files as part of the package, we store them remotely, and provide a function in cptac, cptac.download, for downloading the dataset files that the user wants to work with. The files are stored on [Zenodo.org] (https://zenodo.org/record/7897498).

Once a user has downloaded the data source files for their desired dataset, they can use cptac to access and work with that data in a Python interpreter or script. We chose to have users interface with each dataset through a Python class representing that dataset. Each dataset class reads and parses the dataset's files into dataframes when it's initialized, stores the dataframes in a private dictionary, and provides "get" methods to allow the user to access the dataframes. It also provides functions that handle complex joining between different dataframes.

For example, this is how a user would download the endometrial dataset, load it in their python interpreter, and then get the proteomics dataframe:

![basic_loading_demo](imgs/basic_loading_demo.png)

<br><br><br>

# How We Did it

In the following sections, we'll go through how we implemented the package.

<br><br><br>


# Data Storage and Versioning

The data for each dataset are currently stored on [Zenodo.org] (https://zenodo.org/record/7897498). Zenodo provides a REST API to interact with the datafiles, which we access via the python `requests` package. As Zenodo does not support a filesystem (that is--without zipping all the data into one file), each file has a prefix specifying the source, cancer, and datatype it is associated with. For instance, the file that stores the Breast cancer protemic data from the University of Michigan is called `Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`. Accordingly, the file is stored on Zenodo as `umich-brca-proteomics-Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`. When the data is downloaded via the `cptac` package, that prefix is removed and used to build out the path to the file. For instance, the above file would be saved with the path `~/data/umich-brca/Report_abundance_groupby=protein_protNorm=MD_gu=2.tsv.gz`, where `~` is the base directory for cptac.

As you can see, we refer to each dataset by an acronym or short name. BRCA is for breast cancer, CCRCC is for clear cell renal cell carcinoma, UCEC is for uterine corpus endometrial carcinoma (endometrial cancer), and so on. Some files contain data for multiple cancer types at once. These files have been given the short name ALL_CANCERS. For instance, the file containing the clinical metadata, `clinical_Pan-cancer.May2022.tsv.gz` is stored on Zenodo as `mssm-all_cancers-clinical-clinical_Pan-cancer.May2022.tsv.gz`.

Previous cptac users may recall that the cptac package employed extensive data versioning logic. This is because in earlier releases of this package were developed while data was still being processed and updated. However, at the time of this release, each data file is fully finished and will no longer be modified. Though new datasets may (and will hopefully) continue to be added, any file that is currently accessible will not be changed. As a result, there is no reason to not keep users updated with the most recent version. The package updates now use Zenodo's simpler versioning system, and always stay up-to-date with the most recent data available.

<br><br><br>

# Downloading the Data

When a user installs the cptac package, it is installed to a package installation directory at one of the following locations (if they're using Anaconda):

`~/anaconda3/lib/python3.7/site-packages/cptac` if they've installed it to their base environment, or 

`~/anaconda3/envs/[ENVIRONMENT NAME]/lib/python3.7/site-packages/cptac` if they've installed it to a different environment.

 Within the package installation directory after the package is first installed, the file structure for the home directory looks like this:

```
cptac/
    cancers/
    tools/
    utils/
    __init__.py
    exception.py
    version.py
```

After a user first imports cptac into python (with an internet connection), a new `data` folder is added, changing the home directory to:
```
cptac/
    cancers/
    data/
    tools/
    utils/
    __init__.py
    exception.py
    version.py
```

From this point onwards, the home `cptac` directory will be referred to as `~`. A description of these files is as follows:
*   `cancers/`: Contains the bulk of the essential code for this package, such as for loading, saving, and accessing datasets. The organization is explained in more detail in the next section.
*   `data/`: Contains the data files as well as an index file used to look them up. We use `~/data/index.tsv` to associate a specific combination of cancer_name, source, and datatype with the file name on zenodo, and also to store a checksum to make sure that the local and remote files are identical. The subfolders in `~/data` are structured in the format `{source}-{cancer_type}`. For instance, the Breast cancer (brca) data from the University of Michigan (umich) is stored in the folder `~/data/umich-brca`.
*   `tools/`: Contains modules essential for the under-the hood working of cptac, but that is generally not used by the user. Examples include implementing downloading a file via the requests package, modifying the standard output for error messages, or dataframe operations that are frequently used when reading in data files.
*   `utils/`: Contains modules and functions designed for the user to simplify common operations, such as performing bulk t-tests across multiple genes, simplifying working with a multi-index, or associating gene/protein data with specific metabolic pathways.
*   `__init__.py`: Initializes the package. Sets a custom exception hooks, creates `~/data/index.tsv`, and makes the cancer objects accessible to the user.
*   `exception.py`: Defines exceptions that the exception hook will catch.
*   `version.py`: Specifies the version of the package. Used by PyPi.

<br><br><br>

# Structure of `~/cancers/`
As mentioned above, the `~/cancers/` folder contains the bulk of the essential code for this package. There are two type of custom classes in this folder: **cancers** and **sources**.

**Cancers** all inherit from `~/cancers/cancer.py`, and are the primary object that users interact with. Cancers provide the functionality displaying, joining, and requesting dataframes in an easy-to-use interface. Cancers store a dictionary of sources (*self._sources*), which they use to get the data that they then use *(described below)*. All `.py` files in `~/cancers/` specify cancer objects with the exception of `~/cancers/source.py`. In addition, almost all of the functionality of these cancer objects is stored in the parent class `~/cancers/cancer.py`, with the other cancer objects only filling their *self._sources* attribute with their cancer-specific sources. As such, **one of the best places to start for a developer is with `~/cancers/cancer.py`**.
<br><br>
**Sources** all inherit from `~/cancers/source.py`, and are responsible for locating the data files, reading them in, and then modifying the dataframes to have a consistent formatting. Other than the parent class, all source classes are tied to a specific data source and cancer type, and are stored within folders associated with that source. For instance, the code for reading in data files from the University of Michigan *(umich)* is stored in `~/cancers/umich/`, and the code for reading in Washington University's *(washu)* Breast cancer *(brca)* files is stored in `~/cancers/washu/washubrca.py`. These files are primarily comprised of *load_functions*, which specify how each dataframe should be loaded, and need to be custom-written every time a new dataset is added  *(see [02_add_new_dataset.md](02_add_new_dataset.md#adding-the-code-for-a-new-dataset))*. Common operations like locating a file, saving a dataframe, or requesting a download are stored in the parent class *(`~/cancers/source.py`)*, and generally do not need to be rewritten,
<!--In the next update, we will likely have a class in the middle that is specific to a source but not a cancer, such as `~/cancers/washu/washu.py`. This is because many of the load functions across all washu's cancer types are identical (i.e. The proteomic data is loaded the same regardless of source). Putting this repeated code in a class that all can inherit will reduce code bloat.-->


<br><br><br>

# Loading and Accessing Datasets
Once a user has downloaded a dataset, they can access it in their Python interpreter, or in a script. Here, we will provide a step-by-step explanation of what happens as the user executes basic commands for accessing and working with a dataset, using the following basic code as an example:
```
import cptac
en = cptac.Ucec()
prot = en.get_proteomics('umich')
```

*   `import cptac`
    *   This imports the package into the current namespace and runs everything in `~/__init__.py`, which does the following:
        - Defines some helper functions (cptac.list_datasets, cptac.how_to_cite, etc.)
        - Sets up error and warning handling for cptac-generated errors and warnings (more on that later)
        - Checks that the package is up-to-date and re-creates the index file (`~/data/index.tsv`, more on that later)
*   `en = cptac.Ucec()`
    *   This creates an instance of the Endometrial class, which is defined in `~/cancers/ucec.py`, and calls the class's `__init__` function. 
    *   The Endometrial `__init__` function first calls the `__init__` function of its parent class in `~/cancers/cancer.py`. This initializes several private variables, such as the valid dataframe names. <!--This private variable may soon become obselete, as it is handled by the index.tsv initialization. -->
    *   Then, the Endometrial `__init__` function defines two essential private variables: `self._data`, which is the dictionary that holds the dataset's dataframes once they have been loaded; and   `self._load_functions`, which define how data should be fetched if it is not already loaded into `self._data`. More details are given in the file [02_add_new_dataset.md](02_add_new_dataset.md).
*   `prot = en.get_proteomics('umich')`
    *   Like all the other `get` functions, `get_proteomics` is little more than an alias for the `get_dataframe` function in `~/cancers/cancer.py`. This call is identical to `en.get_dataframe('proteomics', 'umich')`.
    *   `get_dataframe` checks to make sure that the <u>datatype</u> and <u>source</u> combination is valid, then calls the requested source's `get_df()` function, as so: `self._sources['umich'].get_df('proteomics')`. 
    *   `get_df` checks whether the data has already been loaded by checking the `self._data` variable. If it has not been, it calls the associated load function to load that data, downloading if needed. Finally, that dataframe is returned, and passes up the chain until it is eventually stored in the user's `prot` variable. 


# Using the package without an internet connection

Wherever possible, cptac is written so it can still be used without an internet connection. Obviously, you can't download data files without internet, but if a user has already downloaded the files for a dataset, they load and work with the dataset just fine.


# Multi-level column indices

We take great care to ensure that none of the dataframes have duplicate columns. Most of the omics dataframes have gene names as the column headers, but sometimes there are duplicated gene names, obviously in phosphorylation and acetylation dataframes, but also sometimes in other omics dataframes. To solve this problem, we use a multi-level index as the column headers, where needed. We implement this with a pandas MultiIndex (see [this page for documentation from pandas](https://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html)). 

Briefly, a multi-level index is where each data series has multiple keys associated with it, kind of like having a named tuple for the index. This allows us to include the keys necessary to uniquely identify each column. Every data series will have the same number of keys, and each key belongs to a "level" of the multiindex--the first keys for all the data series make up the first level of the multiindex, the second keys make up the second level, and so on.

The phosphorylation and acetylation dataframes have four levels in their column index: Name (gene name), Site, Peptide, and Database_ID. Other omics dataframes with multi-level column indices have some or possibly all of these levels, as needed to uniquely identify each column. 

First, we will show what this sort of multiindex would look like for a simple mockup dataframe, for illustration:


![multiindex_mockup](imgs/multiindex_mockup.png)


Each column has its own key value for each level (Gene, Site, Peptide, Database_ID) of the multiindex. In other words, each column has its own gene, site, peptide, and database ID associated with it. However, note that when several adjacent columns have the same key for a particular level of the column index, that key is only printed for the first column, to improve readability. For example, in the table above, the first 4 columns all are for gene1, but that label is only shown for the first column. The first 3 columns also all have the same site, so the site is only shown for the first of the 3. The 4th column has a different site, so that's shown.

To help illustrate the concept of a multiindex, though, we show the image below, which is the same mockup dataframe as above, but showing every key for every column--in other words, not hiding repeat keys. Thus, the first 4 columns all have "gene1_phos" as their key for the "Name" level, and the first 3 columns all also have "S12" as their key for the "Site" level. In normal printing, the repeated key would only be shown for the first column in the group.

![multiindex_verbose_mockup](imgs/multiindex_verbose_mockup.png)


Below, we print an example of what an actual dataframe with a multi-level column index looks like. This is the phosphoproteomics dataframe from the clear cell renal cell carcinoma (CCRCC) dataset.


![ccrcc_multiindex](imgs/ccrcc_multiindex.png)


Multi-level indexing has many advantages as a method for uniquely identifying columns. It makes it simple to select columns based on any combination of keys, because each key is a different element in the tuple, so it's easy to look at just a single key, or multiple keys, when selecting data. For example, we could select all columns with a particular gene name, or look at columns with a particular type of site in the "Site" level of the index. 

We provide a helper function, DataSet.reduce_multiindex, to make it easier for users to work with dataframes that have multi-level column indices. It has a `levels_to_drop` parameter that allows a user to drop one or multiple levels from a column multiindex, and it warns the user if such an operation creates any duplicated column headers. It also has a `flatten` parameter that provides the option to concatenate the values from each level of the index into a single string for each column header, thus create a single-level index, but avoiding duplicate column headers. By default, it uses an underscore as the separator character between the values joined from different levels, but the user can specify a different character.


# Merging dataframes

cptac provides five built-in functions for joining different dataframe types:



*   join_omics_to_omics
*   join_omics_to_mutations
*   join_metadata_to_metadata
*   join_metadata_to_omics
*   join_metadata_to_mutations

These functions handle numerous details that would normally make joining tables a complicated process for users. As a result, it is very easy for users to join columns from multiple tables, for their analyses. Here are some features of the join functions:



*   When a column is selected from a dataframe, the name of the dataframe it came from is appended to the value of the Name level of the column header, to avoid confusion about which columns are which in the joined dataframe.
*   When two dataframes are joined that have a different number of levels in their column indices, levels filled with NaNs are added to the column index of the dataframe with less levels, so they are compatible for joining.
*   The join functions use full outer joins. If this causes any cells to be filled with NaN because there wasn't a value in that column for a particular sample that did exist in the other dataframe, the function warns the user.
*   When joining an omics or metadata dataframe to the mutations dataframe, the join function automatically handles the work of selecting the specified mutations from the somatic mutation dataframe, and reformatting the table to be joinable with an omics or metadata dataframe.
*   When there are multiple mutations for a single sample in a particular gene, the join_omics_to_mutations and join_metadata_to_mutations functions have a `mutations_filter` parameter that allows you to filter those multiple mutations down to a single mutation, using either a default priority, or using a priority the user supplies.

(Note that some metadata tables, such as treatment and medical_history, have multiple rows for each sample, and as a result cannot be used in join functions.)


# Package checks if it's up-to-date

As noted earlier, cptac automatically checks that it's up-to-date every time a user imports it. This is done in the `~/__init__.py` file, using the function `init_files()`. This function requests the record from Zenodo that is specified by the hard-coded concept-doi, which automatically points to the most upd-to-date version of the data. Zenodo's response is used to generate a new index file on-the-fly, using the most-up-to-date data.


# Exceptions and warnings

There are various reasons that an exception or warning might need to be raised by cptac. For example, a user may pass the wrong type of parameter to a function, or request a dataframe or column that doesn't exist in that particular dataset, or create duplicate levels in a column multiindex when they request to drop a particular level. If a problem is severe enough that the entire operation needs to be cancelled, cptac does so by raising an exception. If the issue is not fatal, and is merely something the user should be aware of, cptac issues a warning, using the standard Python warnings module. This design is clean and simple. Additionally, since all warning and exception messages are sent to stderr, not stdout, it is easy for a user to separate the output of a command or script from the warning or error messages--they simply just pipe stderr and stdout to different output files.

The disadvantage to exceptions and warnings is that their syntax can be confusing or intimidating to users with less coding experience, and that includes many biologists who we hope will use this package. To solve this problem, we created custom exception and warning classes for all exceptions and warnings generated by cptac. They are defined in `cptac/exceptions.py`. Then, in the package initialization code in `cptac/__init__.py`, which is run immediately when a user imports the package, we set a custom exception hook (sys.excepthook) and warning displayer (warnings.showwarning). When an exception or warning is routed through one of these, they check whether it inherits from CptacException or CptacWarning, the base classes for all the exceptions and warnings defined in `cptac/exceptions.py`. If it is, the exception or warning message is parsed and printed in a simple, non-intimidating format, and the user can go on their merry way. If the exception or warning originated from another source, however, it is printed in the default format. We've included an example below.


![errors_demo](imgs/errors_demo.png)


# Utils sub-Module

cptac includes a sub-module, cptac.utils, which provide various utility functions for various tasks such as performing t-tests, getting interacting proteins, and finding which genes are most frequently mutated in a dataset.