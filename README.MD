
For a command line usage : python bio_to_upset_plots.py DIR

-  where DIR is the directory where the aligned result files (BIO format) are stored

For an external usage
- call the files_2_cat function with a path containg BIO files
- by default the names displayed on the graphs are derived from the filenames
- Options:
 - img_format : specifying another image format (defaul : png)
 - dic_names  : alternate names to display for the graph
  - Works IF
    - dic_names associates a name to each path
    - OR dic_names associates a name to each filename
  - else:
    - use default names (e.g. filenames)


TODO:
- better handling of BIO format
- offer the possibility to select the appropriate columns
- display messages  regarding the quality of the data
