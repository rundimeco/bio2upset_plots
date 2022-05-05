### Creating upset plots from BIO (CONLL) files

Upset plots are a clever way to present intersection and differences between sets. It is pretty useful for comparing different systems that produce results in BIO format (NER for instance). Contrary to Venn Diagrams, upset plots allows to represent clearly more than three sets intersections and differences.

The upset plots themselves are produced using the [upsetplot python Library](https://upsetplot.readthedocs.io/en/stable/)

### Acknowledgements
 Special thanks to *Jean-Baptiste Tanguy* (Sorbonne Université) for bringing this visualization to my attention. Thanks also to *Steve Mutuvi* (Multimedia University of Kenya, La Rochelle University), Caroline Koudoro-Parfait and Richy Nutt (Sorbonne Université) for testing the different versions and bringing useful ideas and comments.


### Usage

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
