OlfactorySystem
===============

Olfactory system simulation 


  For the distribution of Olfactory Receptor affinities:
    The data in the Haddad_data folder is from http://www.nature.com/nmeth/journal/v5/n5/extref/nmeth.1197-S3.xls [Haddad 2008 "A metric for odorant comparison", Nature Methods] 
    The script cluster_odorant_space.py computes the distance between the
    real-world data and the virtual olfactory receptors (=centroids after
    k-means clustering) for many trials and for various numbers of centroids.
    
    Distances between ORs and odorants is pooled by
    average_OR_affinity_distributions.py over many trials.
    It tries to fit a distribution to the data and writes the fit parameters
    to a file, which can be displayed by plot_OR_placement_fit_params.py.
    

     
