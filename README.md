EPA-classifier
==============

(1) What does EPA-classifier do?

    Classification of sequences into taxonomic ranks using EPA.
    
(2) Input.

    Reference database in EPA-classifier json format, and query 
    sequences in fasta format.
    
(3) Output.

    query_sequence_name \tab taxonomic_ranks \tab confidence \tab remark
    remark can be one of the three:
    ?: means can not determine the taxonomic ranks
    *: means the query sequence can not be assigned to genus level, and 
       is likely to come from a genus that does not exist in the 
       reference database.
    o: means the if the query sequence is not assigned to genus level is
       likely due to the incomplete information in the reference.  
    

(4) Download and install.

    Linux 64bit (due to precompiled RAxML, muscle, and hmmer binary)
    python>=2.5 (python 3 is not supported)
    

(5) Version

    EPA-classifier is still under development.
