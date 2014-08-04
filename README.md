EPA-classifier
==============

(1) What does EPA-classifier do?

    Classification of sequences into taxonomic ranks using EPA.
    
(2) Input.

    Reference database in EPA-classifier json format, and query 
    sequences in FASTA format.
    
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
    
(5) Basic usage

    Classification:
    
        ./epa_classifier.py -r example/prokaryotes.json -q example/query.fa -T 4
    
    Will classify sequences in query.fa using the standard reference for
    prokaryotes (Bacteria + Archaea) using 4 threads/CPUs. Classification
    results will be saved to example/query.fa.assignment.txt, 
    see (3) for output format.
    
    Training:
    
        ./epa_trainer.py -t example/training_tax.txt -s example/training_seq.fa -r example/myref.json -T 4
        
    Build a reference database from MSA (training_seq.fa) and taxonomic 
    annotations (training_tax.txt). The output file, myref.json, can be 
    subsequently used as reference for the classification (s. above).
    
    For advanced usage, please refer to online help (-h switch) and 
    HOWTO file in this directory.

(6) Version

    EPA-classifier is still under development.
    
    
    
    EPA_PTP.py (Broken now, I will fix it when I have time)This is a pipeline
               that uses evolutionary placement algorithm (EPA) and PTP to 
               count species number when reference data is available. For 
               details of EPA, please read this paper: Performance, accuracy, 
               and Web server for evolutionary placement of short sequence 
               reads under maximum likelihood. Systematic biology, 60(3), 
               291–302. The pipeline will first run USEARCH to remove the 
               chimera sequences, then it will use EPA to place the query 
               reads to optimal position on the reference tree inferred 
               from the reference alignment. PTP will then be applied to 
               the reads been placed on each branch, with a fixed speciation
               rate inferred from the reference data.  
               Similar analysis used in bacterial metagenomics studies are 
               called OTU-picking. For discussions about OTU-picking and EPA_PTP 
               species counting, please have a look at our paper. 
               Type ./EPA_PTP.py for help and instructions.
