'''
Created on Nov 23, 2009
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os
# ropy import
import ropy
# genepy import
from setup import connectionFactory

def main():
    query = ropy.query.QueryProcessor(connection=connectionFactory)
    query.setSQLFilePath(os.path.dirname(__file__) + "/sql/")
    query.addQueryFromFile("organism_query", "get_all_organisms_with_polyseq.sql")
    organism_rows = query.runQuery("organism_query")
    print "%s organism sequences from geneDB." % len(organism_rows)
    
    query.addQueryFromFile("fasta_query", "get_fasta_polyseq_for_organism.sql")
    
    for organism in organism_rows:
        organism_name = organism[1]
        organism_id = organism[0]
        if organism_name == "dummy":
            continue
        
        # Dump sequence of each organism into a fasta file
        print "Extracting %s..." % organism_name
        fasta_rows = query.runQuery("fasta_query", (organism_id, ))
        for row in fasta_rows:
            if not row[0] == None:
                print row[0]

if __name__ == '__main__':
    main()