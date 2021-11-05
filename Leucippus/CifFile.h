#pragma once
/*
* Rachel Alcraft, 04/11/2021
* Generically loads a file fo the mmcif format which is used for pdb files and strcture factors
*/

#include <string>
#include <vector>
#include <map>

using namespace std;

class CifFile
{
private:
    //A dictionary of item name to item, all 1:! the items in the file are aggregated in the same place
    map<string, string> _nonLoopElements;    
    //Whenever we come across a loop, here we have: a dictionar wherer the key is the first item in the header list before the "dot"
    //The keys of all the dictionaries are the headers
    // The lists belonging to each header are the data which are space delim
    map<string, map<string, vector<string> > > _LoopElements;
    /* AN EXAMPLE of the data
    *   #        
        _citation.pdbx_database_id_PubMed 10737790
        _citation.pdbx_database_id_DOI 10.1073/pnas.97.7.3171
        #
        loop_
        _citation_author.citation_id
        _citation_author.name
        _citation_author.ordinal
        primary 'Jelsch, C.' 1
        primary 'Teeter, M.M.' 2
        primary 'Lamzin, V.' 3        
        #
        Would look like
        _nonLoopElements[citation.pdbx_database_id_PubMed] = 10737790
        _nonLoopElements[citation.pdbx_database_id_DOI] = 10.1073/pnas.97.7.3171

        _LoopElements[citation_author][citation_id] = [primary,primary,primary]
        _LoopElements[citation_author][name] = ['Jelsch, C.','Teeter, M.M.','Lamzin, V.']
        _LoopElements[citation_author][ordinal] = [1,2,3]

    */

public:
    CifFile (string fileName);    
};