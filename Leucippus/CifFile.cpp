
#include <fstream>

#include "CifFile.h"
#include <iostream>

CifFile::CifFile(string fileName)
{
    /* AN EXAMPLE of the data
       _citation.pdbx_database_id_PubMed 10737790
       _citation.pdbx_database_id_DOI 10.1073/pnas.97.7.3171

       loop_
       _citation_author.citation_id
       _citation_author.name
       _citation_author.ordinal
       primary 'Jelsch, C.' 1
       primary 'Teeter, M.M.' 2
       primary 'Lamzin, V.' 3

       Would look like
       _nonLoopElements[citation][pdbx_database_id_PubMed] = 10737790
       _nonLoopElements[citation][pdbx_database_id_DOI] = 10.1073/pnas.97.7.3171

       _LoopElements[citation_author][citation_id] = [primary,primary,primary]
       _LoopElements[citation_author][name] = ['Jelsch, C.','Teeter, M.M.','Lamzin, V.']
       _LoopElements[citation_author][ordinal] = [1,2,3]

   */
   //map<string, map<string, string> > _nonLoopElements;    
   //map<string, map<string, vector<string> > > _LoopElements;

    ifstream myfile(fileName.c_str());

    if (myfile.is_open())
    {
        bool on_loop = false;
        bool on_item = false;
        string item_name = "";

        string line = "";
        vector<string> headers;
        while (getline(myfile, line))
        {
            //cout << '\n';
            if (line[0] == '#' && on_item)//then we are closing our current looop or item
            {
                on_item = false;
                on_loop = false;
                item_name = "";
            }
            else if (line[0] == 'l' && !on_item)//we are entering a loop set of data
            {
                on_loop = true;
                on_item = true;
                headers.clear();
            }
            else if (line[0] == '_' && !on_item)//then we are starting a non loop set of data
            {
                on_item = true;
                on_loop = false;
                string id = getIdFromLine(line);
                string nam = getNameFromLine(line);
                string val = getValueFromLine(line);
                NonLoopElements.insert(pair<string, map<string, string> >(id, map<string, string>()));
                NonLoopElements[id].insert(pair<string, string>(nam, val));
            }
            else if (line[0] == '_' && on_item && !on_loop)//then we are continuiing a non loop set of data
            {
                on_item = true;
                on_loop = false;
                string id = getIdFromLine(line);
                string nam = getNameFromLine(line);
                string val = getValueFromLine(line);
                NonLoopElements[id].insert(pair<string, string>(nam, val));
            }
            else if (line[0] == '_' && on_item && on_loop)//then we are adding headers to a loop set
            {
                on_item = true;
                on_loop = true;
                string id = getIdFromLine(line);
                item_name = id;
                string nam = getNameFromLine(line);
                map<string, map<string, vector<string> > >::iterator iter;
                iter = LoopElements.find(id);
                if (iter == LoopElements.end())
                {
                    LoopElements.insert(pair<string, map<string, vector<string> > >(id, map<string, vector<string> >()));
                }
                LoopElements[id].insert(pair<string, vector<string> >(nam, vector<string>()));
                headers.push_back(nam);
            }
            else if (line[0] != '#' && line[0] != '_' && on_item && on_loop)//then we are adding a row to a loop set
            {
                vector<string> vals = getValuesFromLine(line);
                for (unsigned int i = 0; i < vals.size(); ++i)
                {
                    string header = headers[i];
                    string val = vals[i];
                    LoopElements[item_name][header].push_back(val);
                }

            }

        }
        myfile.close();

    }
}



// ~~~~~~~~~~~~~~~~~~~~~~~~ Private helper funtions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
string CifFile::getIdFromLine(string line)
{
    //Data looks like either of these 2
    // _id.name value
    //_citation.pdbx_database_id_PubMed 10737790
    //_citation_author.citation_id
    size_t pos = line.find(".");
    string id = line.substr(0,pos);
    //cout << id << " ";
    return id;
}
string CifFile::getNameFromLine(string line)
{
    //Data looks like either of these 2
    // // _id.name value
    //_citation.pdbx_database_id_PubMed 10737790
    //_citation_author.citation_id
    size_t pos = line.find(".");
    string rhs = line.substr(pos+1);
    pos = rhs.find(" ");
    if (pos > 0)    
        rhs = rhs.substr(0, pos);
    
    //cout << rhs << " ";
    return rhs;
}
string CifFile::getValueFromLine(string line)
{
    //_citation.pdbx_database_id_PubMed 10737790
    // _id.name value
    size_t pos = line.find(".");
    string rhs = line.substr(pos + 1);
    
    pos = rhs.find(" ");    
    if (pos > 0)
    {
        rhs = rhs.substr(pos + 1);
        pos = rhs.find(" ");
    }
    
    while (pos == 0)
    {
        rhs = rhs.substr(pos + 1);
        pos = rhs.find(" ");
    }
    if (pos > 0)
        rhs = rhs.substr(0,pos);

    //cout << rhs << " ";
    return rhs;
}
vector<string> CifFile::getValuesFromLine(string line)
{//primary 'Jelsch, C.' 1
    //space of single quite delim if there are spaces (could do this better I know).
    vector<string> vals;
    string curr_val = "";

    size_t pos_space = line.find(" ");
    size_t pos_quote = line.find("'");
        
    while (line.size() > 0)
    {
        pos_space = line.find(" ");
        while (pos_space == 0 && line.size() > 0) //trim off the spaces
        {
            line = line.substr(1);
            pos_space = line.find(" ");
        }
        if (line.size() > 0)
        {

            pos_space = line.find(" ");
            pos_quote = line.find("'");

            if (pos_quote == 0)
            {//then the quite comes before the space so we have an open quote word so it ends at next quote
                size_t pos_quote_end = line.find("'", pos_quote + 1);
                string val = line.substr(pos_quote+1, pos_quote_end);
                vals.push_back(val);
                line = line.substr(pos_quote_end+1);
                //cout << val << " ";
            }
            else
            {//we have a spaces word so it ends at next sapce unless that gives us a word of only spaces
                string val = line.substr(0,pos_space);
                vals.push_back(val);
                //cout << val << " ";                
                line = line.substr(pos_space);
            }
        }        
    }
    return vals;

}

