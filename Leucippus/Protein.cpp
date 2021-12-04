#include "Protein.h"

Protein::Protein(map<string, vector<string>> atomList)
{
	vector<string> xs = atomList["Cartn_x"];
	vector<string> ys = atomList["Cartn_y"];
	vector<string> zs = atomList["Cartn_z"];
	vector<string> ats = atomList["label_alt_id"]; //EG CA, CB OG1
	vector<string> ts = atomList["group_PDB"]; //eg HETATM or ATOM
	vector<string> occs = atomList["occupancy"];
	vector<string> bfs = atomList["B_iso_or_equiv"];
	MinX = 10000;
	MaxX = -10000;
	MinY = 10000;
	MaxY = -10000;
	MinZ = 10000;
	MaxZ = -10000;
	for (unsigned int v = 0; v < xs.size(); ++v)
	{
		double x = atof(xs[v].c_str());
		double y = atof(ys[v].c_str());
		double z = atof(zs[v].c_str());
		string at = ats[v];
		string t = ts[v];
		string occ = occs[v];
		string bf = bfs[v];
		Atom *atm = new Atom(x,y,z,atof(bf.c_str()), atof(occ.c_str()), at, t);
		Atoms.push_back(atm);
		MinX = min(MinX, x);
		MaxX = max(MaxX, x);
		MinY = min(MinY, y);
		MaxY = max(MaxY, y);
		MinZ = min(MinZ, z);
		MaxZ = max(MaxZ, z);
	}
}
