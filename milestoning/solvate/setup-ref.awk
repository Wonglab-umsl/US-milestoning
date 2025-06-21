/ATOM/ {
	occ=0.0;
	resname=substr($0,18,3);
	aname=substr($0,13,4);
	gsub(" ", "", aname);
	if ((resname!="WAT") && (resname!="Na+") && (resname!="Cl-") && (substr(aname,1,1)!="H")) occ=1.0;
	line=substr($0,1,60) sprintf("%6.2f",occ) substr($0,67,20);
	print line;
}

