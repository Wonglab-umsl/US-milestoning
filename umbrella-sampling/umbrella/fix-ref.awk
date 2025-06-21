#set the occupancy to 1 for the drug, 0 for everything else
#set beta to 1 for the protein's backbone and 0 for everthing else
function trim(field){
   gsub(/^ +| +$/,"", field); 
   return field
   }
/ATOM/ {
	if ((substr($0,18,3)=="10N") || (substr($0,18,3)=="CPD")) occ=1.0; else occ=0.0;
	aname=trim(substr($0,13,4));
	if ((occ==0.0) && ((aname=="N") || (aname=="CA") || (aname=="C"))) beta=1.0; else beta=0.0;
	line=substr($0,1,54) sprintf("%6.2f%6.2f",occ,beta) substr($0,72,9);
	print line;
}

($0!~/ATOM/) {print $0}
