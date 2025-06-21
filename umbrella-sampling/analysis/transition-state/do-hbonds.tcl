#set vmdhome $env(VMDHOME)
set vmdhome "/share/apps/common/VMD_1.9.4/"
#play "/net/home/spiritij/vmd/plugins/noarch/tcl/hbonds1.2/hbonds.tcl"
play "$vmdhome/lib/vmd/plugins/noarch/tcl/hbonds1.2/hbonds.tcl"
#quit
set name [lindex $argv 0]
set dcdname [lindex $argv 1]
set shortname [lindex $argv 2]
#set scratch [lindex $argv 3]
#mol new $psfname-ref.pdb
mol new prmtop/$name.prmtop
mol addfile $dcdname waitfor all
set nf [molinfo top get numframes]




set fname1 "an/hbonds-time-series-$shortname"
set fname2 "an/hbonds-counts-$shortname"
#puts "marker"
set prot [atomselect top "protein"]
set drug [atomselect top "resname CPD"]
hbonds -sel1 $prot -sel2 $drug -upsel no -frames all -dist 3.5 -ang 20 -type unique -writefile yes -outfile $fname1 -detailout $fname2 -plot no
quit
