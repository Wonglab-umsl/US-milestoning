set refname [lindex $argv 0]
set psf [lindex $argv 1]
set traj [lindex $argv 2]
set iframe [lindex $argv 3]
set outfname [lindex $argv 4]

mol new $refname
mol new $psf
mol addfile $traj first $iframe last $iframe waitfor all

set bb [atomselect top "protein and (name N or name CA or name C)"]
$bb num
set bbref [atomselect 0 "protein and (name N or name CA or name C)"]
$bbref num
#quit
set all [atomselect top "all"]
set prot [atomselect top "protein"]
mol list
set matrix [measure fit $bb $bbref]
$all move $matrix
set prot [atomselect top "protein or resname CPD"]
set drugall [atomselect top "resname CPD"]
set com [measure center $drugall weight mass]
puts "Old drug COM: $com"
#need to center the complex on the origin for solvation
#set center [measure center $prot weight mass]
#set trans [vecinvert $center]
#$all moveby $trans

#set newcom [measure center $drugall weight mass]
#puts "New drug COM: $newcom"
$prot writepdb $outfname

quit

