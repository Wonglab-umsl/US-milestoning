set psfname [lindex $argv 0]
set dcdname [lindex $argv 1]
set refname [lindex $argv 2]
set outfname [lindex $argv 3]
set com0 [lrange $argv 4 6]
#set fname "../data/$name-force-250-window-$window-ref.pdb"
mol new $refname waitfor all
mol new $psfname
mol addfile $dcdname  waitfor all
#mol addfile $dcdname first 100000 last 199999 waitfor all
mol list
set nf [molinfo top get numframes]
set bb [atomselect top "protein and (name N or name CA or name C)"]
set bbref [atomselect 0 "protein and (name N or name CA or name C)"]
set drugall [atomselect top "resname CPD"]
set drugallref [atomselect 0 "resname CPD"]
set all [atomselect top "all"]
#set com0 [measure center $drugallref weight mass]
set output [open "$outfname" "w"]
for {set i 0} {$i < $nf} {incr i} {
	#puts $i
	if { ( $i % 1000 ) == 0 } {
		puts $i
	}
	animate goto $i
	set matrix [measure fit $bb $bbref]
	$all move $matrix
	#puts $matrix
	set bbrmsd [measure rmsd $bb $bbref]
	#set drugrmsd [measure rmsd $drug $drugref]
	set com [measure center $drugall weight mass]
	#if { $i==0 } {set com0 $com}
	set d [veclength [vecsub $com $com0]]
	set time [expr ($i + 1) * .001]
	#to be compatible with the input statement in wham3d.f90 and wham3d-1d.f90
	#omit first 10 data points to try to get a little better equilibrium and avoid artifacts
	#if {$i > 1000} {
		puts $output "$time $com $d" 
	#}
	#puts $output "$time $bbrmsd $com $d"
	#if { $i > 5} { quit }
}
close $output
quit
