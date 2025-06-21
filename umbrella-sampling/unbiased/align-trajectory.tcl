set psfname [lindex $argv 0]
set dcdname [lindex $argv 1]
set refname [lindex $argv 2]
set outfname [lindex $argv 3]
#set com0 [lrange $argv 4 6]
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
#set output [open "$outfname" "w"]
for {set i 0} {$i < $nf} {incr i} {
	#puts $i
	if { ( $i % 1000 ) == 0 } {
		puts $i
	}
	animate goto $i
	set matrix [measure fit $bb $bbref]
	$all move $matrix
}
#close $output
puts "writing trajectory"
set nfm1 [expr $nf - 1]
animate write dcd $outfname beg 0 end $nfm1 skip 1 waitfor all top
puts "trajectory written"
quit
