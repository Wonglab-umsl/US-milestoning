set name [lindex $argv 0]
#set force [lindex $argv 1]
set nwindow [lindex $argv 1]
set outfname $name-tramd-windows-rel-to-first
set output [open "$outfname" "w"]
mol new data/$name-tramd-window-1-solvated-ref.pdb

for {set i 1} {$i <= $nwindow} {incr i} {
	mol new data/$name-tramd-window-$i-solvated-ref.pdb
	set bb [atomselect top "protein and (name N or name CA or name C)"]
	set bbref [atomselect 0 "protein and (name N or name CA or name C)"]
	set all [atomselect top "all"]
        set matrix [measure fit $bb $bbref]
        $all move $matrix
        set bbrmsd [measure rmsd $bb $bbref]
	set sel [atomselect top "resname CPD"]
	set com [measure center $sel weight mass]
	puts $output "$i $com $bbrmsd"
	mol delete top
}
