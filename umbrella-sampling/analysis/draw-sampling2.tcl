#set dxfile [lindex $argv 0]
#set outfnamefmt [lindex $argv 1]
set name [lindex $argv 0]
set outfname "images/$name-sampling-rel-to-first.tga"
set outfname2 "images/$name-sampling-rel-to-first-side.tga"
set nwindow [lindex $argv 1]
set freq [lindex $argv 2]
set windowfile [lindex $argv 3]
color Display Background white
display projection orthographic
#set name fak-cpd32-cgenff
#if { "$name" == "fak-cpd32-cgenff" } {
#	set force 250
#} else {
#	set force 350
#}
if { ( "$name" == "pyk2-cpd1-amber" ) || ( "$name" == "pyk2-cpd10-amber" ) } { 
	set refpdb  "../ref/$name-tramd-window-1-ref.pdb"
} else {
	set refpdb  "../ref/$name-window-1-solvated-ref.pdb"
}
mol new $refpdb
	
#mol addfile $dxfile type dx
mol modselect 0 0 "protein and noh"
mol modstyle 0 0 NewCartoon 0.3 20.0 4.1 0
mol modcolor 0 0 ColorID 6

set nm1 [expr $nwindow - 1]
mol new
display update off
#draw color silver
#color scale method rgb
set colorlist [list red red3 orange2 orange orange3 yellow yellow2 yellow3 green green2 lime green3\
        cyan cyan2 cyan3 blue2 blue blue3 violet violet2 purple magenta magenta2 red2]
set ncolor [llength $colorlist]
set input2 [open "$windowfile" "r"]
for {set i 1} {$i <= $nwindow} {incr i} {
        #set cidx [expr int( 1023 * ( ( $i - 1.0 ) / ($nm1) ) + 33) ]
	if { ( "$name" == "pyk2-cpd1-amber" ) || ( "$name" == "pyk2-cpd10-amber" ) } {
		set infname "data-rel-to-first/$name-tramd-com-$i"
	} else {
		set infname "data-rel-to-first/$name-com-$i"
	}
	set color [lindex $colorlist [expr ($i - 1) % $ncolor]]
        puts "$i $color"
        draw color "$color"
	gets $input2 line
	scan $line "%d %g %g %g" ii x y z
	if { $ii == $i } {
		set windowcenter [list $x $y $z]
		draw sphere $windowcenter radius 0.2 resolution 30
	}
	set input [open $infname "r"]
	set counter 1
	while { [gets $input line] >= 0 } {
		set r [expr $counter % $freq]
		if { $r == 0 } {
			scan $line "%g %g %g %g" time x y z
			#puts "$time $x $y $z"
			if { $time >= 2.0 } {
       		 		set com [list $x $y $z]
				draw point $com
			}
		}
		set counter [expr $counter + 1]
	}
	close $input
}

display update on

display resetview
if { "$name" == "fak-cpd32-cgenff" } {
	#scale by 7.00
	rotate y by -90
	rotate z by 30
	rotate x by 30
	rotate y by -15
	#translate by -2.0 0.5 -1.0
} elseif { "$name" == "fak-cpd2-cgenff"} {
	#scale by 7.00
	rotate y by -120
	rotate x by -15
	rotate z by 15
	#translate by -2.0 0.0 -1.0
} elseif { "$name" == "fak-cpd41-cgenff"} {
        #scale by 7.00
        rotate y by -90
        rotate z by 5
        translate by -2.5 0.0 -1.0
} elseif { "$name" == "pyk2-cpd1-cgenff" } {
        #scale by 7.00
        rotate x by 180
        rotate y by 180
        rotate z by -30
        #translate by -2.0 0.0 0.0
} elseif { "$name" == "pyk2-cpd1-amber" } {
        scale by 2.00
        rotate x by 180
        #rotate y by 180
        #rotate z by -30
        rotate y by 90
        rotate x by 45
        rotate y by 30
        #translate by -2.0 0.0 0.0
} elseif { "$name" == "pyk2-cpd10-amber" } {
	scale by 2.00
	#rotate y by 90
	#rotate z by 60
	#rotate x by -30
	#rotate y by 15
	#translate by 0.0 -1.0 0.0
	rotate y by 120
	rotate z by 20
} elseif { "$name" == "pyk2-cpd2-amber" } {
        scale by 1.00
        rotate x by 180
        #rotate y by 180
        #rotate z by -30
        rotate y by 90
        rotate x by 45
} elseif { "$name" == "pyk2-cpd6-amber" } {
        scale by 1.50
        rotate x by 180
        #rotate y by 180
        #rotate z by -30
        rotate y by 90
        rotate x by 45
} elseif { "$name" == "pyk2-cpd8-amber" } {
        scale by 1.50
        rotate x by 180
        #rotate y by 180
        #rotate z by -30
        rotate y by 90
        rotate x by 45
} else {
}

display nearclip set 0.01
axes location off

render TachyonInternal $outfname
display resetview
scale by 1.50
if { "$name" == "pyk2-cpd10-amber" } {
	rotate y by 180
	rotate z by 30
	rotate x by 15
} else {
	rotate x by 180
	rotate y by 180
	rotate z by -30
}
render TachyonInternal $outfname2



quit
