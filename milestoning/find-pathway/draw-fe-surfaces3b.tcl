#set dxfile [lindex $argv 0]
#set outfnamefmt [lindex $argv 1]
set name [lindex $argv 0]
#set force [lindex $argv 1]
#set dxfile "$name-explicit-wham-fe.dx"
set dxfile "$name-explicit-wham-fe-rel-to-first.dx"
#set outfnamefmt "images/$name-explicit-wham-fe-%.1f.tga"
set outfnamefmt "images/$name-explicit-wham-fe-%.1f-rel-to-first.tga"
set outfnamefmt2 "images/$name-explicit-wham-fe-%.1f-rel-to-first-side.tga"
set festart [lindex $argv 1]
set feinc [lindex $argv 2]
set feend [lindex $argv 3]
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
	
mol addfile $dxfile type dx
mol modselect 0 0 "protein and noh"
mol modstyle 0 0 NewCartoon 0.3 20.0 4.1 0
mol addrep 0
mol modselect 1 0 "protein and (resid 429 or resid 431 or resid 509 or resid 512) and noh"
#"protein and (resid 426 or resid 428 or resid 501 or resid 502 or resid 506 or resid 509 or resid 553 or resid 567) and noh"
mol modstyle 1 0 Licorice 0.3 24.0 24.0
mol addrep 0
set rep 2
mol modcolor $rep 0 ColorID 6
for {set fe $festart} {$fe <= $feend} {set fe [expr $fe + $feinc]} {
	display resetview
	display nearclip set 0.01
	display depthcue off
	axes location off
	if { "$name" == "pyk2-cpd1-cgenff" } {
		scale by 7.00
		rotate x by 180
		rotate y by 180
		rotate z by -30
		rotate x by -15
		rotate y by -30
		translate by -2.0 -0.5 0.0
	}  elseif { "$name" == "pyk2-cpd1-amber" } {
	        scale by 6.00
        	rotate x by 180
        	#rotate y by 180
        	#rotate z by -30
		rotate y by 90
		rotate x by 45
		rotate y by 30
		#translate by -2.0 0.0 0.0
       		translate by -0.5 -0.5 -1.0
	} elseif { "$name" == "pyk2-cpd10-amber" } {
        	scale by 5.00
        	#rotate y by 90
        	#rotate z by 60
        	#rotate x by -30
        	#rotate y by 15
        	#translate by 0.0 -1.0 0.0
		rotate y by 120
	        rotate z by 20
		translate by -0.5 -0.25 -1.0
		puts "done resetting view"
	} elseif { "$name" == "pyk2-cpd2-amber" } {
		scale by 5.00
                rotate x by 180
                #rotate y by 180
                #rotate z by -30
                rotate y by 90
                rotate x by 45
	} elseif { "$name" == "pyk2-cpd6-amber" } {
        	scale by 6.00
		rotate x by 180
		#rotate y by 180
		#rotate z by -30
		rotate y by 90
		rotate x by 45
	} else {
	}
	puts $fe
	mol modstyle $rep 0 Isosurface $fe 0 0 0 1 1
	set fname [format $outfnamefmt $fe]
	render TachyonInternal $fname	
	display resetview
        scale by 5.00
	if { "$name" == "pyk2-cpd10-amber" } {
		rotate y by 180
		rotate z by 30
		rotate x by 15
	} else {
        	rotate x by 180
        	rotate y by 180
        	rotate z by -30
	}
	translate by -0.5 -0.5 -2.0
	set fname [format $outfnamefmt2 $fe]
        render TachyonInternal $fname
}


quit
