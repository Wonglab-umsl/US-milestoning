#set dxfile [lindex $argv 0]
#set outfnamefmt [lindex $argv 1]
set name [lindex $argv 0]
set dxfile "$name-explicit-wham-fe-rel-to-first.dx"
set fe1 [lindex $argv 1]
set fe2 [lindex $argv 2]
set fname [format "images/$name-explicit-wham-fe-contours-%.1f-%.1f.tga" $fe1 $fe2]
set fname2 [format "images/$name-explicit-wham-fe-contours-%.1f-%.1f-side.tga" $fe1 $fe2]

color Display Background white
display projection orthographic
#set name fak-cpd32-cgenff
#if { "$name" == "fak-cpd32-cgenff" } {
#	set force 250
#} else {
#	set force 350
#}
if { ( $name == "pyk2-cpd1-amber" ) || ( $name == "pyk2-cpd10-amber" ) } {
	mol new ../ref/$name-tramd-window-1-ref.pdb
} else {
	mol new ../ref/$name-window-1-solvated-ref.pdb
}
mol addfile $dxfile type dx
mol modselect 0 0 "protein and noh"
mol modstyle 0 0 NewCartoon 0.3 20.0 4.1 0
mol modcolor 0 0 ColorID 6
mol addrep 0
mol modstyle 1 0 Isosurface $fe1 0 0 0 1 1
mol modcolor 1 0 ColorID 0
mol addrep 0
mol modstyle 2 0 Isosurface $fe2 0 0 0 1 1
mol modcolor 2 0 ColorID 1
mol modmaterial 2 0 Transparent
#mol addrep 0
#mol modselect 3 0 "protein and (resid 14 or resid 90 or resid 93 or resid 98) and noh"
#mol modstyle 3 0 Licorice 0.3 24.0 24.0

mol list 0

display resetview
if { "$name" == "fak-cpd32-cgenff" } {
	scale by 7.00
	rotate y by -90
	rotate z by 30
	rotate x by 30
	rotate y by -15
	translate by -2.0 0.5 -1.0
} elseif { "$name" == "fak-cpd2-cgenff"} {
	scale by 7.00
	rotate y by -120
	rotate x by -15
	rotate z by 15
	translate by -2.0 0.0 -1.0
} elseif { "$name" == "fak-cpd41-cgenff"} {
	scale by 7.00
	rotate y by -90
	rotate z by 5
	translate by -2.5 0.0 -1.0
} elseif { "$name" == "pyk2-cpd1-amber" } {
        scale by 7.00
        rotate x by 180
        rotate y by 90
        rotate x by 45
        #rotate y by 30
        translate by 0.0 -0.15 0.0
	#scale by 7.00
	#rotate x by 180
	#rotate y by 180
	#rotate z by -30
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
        } elseif { "$name" == "pyk2-cpd2-amber" } {
                scale by 4.00
                rotate x by 180
                #rotate y by 180
                #rotate z by -30
                rotate y by 90
                rotate x by 45
		rotate y by 30
        } elseif { "$name" == "pyk2-cpd6-amber" } {
                scale by 4.00
                rotate x by 180
                #rotate y by 180
                #rotate z by -30
                rotate y by 90
                rotate x by 40
        } elseif { "$name" == "pyk2-cpd8-amber" } {
                scale by 4.00
                rotate x by 180
                #rotate y by 180
                #rotate z by -30
                rotate y by 90
                rotate x by 40
} else {
}


display nearclip set 0.01
axes location off
render TachyonInternal $fname

display resetview
        scale by 4.00
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
        render TachyonInternal $fname2

quit

