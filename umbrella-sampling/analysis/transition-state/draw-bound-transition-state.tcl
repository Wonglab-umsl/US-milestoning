set name [lindex $argv 0]
color Display Background white
display projection orthographic
set dir /lustre/scratch/jmsc87/pyk2/umbrella-amber/transition-state

mol new prmtop/$name.prmtop
mol addfile $dir/$name-bound-state-short.dcd waitfor all
mol new prmtop/$name.prmtop
mol addfile $dir/$name-transition-state-short.dcd waitfor all

mol modcolor 0 0 ColorID 0
mol modstyle 0 0 NewCartoon
mol drawframes 0 0 0:10:1000
mol modcolor 0 1 ColorID 1
mol modstyle 0 1 NewCartoon
mol drawframes 1 0 0:10:1000



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
        scale by 1.50
        rotate x by 180
        #rotate y by 180
        #rotate z by -30
        rotate y by 90
        rotate x by 45
        rotate y by 30
        translate by 0.0 -0.15 0.0
} elseif { "$name" == "pyk2-cpd10-amber" } {
        scale by 2.00
        rotate y by 90
        rotate z by 60
        rotate x by -30
        rotate y by 15
        #translate by 0.0 -1.0 0.0
} elseif { "$name" == "pyk2-cpd2-amber" } {
                scale by 1.50
                rotate x by 180
                #rotate y by 180
                #rotate z by -30
                rotate y by 90
                rotate x by 45
                rotate y by 30
        } elseif { "$name" == "pyk2-cpd6-amber" } {
                scale by 1.50
                rotate x by 180
                #rotate y by 180
                #rotate z by -30
                rotate y by 90
                rotate x by 40
        } elseif { "$name" == "pyk2-cpd8-amber" } {
                scale by 1.50
                rotate x by 180
                #rotate y by 180
                #rotate z by -30
                rotate y by 90
                rotate x by 40

} else {
}

display nearclip set 0.01
axes location off

set outfname "images/$name-bound-transition-state.tga"
render TachyonInternal $outfname

rotate y by 180
set outfname "images/$name-bound-transition-state2.tga"
render TachyonInternal $outfname



quit
	

