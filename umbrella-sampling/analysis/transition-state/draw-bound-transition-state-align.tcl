proc align_traj {m mref alignsel} {
	set ref [atomselect $mref $alignsel]
	mol top $m
	set nframes [molinfo $m get numframes]
	#rms align trajectory -- mol. 0 is the trajecotry, 1 is reference
	for {set i 0} {$i < $nframes} {incr i} {
        	#puts "RMS alignment mol $m frame $i"
        	animate goto $i
        	set meas [atomselect $m $alignsel]
        	set trans [measure fit $meas $ref]
        	set move [atomselect $m "all"]
        	$move move $trans
	}
}


set name [lindex $argv 0]
set alignsel $env(alignsel)
color Display Background white
display projection orthographic
set dir /mnt/stor/ceph/scratch/jmsc87/pyk2/umbrella-amber/transition-state

#for some reason this script crashes if the full trajectory is loaded.  
mol new prmtop/$name.prmtop
mol addfile $dir/$name-bound-state-short.dcd step 100 waitfor all
mol new prmtop/$name.prmtop
mol addfile $dir/$name-transition-state-short.dcd step 100 waitfor all
mol new pdb/$name-ref.pdb

mol list
align_traj 0 2 $alignsel
align_traj 1 2 $alignsel

mol off 2
mol modcolor 0 0 ColorID 0
mol modstyle 0 0 NewCartoon 0.3 20 4.1 0
mol drawframes 0 0 0:1:10
mol addrep 0
mol modcolor 1 0 ColorID 0
mol modstyle 1 0 Licorice 0.3 30 30
mol modselect 1 0 "noh and resname CPD"
mol drawframes 0 1 0:1:10

mol modcolor 0 1 ColorID 1
mol modstyle 0 1 NewCartoon 0.3 20 4.1 0
mol drawframes 1 0 0:1:10
mol addrep 1
mol modcolor 1 1 ColorID 1
mol modstyle 1 1 Licorice 0.3 30 30
mol modselect 1 1 "noh and resname CPD"
mol drawframes 1 1 0:1:10

puts "marker 1"
flush stdout
#display resetview
puts "marker 1.1"
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
	puts "marker 1.2"
	flush stdout
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
		puts "marker 1.2"
		flush stdout
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
		puts "marker 1.2"
		flush stdout
                scale by 1.50
                rotate x by 180
                #rotate y by 180
                #rotate z by -30
                rotate y by 90
                rotate x by 40

} else {
}
puts "marker 2"
display nearclip set 0.01
axes location off

set outfname "images/$name-bound-transition-state.tga"
render TachyonInternal $outfname

rotate y by 180
set outfname "images/$name-bound-transition-state2.tga"
render TachyonInternal $outfname



quit
	

