set cubefile [lindex $argv 0]
set outfname [lindex $argv 1]

color Display Background white
display projection orthographic

mol new $cubefile
mol modstyle 0 0 Licorice 0.3 27.0 27.0
mol addrep 0
mol modstyle 1 0 Isosurface 0.01 0 0 0 1 1
mol modcolor 1 0 ColorID 0
mol addrep 0
mol modstyle 2 0 Isosurface -0.01 0 0 0 1 1
mol modcolor 2 0 ColorID 1

display resetview
rotate z by 90
rotate x by -90

display nearclip set 0.01
axes location off
render TachyonInternal $outfname
quit

