#!/usr/bin/env bash
# Usage: ./airfoil2stl.sh input.dat output.stl [span]
set -euo pipefail
infile="${1:-}"; outfile="${2:-}"; span="${3:-0.1}"
if [[ -z "${infile}" || -z "${outfile}" ]]; then
  echo "Usage: $0 input.dat output.stl [span]"; exit 1
fi
export LC_ALL=C

awk -v span="$span" -v out="$outfile" '
function cross(ax,ay,az, bx,by,bz,   cx,cy,cz){ cx=ay*bz-az*by; cy=az*bx-ax*bz; cz=ax*by-ay*bx; return cx" "cy" "cz }
function norm(x,y,z,   m){ m=sqrt(x*x+y*y+z*z)+1e-30; return x/m" "y/m" "z/m }
function facet(x1,y1,z1, x2,y2,z2, x3,y3,z3,   ux,uy,uz, vx,vy,vz, v, nn) {
  ux = x2-x1; uy = y2-y1; uz = z2-z1;
  vx = x3-x1; vy = y3-y1; vz = z3-z1;
  split(cross(ux,uy,uz, vx,vy,vz), v,  " ");
  split(norm(v[1],v[2],v[3]),      nn, " ");
  printf("  facet normal %.8g %.8g %.8g\n", nn[1], nn[2], nn[3]) >> out;
  print  "    outer loop" >> out;
  printf("      vertex %.10g %.10g %.10g\n", x1,y1,z1) >> out;
  printf("      vertex %.10g %.10g %.10g\n", x2,y2,z2) >> out;
  printf("      vertex %.10g %.10g %.10g\n", x3,y3,z3) >> out;
  print  "    endloop" >> out;
  print  "  endfacet" >> out;
}

BEGIN{ OFS=" "; half = span/2.0; cnt = 0; }

{
  gsub(/\r$/, ""); gsub(/,/, " "); gsub(/\t+/, " ");
  if ($1 ~ /^([+-]?([0-9]*[.])?[0-9]+([eE][+-]?[0-9]+)?)$/ && \
      $2 ~ /^([+-]?([0-9]*[.])?[0-9]+([eE][+-]?[0-9]+)?)$/) {
    x[cnt] = +$1; y[cnt] = +$2; cnt++;
  }
}

END{
  if (cnt < 3) { print "Error: need >=3 numeric points"; exit 1 }

  # Close loop if not already closed
  if (x[0] != x[cnt-1] || y[0] != y[cnt-1]) { x[cnt]=x[0]; y[cnt]=y[0]; cnt++ }

  print "solid airfoil" > out;

  # Side walls (skip degenerate segments)
  for (i=0; i<cnt-1; i++) {
    x1=x[i]; y1=y[i]; x2=x[i+1]; y2=y[i+1];
    dx=x2-x1; dy=y2-y1; if (dx*dx+dy*dy < 1e-18) continue;
    facet(x1,y1, half,  x2,y2, half,  x2,y2,-half);
    facet(x1,y1, half,  x2,y2,-half,  x1,y1,-half);
  }

  # Caps: fan around centroid
  m = cnt-1; cx=cy=0.0;
  for (i=0;i<m;i++){ cx+=x[i]; cy+=y[i] } cx/=m; cy/=m;

  # Top cap (+z): include wrap-around triangle (last -> first)
  for (i=0;i<m-1;i++){
    facet(cx,cy, half,  x[i],y[i], half,  x[i+1],y[i+1], half);
  }
  facet(cx,cy, half,  x[m-1],y[m-1], half,  x[0],y[0], half);

  # Bottom cap (-z): reversed winding, include wrap-around
  for (i=0;i<m-1;i++){
    facet(cx,cy,-half,  x[i+1],y[i+1],-half,  x[i],y[i],-half);
  }
  facet(cx,cy,-half,  x[0],y[0],-half,  x[m-1],y[m-1],-half);

  print "endsolid airfoil" >> out; close(out);
}
' "$infile"

echo "✅ Wrote STL: $outfile"
