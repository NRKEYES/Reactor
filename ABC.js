
var flag = true ;
var a = new Array() ; // matrix to diagonalize starting at i=1
var r = new Array() ;
var n = 3 ; // order of matrix
var x = new Array() ;
var y = new Array() ;
var z = new Array() ;
var ele = new Array() ;
var m = new Array() ;
var Names = new Array("H","D","Li","B","C","N","O","F","P","S","Cl","K","Ca",
   "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Br","I") ;
var mass = new Array(1.0078,2.0140,7.0160,11.0093,12.0,14.0031,15.9949,18.9984,30.9738,
    31.9721,34.9688,38.9637,39.9626,44.9559,47.9479,50.944,51.9405,54.9380,55.9349,
    58.9332,57.9353,62.9296,63.9291,68.9256,73.9211,74.9216,78.9183,126.9045) ;
var nnames = 28 ;
// fundmental constants
var h = 6.62608E-34 ;
var h2ok = 3.18002e-44 ; // h^2/k
var clight = 2.997925e8 ;
var Tc = 298.15 ;
var RJ = 8.31451 ;
var NA = 6.02214E+23 ;
var pi = 3.141592654 ;
var kB = 1.38066E-23 ;
var Po = 1.0e5 ; // use Po = 101325 for 1 atm standard state
var lamdainv = Math.pow(2*pi*kB/1000/NA,1.5)/Math.pow(h,3)*RJ/Po/NA ;
var natm = 0 ;
var comment = "" ;
var MM = 0.0
var Inertia = new Array(3) ;
var group = "forced" ;
var sigma = 1 ; // rotational symmetry number=order of rotational subgroup
var symtol = 0.02   ; // tolerance for symmetry matching
var symrange = 25.0 ; // factor for increasing symtol for large molecules
var TK = 0 ; // temperature

function getinput(form) {
var line = new Array() ;
var instr = new Array(4) ;
var flag = true ;
var rews = /[\s\t]+/g ; // regular expression matching white space
var resws = /^\s+/ ; // matches only white space at the beginning of the line
var Ae = 0.0 ; var Be = 0.0 ; var Ce = 0.0 ; var qr = 0.0 ;
var Sr = 0.0 ; var Str = 0.0 ; var Gibbs = 0.0 ;
var lstart = 1 ; var lend = 4 ; var iline = -1 ;

MM = 0.0 ; group = "forced" ; sigma = 1 ;
line = form.XYZ.value.split("\n") ;
if ( line.length <= 1 ) line = form.XYZ.value.split("\r") ;
comment = line[0].substr(0,line[0].length-1) ; // for mol files
for ( var i=0 ; i<5 ; i++ ) {
  iline++ ; natm = parseInt(line[iline]) ;
  if ( natm>0 ) break ;
}
if ( isNaN(natm)||(natm<2) ) { alert("too few atoms") ; return }
if ( iline==0 ) {
// for XYZ files
iline++ ; comment = line[iline].substr(0,line[iline].length-1) }
iline++ ;
for (i=0 ; i<natm ; i++ ) {
instr = line[i+iline].replace(resws,"") ;
instr = instr.split(rews,4) ;
ele[i] = instr[0] ;
if ( isNaN(parseFloat(ele[i])) ) {
x[i] = parseFloat(instr[1]) ;
y[i] = parseFloat(instr[2]) ;
z[i] = parseFloat(instr[3]) }
  else {
   // mol format
   ele[i] = instr[3] ;
   x[i] = parseFloat(instr[0]) ;
   y[i] = parseFloat(instr[1]) ;
   z[i] = parseFloat(instr[2]) }

  for (j=0 ; j<nnames ; j++ ) { if ( ele[i]==Names[j] ) break ; }
  if ( j == nnames ) { alert("Unknown element "+ele[i]) ; return }
  m[i] = mass[j] ;
  MM += m[i] ;
} ; // end for i

form.Out.value = "" ;
//form.Out.value = "Input coordinates\n" ;
//for (i=0 ; i<natm ; i++ ) {
//form.Out.value += ele[i]+" "+m[i]+" "+x[i]+" "+y[i]+" "+z[i]+"\n" ;}

moment() ;

eigen(a,r,n) ;
//vectors(form) ;

//form.Out.value += "___________________________________________\n" ;
form.Out.value += " Moments of Inertia\n" ;
for ( i=1 ; i<4 ; i++ ) {
if ( i == 1 ) str = "Ix" ;
if ( i == 2 ) str = "Iy" ;
if ( i == 3 ) str = "Iz" ;
Inertia[i] = a[i+(i*i-i)/2] ;
form.Out.value += str+" = "+fmt(Inertia[i],'####.####')+" g mol-1 Ang2 " ;
ISI = trunc(Inertia[i]/NA/1000*1e-20/1e-45,5) ;
form.Out.value += " or "+ISI+"e-45 kg m2\n" ;
}

form.Out.value += "___________________________________________\n Rotational constants\n" ;
  Ae = getBe(form,Inertia[1],"T") ;
  Be = getBe(form,Inertia[2],"T") ;
  Ce = getBe(form,Inertia[3],"T") ;
   lstart = 1 ; lend = 4 ;
   if ( (Ce<0.001)||(isNaN(Ce)) ) { lstart = 2 ; lend = 3 } ; // catch linear molecules
for ( i=lstart ; i<lend ; i++ ) {
if ( i == 1 ) str = "Ae" ;
if ( i == 2 ) str = "Be" ;
if ( i == 3 ) str = "Ce" ;
form.Out.value += str+" = "+trunc(getBe(form,Inertia[i],"C"),5)+" cm-1" ;
form.Out.value += " or "+trunc(getBe(form,Inertia[i],"G"),5)+" GHz" ;
form.Out.value += " or "+trunc(getBe(form,Inertia[i],"T"),5)+" K\n" ;
}

pcs() ; // rotate to principal coordinates frame

// symmetry determination
symtol = form.symtol.value ;
sigma = form.sigma.value ;
if ( isNaN(sigma) || (sigma<1)) {
sigma = getsigma(Ae,Be,Ce,form) ; // find symmetry group and rotational order
}
form.Out.value += "Point group "+group+"\n" ;
form.Out.value += "___________________________________________\n Thermodynamics\n" ;
   TK = form.TK.value ;
qr = getqr(Ae,Be,Ce,TK,form) ;
form.Out.value += "qr="+trunc(qr,5)+" with sigma="+sigma+"\n" ;
  Sr = getSr(Ae,Be,Ce,TK,form) ;
form.Out.value += "Rotational Entropy="+trunc(Sr,5)+" J mol-1 K-1\n" ;
  Str = getStr(MM,TK,form) ;
form.Out.value += "Translational Entropy="+trunc(Str,7)+" J mol-1 K-1\n" ;
  if ( (group != "forced")||(sigma>2) ) {
   if ( (group=="Cinfv")||(group=="Dinfh") ) { Gibbs = (7*RJ/2-Sr-Str)*TK/1000 }
    else { Gibbs = (4*RJ-Sr-Str)*TK/1000 }
form.Out.value += "Trans+Rotation Gibbs Free Energy="+trunc(Gibbs,6)+" kJ mol-1\n" ;
  }


form.Out.value += "=============================================\n" ;
form.Out.value += "\n Coordinates in Principal Coordinates Frame\n" ;
form.Out.value += "___________________________________________\n" ;
var frmt = "####.####" ;
form.Out.value += natm+"\n" ;
form.Out.value += comment+" in pcf\n" ;
for (i=0 ; i<natm ; i++ ) {
  str = ( ele[i].length == 1 ) ? ele[i]+" " : ele[i] ;
form.Out.value += str+" "+fmt(x[i],frmt)+" "+fmt(y[i],frmt)+" "+fmt(z[i],frmt)+"\n" ;}

} ; // end fcn getinput

function moment() {
// calculate the moment of inertia matrix
var sxx = 0.0 ;
var syy = 0.0 ;
var szz = 0.0 ;
var sxy = 0.0 ;
var sxz = 0.0 ;
var syz = 0.0 ;
var xcm = 0.0 ;
var ycm = 0.0 ;
var zcm = 0.0 ;
var mt = 0.0 ;

for (i=0 ; i<natm ; i++ ) {
xcm += m[i]*x[i] ; ycm += m[i]*y[i] ; zcm += m[i]*z[i] ;
mt += m[i] ;
} ; // end for i
xcm = xcm/mt ; ycm = ycm/mt ; zcm = zcm/mt ;
for (i=0 ; i<natm ; i++ ) {
x[i] = x[i]-xcm ; y[i] = y[i]-ycm ; z[i] = z[i]-zcm ;
} ; // end for i

for (i=0 ; i<natm ; i++ ) {
sxx += m[i]*x[i]*x[i] ; syy += m[i]*y[i]*y[i] ; szz += m[i]*z[i]*z[i] ;
sxy += m[i]*x[i]*y[i] ; sxz += m[i]*x[i]*z[i] ; syz += m[i]*y[i]*z[i] ;
} ; // end for i

a[1] = syy+szz ; a[2] = -sxy ; a[3] = sxx+szz ; a[4] = -sxz ; a[5] = -syz ;
a[6] = sxx + syy ;
} ; // end fcn moment

function getBe(form,Inert,units) {
// calculate Be in different units
var str = "" ;
if ( Inert>0.001 ) {
  Inert = Inert/NA/1000*1e-20 ;
  var Be = h/8/pi/pi/Inert/clight/100 ; // cm-1
  if ( units == "G" ) Be = Be*clight*100*1e-9 ;
  if ( units == "T" ) Be = Be = h2ok/8/pi/pi/Inert ; // K
}
return Be
} ; // end function getBe

function getStr(MM,TK,form) {
// calculate translational entropy
qt = lamdainv*Math.pow(TK,2.5)*Math.pow(MM,1.5) ;
Str = RJ*Math.log(qt) + 2.5*RJ ;
return Str ;
} ; // end fcn getStr

function getqr(Ae,Be,Ce,TK,form) {
// calculate rotational partition function
var qr = 0.0 ;
if ( (Ce < 0.001)||(isNaN(Ce)) ) { qr = TK/Be/sigma }
  else { qr = Math.sqrt(pi*TK/Ae*TK/Be*TK/Ce)/sigma }
return qr
} ; // end fcn getqr

function getSr(Ae,Be,Ce,TK,form) {
// calculate rotational contribution to the entropy
var Sr = 0.0 ;
if ( (Ce < 0.001)||(isNaN(Ce)) ) { Sr = RJ*Math.log(TK/Be/sigma)+RJ }
  else { Sr = RJ*Math.log(Math.sqrt(pi*TK/Ae*TK/Be*TK/Ce)/sigma) + 1.5*RJ }
return Sr
} ; // end fcn getSr

function pcs() {
// rotate molecule into principal coordinate frame
// new X = R(tr)*X, X = (x(k),y(k),z(k)) for atom k
var t = new Array(3) ;
for (k=0 ; k<natm ; k++ ) {
   for ( i=1 ; i<=3 ; i++ ) {
    t[i] = r[1+(i-1)*3]*x[k] + r[2+(i-1)*3]*y[k] +r[3+(i-1)*3]*z[k] ;
   } ; // end for i
x[k] = t[1] ; y[k] = t[2] ; z[k] = t[3] ;
} ; // end for k
} ; // end fcn pcs

function getsigma(Ae,Be,Ce,form) {
// get the rotational symmetry number from the point group
var dAB = 0.0 ; var dBC = 0.0 ;
var Btol = 0.01 ; // tolerance for B constant matching *100=%
var order = 1 ;
var nC2 = 0 ;
var prinax = "x" ;
var hplane = "yz" ;
var vplane1 = "xy" ;
var vplane2 = "xz" ;
var subgroup = "" ;
sigma = 1 ; group = "forced" ;
if ( symtol>0.2 ) Btol = 0.1 ; // loosen up symmetric top for very coarse symmetry
if ( symtol>0.2 ) symrange = 12.0 ; // loosen up symtol scaling for very coarse sym

//linear 0xC2 axis Cinfv, 1xC2 axis Dinfh 
if ( (Ce<0.003)||(isNaN(Ce)) ) {
 group = "Cinfv" ;
 if ( symcmp("x",180.0) ) {sigma = 2 ; group = "Dinfh" }
} ; // end if linear

dAB = Math.abs(Ae - Be)/Ae ; dBC = Math.abs(Be - Ce)/Be ;
// asymmetric top 1xC2 axis -> C2x, 3xC2 axes -> D2x
if  ( (dAB>Btol)&&(dBC>Btol) ) {
  if ( symcmp("x",180.0) ) { nC2++ }
  if ( symcmp("y",180.0) ) { nC2++ ; hplane = "xz" ; vplane2 = "yz" }
  if ( symcmp("z",180.0) ) { nC2++ ; hplane = "xy" ; vplane2 = "yz" }
 if ( nC2>0 ) order = 2 ;
 group = ( nC2>1 ) ? "D" : "C" ;
 subgroup = ( symcmp(vplane1,0.0)||symcmp(vplane2,0.0) ) ? "v" : "" ;
 subgroup = ( symcmp("i",0.0) ) ? "i" : subgroup ;
 subgroup = ( symcmp("xyz",90.0) ) ? "d" : subgroup ;
 subgroup = ( symcmp(hplane,0.0) ) ? "h" : subgroup ;
 if ( (group=="C")&&(subgroup=="d") ) { group = "S" ; subgroup = "" }
 group += ( group=="S" ) ? 2*order : order ;
 group += subgroup ;
   if ( (group=="C1v")||(group=="C1h") ) group = "Cs" ;
   if ( group=="C1i" ) group = "Ci" ;
 sigma = ( nC2>1 ) ? 2*order : order ;
} ; // end if asymmetric top

// spherical top 0xi -> Td, 0xC5 -> Oh, 1xC5 -> Ih
if ( (dAB<=Btol)&&(dBC<=Btol) ) {
     order = aligna(form) ; alignv("x",form) ;
     if ( (order==2)||(order==3) ) { order = aligndi(order) } ; // try to find C4 or C5
     if ( order==5 ) { sigma = 60 ; group = symcmp("i",0.0) ? "Ih" : "I" }
     if ( order==4 ) { sigma = 24 ; group = symcmp("i",0.0) ? "Oh" : "O" }
     if ( order==3 ) { sigma = 12 ; 
       if ( symcmp("i",0.0) ) { group = "Th" }
        else { group = symcmp("xy",0.0) ? "Td" : "T" }
      } ; // end if order==3
} ; // end if symmetric top

// symmetric top 0xC2(+) -> Cnx, >0xC2(+) -> Dnx
if ( ((dAB<=Btol)&&(dBC>Btol))||((dAB>Btol)&&(dBC<=Btol)) ) {
 for ( nord=6 ; nord>1 ; nord-- ) { if ( symcmp("x",360.0/nord) ) break }
   order = nord ; prinax = "x" ;
 for ( nord=6 ; nord>1 ; nord-- ) { if ( symcmp("z",360.0/nord) ) break }
   if ( nord > order ) { order = nord ; prinax = "z" } ;
   if ( (order < 2)&&(dAB<Btol) ) prinax = "z" ; // for D2 A=B>C
// check for C2
  if ( alignC2(prinax,form) ) { 
          if ( order == 1 ) { order = 2 ; prinax = "y" }
           else { nC2++ } }
   if ( prinax == "y" ) { hplane = "xz" ; vplane1 = "xy" ; vplane2 = "yz" }
   if ( prinax == "z" ) { hplane = "xy" ; vplane1 = "yz" }
  if ( (prinax!="x")&&symcmp("x",180.0) ) nC2++ ;
  if ( (prinax!="y")&&symcmp("y",180.0) ) nC2++ ;
  if ( (prinax!="z")&&symcmp("z",180.0) ) nC2++ ;
 if ( (nC2>0)&&(order<2) ) { alert("oops, Symmetry tolerance too coarse") ; return }
 group = ( nC2>0 ) ? "D" : "C" ;
//  get subgroup
if ( group!="D" ) {
 if ( symcmp(vplane1,0.0)||symcmp(vplane2,0.0) ) { subgroup = "v" }
   else { subgroup = ( alignv(prinax,form) ) ? "v" : "" }
} ; // end if group != "D"
 subgroup = ( symcmp(prinax+hplane,360.0/2/order) ) ? "d" : subgroup ;
 subgroup = ( symcmp(hplane,0.0) ) ? "h" : subgroup ;
 if ( (group=="C")&&(subgroup=="d") ) { group = "S" ; subgroup = "" }
 group += ( group=="S" ) ? 2*order : order ;
 group += subgroup ;
 sigma = ( nC2>0 ) ? 2*order : order ;
} ; // end if symmetric top

return sigma
} ; // end fcn getsigma

function symcmp(axis,ang) {
// rotate point x,y,z around axis=x,y,z by ang
// reflect point x,y,z across plane=axis=xz,xy,yz
// improper rotation axis-xyz,yxz,zxy
// invert point x,y,z through cm axis=i
// new X = T*X, X = (x(k),y(k),z(k)) for atom k
// symtol = tolerance for symmetry matching
// symrange = scaling factor for increasing stol for large molecules
var plane = axis ;
var xt = 0.0 ; var yt = 0.0 ; var zt = 0.0 ;
var cosa = Math.cos(ang/360*2*pi) ;
var sina = Math.sin(ang/360*2*pi) ;
var dx = 0.0 ; var dy = 0.0 ; var dz = 0.0 ;
var d2 = 0.0 ; var stol = 0.0 ;
  if ( axis.length > 2 ) {
    plane =  axis.substring(1,3) ; axis =  axis.substring(0,1) }

// test each atom
for (var k=0 ; k<natm ; k++ ) {
    xt = x[k] ; yt = y[k] ; zt = z[k] ;
  if ( axis == "x" ) { yt = cosa*y[k] - sina*z[k] ; zt =  sina*y[k] + cosa*z[k] }
  if ( axis == "y" ) { xt = cosa*x[k] + sina*z[k] ; zt = -sina*x[k] + cosa*z[k] }
  if ( axis == "z" ) { xt = cosa*x[k] - sina*y[k] ; yt =  sina*x[k] + cosa*y[k] }
  if ( axis == "i" ) { xt = -xt ; yt = -yt ; zt = -zt }
  if ( plane == "xz" ) { yt = -yt }
  if ( plane == "xy" ) { zt = -zt }
  if ( plane == "yz" ) { xt = -xt }
 // find symmetry related atom
     // loosen tol for atoms far from cm
     d2 = xt*xt+yt*yt*zt*zt ;
     stol = symtol*(1.0+d2/symrange) ;
 for (var j=0 ; j<natm ; j++ ) {
    if ( m[j]!= m[k] ) continue ;
    dx = Math.abs(x[j]-xt) ;
    dy = Math.abs(y[j]-yt) ;
    dz = Math.abs(z[j]-zt) ;
    if ( (dx<stol)&&(dy<stol)&&(dz<stol) ) break ;
 } ; // end for j
   if ( j==natm ) return false ;
} ; // end for k
return true
} ; // end fcn symcmp

function alignC2(prinax,form) {
// Rotate molecule into unique principal coordinate frame for pc(y).
// The pcf for moments along y and z are not unique for degenerate eigenvalues.
// This routine finds atom pair centroids
// and rotates around the principal to place this C2 along the y-axis.
// new X = R*X, X = (x(k),y(k),z(k)) for atom k
var yc2 = 0.0 ; // point on proposed C2 axis
var xc2 = 0.0 ; // xcoor on C2 axis
var zc2 = 0.0 ; // zcoor on C2 axis
var r = 0.0 ;
//form.Out.value += "aligning C2 axes...\n" ;
for (var j=0 ; j<natm ; j++ ) {
 for (var k=j+1 ; k<natm ; k++ ) {
  if ( m[j]!= m[k] ) continue ;

if ( prinax == "x" ) {
   if ( Math.abs(x[k]+x[j])>symtol ) continue ; // either side of horizontal mirror
  yc2 = (y[k]+y[j])/2 ;
  zc2 = (z[k]+z[j])/2 ;
   if ( (Math.abs(yc2)<symtol)&&(Math.abs(zc2<symtol)) ) continue ; // related by i
  r = Math.sqrt(yc2*yc2+zc2*zc2) ; cosa = yc2/r ; sina = zc2/r ;
} ; // end if x
if ( prinax == "z" ) {
  if ( Math.abs(z[k]+z[j])>symtol ) continue ; // either side of horizontal mirror
  yc2 = (y[k]+y[j])/2 ;
  xc2 = (x[k]+x[j])/2 ;
   if ( (Math.abs(yc2)<symtol)&&(Math.abs(xc2<symtol)) ) continue ; // related by i
  r = Math.sqrt(yc2*yc2+xc2*xc2) ; cosa = yc2/r ; sina = xc2/r ;
} ; // end if z
  rot(prinax,cosa,sina) ;
   if ( symcmp("y",180.0) ) return true ;
 } ; // end for k
} ; // end for j

return false ;
} ; // end fcn alignc2

function alignv(rotax,form) {
// Rotate molecule into unique principal coordinate frame for pc(2).
// The pcf is not unique for degenerate eigenvalues.
// This routine rotates around axis=rotax= x or z to place an atom that sits on
// a vertical mirror plane on the y-axis.
// new X = R*X, X = (x(k),y(k),z(k)) for atom k
var r = 0.0 ;
var cosa = 0.0 ; var sina = 0.0 ;
var vplane1 = ( rotax == "x" ) ? "xy" : "yz" ;
//form.Out.value += "aligning vertical mirror...\n" ;
for (var i=0 ; i<natm ; i++ ) {
 // dont' use atoms that are on principal axis
 if ( (rotax == "x")&&(Math.abs(y[i])<symtol)&&(Math.abs(z[i])<symtol) ) continue ;
 if ( (rotax == "z")&&(Math.abs(y[i])<symtol)&&(x[i]<symtol) ) continue ;
 if ( rotax=="x" ) { r = Math.sqrt(y[i]*y[i]+z[i]*z[i]) ; cosa = y[i]/r ; sina = z[i]/r }
 if ( rotax=="z" ) { r = Math.sqrt(y[i]*y[i]+x[i]*x[i]) ; cosa = y[i]/r ; sina = x[i]/r }
  rot(rotax,cosa,sina) ;
 if ( symcmp(vplane1,0.0) ) { return true }
} ; // end for i
return false ;
} ; // end fcn alignv

function aligna(form) {
// Rotate molecule into unique principal coordinate frame for pc(x).
// The pcf is not unique for degenerate eigenvalues.
// This routine is specific for cubic groups.
// This routine rotates atoms to the x axis and finds maximum order,
// failing that rotates atom pair centroids to the x axis and fins maximum order.
// order=2,3(Th) ; 2,3,4(Oh) ; 2,3,5(Ih) but the 5-fold will not be on an atom
// new X = R*X, X = (x(k),y(k),z(k)) for atom k
var mxorder = 1 ;
var mxi = 0 ; var mxj = 0 ;

// rotate atoms to the x axis
for (var ia=0 ; ia<natm ; ia++ ) {
 if ( (Math.abs(y[ia])>symtol)||(Math.abs(z[ia])>symtol) ) {
  alignx(x[ia],y[ia],z[ia]) ;
 for ( nord=5 ; nord>1 ; nord-- ) { if ( symcmp("x",360.0/nord) ) break }
   if ( nord > mxorder ) { mxorder = nord ; mxi = ia }
   if ( mxorder>3 ) break ; // only Oh has a C4
 } ; // end if y!=0 z!=0
} ; // end for ia
alignx(x[mxi],y[mxi],z[mxi]) ;

if ( mxorder<4 ) {
// rotate pair centroids to x axis
//form.Out.value += "aligning centroids ...\n" ;
for (var ia=0 ; ia<natm ; ia++ ) {
 for (var ja=ia+1 ; ja<natm ; ja++ ) {
  if ( m[ia]!= m[ja] ) continue ;
  xc2 = (x[ia]+x[ja])/2 ;
  yc2 = (y[ia]+y[ja])/2 ;
  zc2 = (z[ia]+z[ja])/2 ;
   if ( (Math.abs(yc2)<symtol)&&(Math.abs(zc2<symtol)) ) continue ; // already tested
  alignx(xc2,yc2,zc2) ;
 for ( nord=5 ; nord>1 ; nord-- ) { if ( symcmp("x",360.0/nord) ) break }
   if ( nord >= mxorder ) { mxorder = nord ; mxi = ia ; mxj = ja }
   if ( mxorder>3 ) break ; // only Oh has a C4
 } ; // end for ja
 if ( mxorder>3) break ;
} ; // end for ia
alignx((x[mxi]+x[mxj])/2,(y[mxi]+y[mxj])/2,(z[mxi]+z[mxj])/2) ;
} ; // end if mxorder<4

return mxorder ;
} ; // end fcn aligna

function alignx(xp,yp,zp) {
// This routine rotates point x,y,z to the x axis.
var r2 = 0.0 ;
var cosa = 0.0 ; var sina = 0.0 ;
 if ( (Math.abs(yp)<symtol)&&(Math.abs(zp)<symtol) ) return ;
  r2 = yp*yp+zp*zp ;
  cosa = yp/Math.sqrt(r2) ;
  sina = zp/Math.sqrt(r2) ;
  rot("x",cosa,sina) ;
  cosa = xp/Math.sqrt(xp*xp+r2) ;
  sina = -Math.sqrt((r2)/(xp*xp+r2)) ;
  rot("z",cosa,sina) ;
} ; // end fcn alignx

function aligndi(order) {
// Rotate molecule into unique principal coordinate frame for pc(x).
// The pcf is not unique for degenerate eigenvalues.
// THis routine rotates the dihedral plane perpendicular to the calculated C5
// axis for I and Ih to the x axis and lloks for a C5. If this fails,
// this routine rotates the dihedral plane perpendicular to (1,1,1) to the x axis
// and finds maximum order. This routine is really only necessary if atom
// alignment produces an order of 2 or 3 for spherical tops. Use aligna and alignv
// before calling this routine.
// order=3(Td) ; 2,3(Th) ; 3,4(Oh) ; 3(Ih) after aligna()
// new X = R*X, X = (x(k),y(k),z(k)) for atom k
var mxorder = 1 ;
var cosa = 0.0 ; var sina = 0.0 ;
// aligned along trigonal axis rotate to C5 axis
if ( order==3 ) {
nord = 5 ;
  cosa = Math.cos(2*pi*37.38/360) ;
  sina = Math.sin(2*pi*37.38/360) ;
  rot("z",cosa,-sina) ;
if ( symcmp("x",360.0/nord) ) return nord ;
  rot("z",cosa,sina) ;
  rot("z",cosa,sina) ;
if ( symcmp("x",360.0/nord) ) return nord ;
  rot("z",cosa,-sina) ; // return to original
}
// aligned along trigonal axis rotate to C4 axis
// aligned along C4 axis for Oh or C2 for Th
if ( (order==4)||(order==2) ) {
  cosa = 1/Math.sqrt(2) ;
  rot("x",cosa,-cosa) ;
 }
  cosa = 1/Math.sqrt(3) ;
  sina = Math.sqrt(2)*cosa ;
  rot("z",cosa,-sina) ;
 for ( nord=4 ; nord>1 ; nord-- ) { if ( symcmp("x",360.0/nord) ) break }
  if ( order>nord ) {
     // go back to original alignment
     rot("z",cosa,sina) ;
     if ( (order==4)||(order==2) ) {
       cosa = 1/Math.sqrt(2) ;
       rot("x",cosa,cosa) ;
     }
   nord = order ;
   } ; // end if order>nord
return nord ;
} ; // end fcn aligndi


function rot(rotax,cosa,sina) {
// rotate molecule around axis=rotax= x or z
var yt = 0.0 ;
if ( rotax == "x" ) {
  for (var k=0 ; k<natm ; k++ ) {
    yt = cosa*y[k] + sina*z[k] ; z[k] = -sina*y[k] + cosa*z[k] ;
    y[k] = yt ;
  } ; // end for k
} ; // end if x

if ( rotax == "z" ) {
  for (k=0 ; k<natm ; k++ ) {
    yt = sina*x[k] + cosa*y[k] ; x[k] = cosa*x[k] - sina*y[k] ;
    y[k] = yt ;
  } ; // end for k
} ; // end if z
} ; // end fcn rot


function eigen(a,r,n) {
// Diagonalize a real symmetric matrix
// calculate eigenvalues and vectors ;
// The matrix should be in a 1-D array as the lower diagonal portion
// of the matrix, by rows. Coverted from the classic fortran routine into
// JavaScript by Tom Shattuck, Dept. of Chemistry, COlby College
// twshattu@colby.edu. Please feel free to grab this application,
// but please leave the acknowledgements.
// n = dimension of the square symmetric matrix
var x = 0 ; var y = 0 ; var mm = 0 ; var mq = 0 ; var lq = 0 ; var lm = 0 ;
var l = 0 ; var m = 0 ; var anorm ; var anrmx ; var range ;
range=1.0e-12 ;
iq=-n ;
   for (var j=1 ; j<=n ; j++ ) {
	iq += n ;
   for (var i=1 ; i<=n ; i++ ) {
	ij=iq+i ;
	r[ij]=0.0 ;
	if (i == j) r[ij]=1.0 ;
		} ; // end for i
	} ; // end for j
anorm=0.0 ;
   for ( i=1 ; i<=n ; i++ ) {
      for (j=1 ; j<=n ; j++ ) {
	if (i != j) {
		ia=i+(j*j-j)/2 ;
		anorm=anorm+a[ia]*a[ia] ;
				} ; // end if not diagonal element
			} ; // end for j
		} ; // end for i
if ( anorm > 0 ) {
   anorm=1.414* Math.sqrt(anorm) ;
   anrmx=anorm*range/n ;

// initialize indicators and compute threshold,thr ;
redo = false ;
thr=anorm ;
// compare threshold with final norm ;
while ( thr > anrmx ) {
   thr=thr/n ;
  for ( pass=0 ; pass<1000 ; pass++ ) {
   l=1 ;
//  test for l=second from last column ;
   while ( l <= (n-1) ) {
	m=l+1 ;
//   test for m=last column
     while ( m <= n ) {
	mq=(m*m-m)/2 ;
	lq=(l*l-l)/2 ;
	lm=l+mq ;
         if ( Math.abs(a[lm]) >= thr ) {
	redo = true ;
	ll=l+lq ;
	mm=m+mq ;
	x=0.5*(a[ll]-a[mm]) ;
	y=-a[lm]/ Math.sqrt(a[lm]*a[lm]+x*x) ;
	if ( x < 0 ) y=-y ;
	sinx=y/Math.sqrt(2.0*(1.0+( Math.sqrt(1.0-y*y)))) ;
	sinx2=sinx*sinx ;
	cosx=Math.sqrt(1.0-sinx2) ;
	cosx2=cosx*cosx ;
	sincs=sinx*cosx ;
//	rotate l and m columns ;
	ilq=n*(l-1) ;
	imq=n*(m-1) ;
	for ( i=1 ; i <=n ; i++ ) {
	     iq=(i*i-i)/2 ;
	     if ( i != l ) {
		if ( i != m ) {
		     if ( i < m ) { im=i+mq } else {im=m+iq }
		     if ( i < l ) { il=i+lq } else { il=l+iq }
			 x=a[il]*cosx-a[im]*sinx ;
			a[im]=a[il]*sinx+a[im]*cosx ;
			a[il]=x ;
		} ; // end if i != m
	} ; // end if i != l
		ilr=ilq+i ;
		imr=imq+i ;
		x=r[ilr]*cosx-r[imr]*sinx ;
		r[imr]=r[ilr]*sinx+r[imr]*cosx ;
		r[ilr]=x ;
	} ; // end for i
	x=2.0*a[lm]*sincs ;
	y=a[ll]*cosx2+a[mm]*sinx2-x ;
	x=a[ll]*sinx2+a[mm]*cosx2+x ;
	a[lm]=(a[ll]-a[mm])*sincs+a[lm]*(cosx2-sinx2) ;
	a[ll]=y ;
	a[mm]=x ;
             } ; // end if a(lm) >= thr

// test for completion
	   m++ ;
         } ; // end while m<= n
	l++ ;
      } ; // end while l <= n-1
		if ( ! redo ) break ;
		redo = false ;
      } ; // end for pass
    } ; // end while thr>anrmx
  } ; // end if anorm>0
// sort eigen values and vectors ;
iq=-n ;
for ( i=1 ; i<=n ; i++ ) {
   iq=iq+n ;
   ll=i+(i*i-i)/2 ;
   jq=n*(i-2) ;
for ( j=i ; j<=n ; j++ ) {
	jq=jq+n ;
	mm=j+(j*j-j)/2 ;
	if( a[ll] < a[mm] ) {
	  x=a[ll] ;
	  a[ll]=a[mm] ;
	  a[mm]=x ;
	      for ( k=1; k<=n ; k++) {
		ilr=iq+k ;
		imr=jq+k ;
		x=r[ilr] ;
		r[ilr]=r[imr] ;
		r[imr]=x ;
	         } ; // end for k
	   } ; // end if a[ll] < a[mm]
	} ; // end for j
    } ; // end for i
} ; // end eigen

// print out eigenvectors--------------------
function vectors(form) {

form.Out.value = "" ;
for ( i=1 ; i<=n ; i++ ) {
form.Out.value += "Eigenvector "+i+": E="+trunc(a[i+(i*i-i)/2],6)+"\n" ;
 for ( j=1 ; j<=n ; j++ ) {
	form.Out.value += ""+trunc(r[j+(i-1)*n],6)+"\n" 
	}
form.Out.value += "-------------\n" 
   }

} ; // end vectors

function examples() {
// fetch an example molecule from exampleABC.html
Example = window.open('exampleABC.html','example',
 'TOOLBAR=no,MENUBAR=no,SCROLLBARS=yes,RESIZABLE=yes,STATUS=yes,WIDTH=700,HEIGHT=400') ;
} ; // end examples

function fmt(xfx,format) {
// typical call example fmt(x,'#####.####')
// note special rounding for mass spec assuming resolution about 0.2 amu
var spc = '          ' ;
var idecfmt = format.indexOf('.',0) ;
var iendfmt = format.length-idecfmt-1 ;
var ipstdecpt = format.length-idecfmt-1 ;
 xfx += 0.19*Math.pow(10,-ipstdecpt) ; // rounding
var xstr = ''+xfx ; // casting
 idec = xstr.indexOf('.',0) ;
   if ( idec < 0 ) {
       xstr = xstr + '.' ;
       idec = xstr.length-1 ;
    }
 iend = xstr.length-idec-1 ;
   trnc = xstr.length-iend+iendfmt ;
   zeros = idecfmt-iend-1 ;
 xstr = (iend>iendfmt) ? xstr.substring(0,trnc) : xstr+spc.substring(0,zeros) ;
 str = spc.substring(0,idecfmt-idec)+xstr ;
 return str ;
} ; // end fmt

// Significant figure functions
function ord(x) {
   return Math.floor(Math.log(Math.abs(x+1e-35))/2.303)
}

// Truncate to n sign. figures
function trunc(x,n) {
var c= Math.floor(x*Math.pow(10,-ord(x)+n-1)+.5)/Math.pow(10,-ord(x)+n-1)
   c = ( Math.abs(c)<1e-13) ? 0 : c ;
   return c
}

