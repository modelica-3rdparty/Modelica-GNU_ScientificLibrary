within GNU_ScientificLibrary.Blocks.Physics;
block DebyeEoS "Debye Model EoS"
  // block by tommy.burch@physik.uni-r.de
  extends GNU_ScientificLibrary.Icons.Block;
  import const = Modelica.Constants;
  parameter Real T_D = 343 "K ; Debye Temperature";
  parameter Real eff_am = 63.546 "Da ; Effective atomic mass";
  parameter Real mdens = 8935 "kg/m^3 ; Mass density";
  final parameter Real ndens = mdens/(eff_am*1.66053906892e-27) "1/m^3 ; Number density of atoms";
  final parameter Real c_s = (2 * T_D * const.k / const.h) * (const.pi/(6*ndens))^0.333333333333 "m/s ; Sound speed";
  Real x,Bx,D3x;
  Modelica.Blocks.Interfaces.RealInput T annotation(
    Placement(transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealOutput edens annotation(
    Placement(transformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput c_v annotation(
    Placement(transformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}})));
equation
  x = T_D / T;
  Bx = x / (exp(x) - 1);
  D3x = GNU_ScientificLibrary.Functions.specfunc.debye_D3(x);
  // energy density:
  edens = 3*ndens*const.k*T*D3x;
  // specific heat capacity:
  //  using...  x D'_n(x) = n ( B(x) - D_n(x) )
  // c_v = ( edens/T - 9*ndens*const.k*(Bx - D3x) ) / mdens;
  c_v = ndens*const.k*(12*D3x - 9*Bx) / mdens;
  annotation (Icon(graphics = {Text( origin = {4, -8}, extent = {{-84, 48}, {76, -33}}, textString = "Debye
EoS")}),
  Documentation(info = "<html><head></head><body>Debye-model Equation of State (EoS), appropriate for the phonon contribution in solids.<div><br></div><div>Energy density:</div><div><br></div><div>e = 3 n k T D<sub>3</sub>(x)</div><div><br></div><div>where n is the number density of atoms,&nbsp;D<sub>3</sub>(x) is the third Debye function, x=T<sub>D</sub>/T, and T<sub>D</sub> is the Debye temperature:</div><div><br></div><div>T<sub>D</sub> = [h c<sub>s</sub> / (2 k)] (6 n / π)<sup>1/3</sup>&nbsp;</div><div><br></div><div>c<sub>s</sub> being the sound speed in the material (value obtainable via the parameter, *.c_s).</div><div><br></div><div>With proper values of T<sub>D</sub>&nbsp;(in K), the atomic mass (in Da), and the mass density (as parameters), this block will provide the energy density (in J m<sup>-3</sup>) as a function of the input temperature (in K): e(T). The specific heat capacity, c<sub>v</sub>(T) (in J kg<sup>-1</sup> K<sup>-1</sup>), is also calculated from the derivative wrt T:</div><div><br></div><div>c<sub>v</sub> = [e/T - 9 n k (B(x)-D<sub>3</sub>(x))] / ρ = n k [12 D<sub>3</sub>(x) - 9 B(x)] / ρ&nbsp;</div><div><font face=\"DejaVu Sans Mono\"><br></font></div><div><font face=\"DejaVu Sans Mono\">where ρ is the mass density, B(x)=x/(exp(x)-1), and we have used xD'<sub>3</sub>(x)=3[B(x)-D<sub>3</sub>(x)] .</font></div><div><font face=\"DejaVu Sans Mono\"><br></font></div><div><font face=\"DejaVu Sans Mono\">The default parameter values are those appropriate for copper.</font></div></body></html>"));
end DebyeEoS;