within GNU_ScientificLibrary.Blocks.Physics;

block FermiDiracEoS "Fermi-Dirac EoS"
  // block by tommy.burch@physik.uni-r.de
  extends GNU_ScientificLibrary.Icons.Block;
  import const = Modelica.Constants;
  final constant Real m_e = 9.10938e-31 "mass of the electron in kg";
  final constant Real k_dens = 2 * (2*const.pi*m_e*const.k/const.h^2)^1.5;
  final constant Real k_pres = k_dens * const.k;
  parameter Real eff_am = 63.546 "Da ; Effective atomic mass";
  parameter Integer N_e = 1 "# conduction e-'s per atom";
  parameter Real mdens = 8935 "kg/m^3 ; Mass density";
  final parameter Real ndens = N_e*mdens/(eff_am*1.66053906892e-27) "1/m^3 ; Number density of e-";
  final parameter Real E_Fermi = 0.5*(0.5*const.h/const.pi)^2 * (3*const.pi^2)^0.666666666667 * ndens^0.666666666667 / m_e "J ; Fermi Energy";
  final parameter Real edens0 = 0.6*ndens*E_Fermi "J/m^3 ; Energy density at T=0";
  Real xmin, xmax, x, y, F3half, edens, mu, Fmhalf, dxdT;
  //
  function FhalfMinusF
    extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
    input Real F;
  algorithm
    y := GNU_ScientificLibrary.Functions.specfunc.fermidirac_Fhalf(u) - F;
  end FhalfMinusF;
  //
  Modelica.Blocks.Interfaces.RealInput T annotation(
    Placement(transformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealOutput e_minus_e0 annotation(
    Placement(transformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput c_v annotation(
    Placement(transformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput mu_minus_Ef annotation(
    Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
equation
  y = ndens * T^(-1.5) / k_dens;
if N_e <= 0 then
  x = 0;
  xmin = 0;
  xmax = 0;
  F3half = 0;
  Fmhalf = 1;
else
  if y < 1.0 then
// Maxwell-Boltzmann asymp. expansion
    xmin = log(y);
    xmax = xmin + 0.5;
  elseif y > 3.0 then
// degenerate asymp. expansion
    xmax = (0.75*y*sqrt(const.pi))^0.6666666666667;
    xmin = xmax - 0.5;
  else
    xmin = log(y);
    xmax = (0.75*y*sqrt(const.pi))^0.6666666666667;
  end if;
  x = Modelica.Math.Nonlinear.solveOneNonlinearEquation(function FhalfMinusF(F = y), xmin, xmax);
  F3half = GNU_ScientificLibrary.Functions.specfunc.fermidirac_F3half(x);
  Fmhalf = GNU_ScientificLibrary.Functions.specfunc.fermidirac_Fmhalf(x);
end if;
  mu = x * const.k * T;
  mu_minus_Ef = mu - E_Fermi;
  edens = 1.5 * k_pres * T^2.5 * F3half;
  e_minus_e0 = edens - edens0;
  dxdT = -1.5 * ndens / (k_dens * T^2.5 * Fmhalf);
  c_v = ( 2.5 * edens / T + 1.5 * k_pres * T^2.5 * y * dxdT ) / mdens;
  annotation(
    Icon(graphics = {Text(origin = {4, -8}, extent = {{-84, 48}, {76, -33}}, textString = "FD
EoS")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
    Documentation(info = "<html><head></head><body><div>Fermi-Dirac Equation of State (EoS), appropriate for a gas of fermions, like the conduction electrons in a metal or plasma.</div><div><br></div><div>Parameters for the atomic mass, # of conduction e-'s per atom, and mass density (ρ) set the number density (n) of e-. Together with the temperature, T (in K), the relationship with the chemical potential (x = μ / kT) can be established:</div><div><br></div><div>n / T<sup>3/2</sup> ~ F<sub>1/2</sub>(x) = y</div><div><br></div>where F<sub>j</sub>(x) are the complete Fermi-Dirac integrals.<div>Numerical inversion of F<sub>1/2</sub>(x)=y: &nbsp;x = F<sub>1/2</sub><sup>-1</sup>(y). Then we can determine the energy density:<div><br></div><div>e ~ T<sup>5/2</sup> F<sub>3/2</sub>(x)</div><div><br></div><div>which also gives us the pressure, P = 2e/3. The specific heat capacity is calculated via&nbsp;</div><div><br></div><div>c_v = de/dT ~ e [2.5 / T&nbsp;+ F<sub>1/2</sub>(x) dx/dT /&nbsp;F<sub>3/2</sub>(x)] / ρ&nbsp;</div><div><br></div><div>where we have used F'<sub>3/2</sub>(x)=F<sub>1/2</sub>(x) and the derivative of the chemical potential (dx/dT) is calculated via&nbsp;</div><div><br></div><div>dx/dT ~ -1.5 n / (T<sup>5/2</sup> F<sub>-1/2</sub>(x))&nbsp;</div><div><br></div><div>thanks to the first relation above and F'<sub>1/2</sub>(x)=F<sub>-1/2</sub>(x). The chemical potential, μ = x k T, is also provided as an output.</div><div><br></div><div>The energy density (in J m<sup>-3</sup>) and chemical potential (in J) outputs are shifted relative to the values at absolute zero:&nbsp;</div><div><br></div><div>e - e<sub>0</sub> &nbsp;, where e<sub>0</sub> = 3 n E<sub>F</sub> / 5&nbsp;</div><div><br></div><div>μ - E<sub>F</sub> &nbsp;, where E<sub>F</sub> = h<sup>2</sup> (3 π<sup>2</sup> n)<sup>2/3</sup> / (8 π<sup>2</sup> m<sub>e</sub>)&nbsp;</div><div><br></div><div>E<sub>F</sub> and e<sub>0</sub> are available as parameters (*.E_Fermi and *.edens0).&nbsp;</div><div><br></div><div><span style=\"font-family: 'DejaVu Sans Mono';\">The default parameter values are those appropriate for copper.</span></div></div></body></html>"));
end FermiDiracEoS;