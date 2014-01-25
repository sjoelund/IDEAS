within IDEAS.Buildings.Components.Shading;
model SideFins "Vertical side fins next to windows"
  extends IDEAS.Buildings.Components.Interfaces.StateShading(controled=false);

  // Window properties
  parameter Modelica.SIunits.Length hWin "Window height";
  parameter Modelica.SIunits.Length wWin "Window width";
  final parameter Modelica.SIunits.Area aWin = hWin*wWin "Window area";

  // Sidefin properties
  parameter Modelica.SIunits.Length hFin "Height of side fin above window";
  parameter Modelica.SIunits.Length dep
    "Overhang depth perpendicular to the wall plane";
  parameter Modelica.SIunits.Length gap
    "Distance between window upper edge and overhang lower edge";

  Real fraSun(final min=0,final max=1, final unit="1")
    "Fraction of window area exposed to the sun";

protected
  final parameter Modelica.SIunits.Area AWin= hWin*wWin "Window area";
  final parameter Modelica.SIunits.Length tmpH[4] = {hFin+hWin,hFin,hFin+hWin,hFin}
    "Height rectangular sections used for superposition";
  final parameter Modelica.SIunits.Length tmpW[4] = {gap+wWin,gap+wWin,gap,gap}
    "Width rectangular sections used for superposition";

  Modelica.SIunits.Length x1[4]
    "Horizontal distance between side fin and point where shadow line and window lower edge intersects";
  Modelica.SIunits.Length x2
    "Horizontal distance between side fin and shadow corner";
  Modelica.SIunits.Length x3[4] "Window width";
  Modelica.SIunits.Length y1[4] "Window height";
  Modelica.SIunits.Length y2
    "Vertical distance between window upper edge and shadow corner";
  Modelica.SIunits.Length y3[4]
    "Vertical distance between window upper edge and point where shadow line and window side edge intersects";
  Modelica.SIunits.Area area[4]
    "Shaded areas of the sections used in superposition";
  Modelica.SIunits.Area shdArea "Shaded area";
  Modelica.SIunits.Area crShdArea "Final value of shaded area";
  Modelica.SIunits.Area crShdArea1
    "Shaded area, corrected for the sun behind the surface/wall";
  Modelica.SIunits.Area crShdArea2
    "Shaded area, corrected for the sun below horizon";
  Modelica.SIunits.Length minX[4];
  Modelica.SIunits.Length minY[4];
  Modelica.SIunits.Length minX2X3[4];
  Modelica.SIunits.Length minY2Y3[4];
  Real deltaL=1e-6 "Small number to avoid division by zero";
  Modelica.SIunits.Angle alt = (Modelica.Constants.pi/2) - angZen;

  Real verAzi;
  Real lambda;

equation
  lambda = tan(alt) / cos(verAzi);
  verAzi = Modelica.Math.acos(cos(angInc)/cos(alt));
  y2*Modelica.Math.cos(verAzi) = dep*Modelica.Math.tan(alt);
  x2 = dep*Modelica.Math.tan(verAzi);

  for i in 1:4 loop
    x1[i] = tmpH[i]/lambda;
    x3[i] = tmpW[i];
    y1[i] = tmpH[i];
    y3[i] = tmpW[i]*lambda;
    minX2X3[i] = IDEAS.BaseClasses.Math.MinSmooth(u1=x2,u2=x3[i],delta=deltaL);
    minX[i] = IDEAS.BaseClasses.Math.MinSmooth(u1=x1[i],u2=minX2X3[i],delta=deltaL);
    minY2Y3[i] = IDEAS.BaseClasses.Math.MinSmooth(u1=y2,u2=y3[i],delta=deltaL);
    minY[i] = IDEAS.BaseClasses.Math.MinSmooth(u1=y1[i],u2=minY2Y3[i],delta=deltaL);
    area[i] = tmpH[i]*minX[i] - minX[i]*minY[i]/2;
  end for;
  shdArea = area[4] - area[3] - area[2] + area[1];
  // correction case: Sun not in front of the wall
  crShdArea1 = Modelica.Media.Air.MoistAir.Utilities.spliceFunction(pos=shdArea,neg=AWin,x=(Modelica.Constants.pi/2)-verAzi,deltax=0.01);
  // correction case: Sun not above horizon
  crShdArea2 = Modelica.Media.Air.MoistAir.Utilities.spliceFunction(pos=shdArea,neg=AWin,x=alt,deltax=0.01);
  crShdArea=IDEAS.BaseClasses.Math.MaxSmooth(u1=crShdArea1,u2=crShdArea2,delta=0.01);
  fraSun = IDEAS.BaseClasses.Math.MinSmooth( u1=IDEAS.BaseClasses.Math.MaxSmooth(u1=1-crShdArea/AWin,u2=0,delta=0.01),u2=1.0,delta=0.01);

  iSolDir = solDir * fraSun;

  connect(solDif, iSolDif) annotation (Line(
      points={{-60,10},{40,10}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(angInc, iAngInc) annotation (Line(
      points={{-60,-50},{-14,-50},{-14,-70},{40,-70}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics), Documentation(info="<html>
<p><h4><font color=\"#008000\">General description</font></h4></p>
<p><h5>Goal</h5></p>
<p>The <code>Overhang.mo</code> model describes the transient behaviour of solar irradiance on a window below a non-fixed horizontal or vertical overhang.</p>
</html>"));
end SideFins;
