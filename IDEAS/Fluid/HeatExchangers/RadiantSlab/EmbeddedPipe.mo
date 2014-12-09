within IDEAS.Fluid.HeatExchangers.RadiantSlab;
model EmbeddedPipe
  "Embedded pipe model based on prEN 15377 and (Koschenz, 2000), water capacity lumped to TOut"
  extends IDEAS.Fluid.HeatExchangers.Interfaces.EmissionTwoPort;
  replaceable parameter
    IDEAS.Fluid.HeatExchangers.RadiantSlab.BaseClasses.RadiantSlabChar RadSlaCha constrainedby
    IDEAS.Fluid.HeatExchangers.RadiantSlab.BaseClasses.RadiantSlabChar
    "Properties of the floor heating or TABS, if present"
    annotation (choicesAllMatching=true);
  extends IDEAS.Fluid.Interfaces.Partials.PartialTwoPort(
    final m=A_pipe*L_r*rho_default,
    vol(nPorts=2),
    final allowFlowReversal = false);
  extends IDEAS.Fluid.Interfaces.TwoPortFlowResistanceParameters(
    computeFlowResistance=false,
    final dp_nominal=Modelica.Fluid.Pipes.BaseClasses.WallFriction.Detailed.pressureLoss_m_flow(
      m_flow=m_flow_nominal/nParCir,
      rho_a=rho_default,
      rho_b=rho_default,
      mu_a=mu_default,
      mu_b=mu_default,
      length=pipeEqLen/nParCir,
      diameter=pipeDiaInt,
      roughness=roughness,
      m_flow_small=m_flow_small/nParCir));

  // General model parameters ////////////////////////////////////////////////////////////////
  parameter Modelica.SIunits.Length roughness(min=0) = 2.5e-5
    "Absolute roughness of pipe, with a default for a smooth steel pipe"
    annotation(Dialog(tab="Pressure drop"));
  parameter Modelica.SIunits.Length L_floor = A_floor^(1/2)
    "Floor length - along the pipe direction"
    annotation(Dialog(tab="Pressure drop"));
  parameter Real N_pipes = A_floor/L_floor/RadSlaCha.T - 1
    "Number of parallel pipes in the slab"
annotation(Dialog(tab="Pressure drop"));
  parameter Modelica.SIunits.Length pipeBendEqLen = 2*(N_pipes-1)*(2.48*RadSlaCha.T/2/pipeDiaInt+3.20)*pipeDiaInt
    "Pipe bends equivalent length, default according to Fox and McDonald"
annotation(Dialog(tab="Pressure drop"));
  parameter Modelica.SIunits.Length pipeEqLen = pipeBendEqLen + (L_floor-2*RadSlaCha.T)*N_pipes
    "Total pipe equivalent length, default assuming 180 dg turns starting at RadSlaCha.T from the end of the slab"
annotation(Dialog(tab="Pressure drop"));
  parameter Modelica.SIunits.MassFlowRate m_flowMin = m_flow_nominal*0.5
    "Minimal flowrate when in operation - used for determining required series discretisation";
  parameter Real nParCir = 1
    "Number of parallel (equally sized) circuits in the tabs";
  parameter Modelica.SIunits.Area A_floor "Floor/tabs surface area";

  //default 15 kg/hm2 -> simplified eqn
  parameter Boolean useSimplifiedRt = m_flowMin/A_floor > 0.000416
    "Use a simplified calculation for Rt"
    annotation(Evaluate=true);
  // parameter 1/U1 + 1/U2 from Koschenz
  parameter Modelica.SIunits.ThermalInsulance R_c = 1/(RadSlaCha.lambda_b/RadSlaCha.S_1 + RadSlaCha.lambda_b/RadSlaCha.S_2)
    "Specific thermal resistivity of (parallel) slabs connected to top and bottom of tabs"
    annotation(Dialog(tab="Advanced", group="Thermal"));

  // Resistances ////////////////////////////////////////////////////////////////

  //For high flow rates see [Koshenz, 2000] eqn 4.37
  // for laminar flow Nu_D = 4 is assumed: between constant heat flow and constant wall temperature
  Modelica.SIunits.ThermalInsulance R_w_val= Modelica.Fluid.Utilities.regStep(x=rey-(reyHi+reyLo)/2,
    y1=RadSlaCha.T^0.13/8/Modelica.Constants.pi*abs((pipeDiaInt/(m_flowSp*L_r)))^0.87,
    y2=RadSlaCha.T/(4*Medium.thermalConductivity(sta_default)*Modelica.Constants.pi),
    x_small=(reyHi-reyLo)/2)
    "Flow dependent resistance value of convective heat transfer inside pipe for both turbulent and laminar heat transfer.";
  final parameter Modelica.SIunits.ThermalInsulance R_w_val_min = Modelica.Fluid.Utilities.regStep(x=m_flow_nominal/nParCir/A_pipe*pipeDiaInt/mu_default-(reyHi+reyLo)/2,
    y1=RadSlaCha.T^0.13/8/Modelica.Constants.pi*abs((pipeDiaInt/(m_flow_nominal/A_floor*L_r)))^0.87,
    y2=RadSlaCha.T/(4*Medium.thermalConductivity(sta_default)*Modelica.Constants.pi),
    x_small=(reyHi-reyLo)/2)
    "Lowest value for R_w that is expected for the set mass flow rate";
  final parameter Modelica.SIunits.ThermalInsulance R_r_val=RadSlaCha.T*log(RadSlaCha.d_a
      /pipeDiaInt)/(2*Modelica.Constants.pi*RadSlaCha.lambda_r)
    "Fix resistance value of thermal conduction through pipe wall * surface of floor between 2 pipes (see RadSlaCha documentation)";
  //Calculation of the resistance from the outer pipe wall to the center of the tabs / floorheating.
  final parameter Modelica.SIunits.ThermalInsulance R_x_val=RadSlaCha.T*(log(RadSlaCha.T
      /(3.14*RadSlaCha.d_a)) + corr)/(2*3.14*RadSlaCha.lambda_b)
    "Fix resistance value of thermal conduction from pipe wall to layer";
  final parameter Real corr = if RadSlaCha.tabs then 0 else
    sum( -(RadSlaCha.alp2/RadSlaCha.lambda_b * RadSlaCha.T - 2*3.14*s)/(RadSlaCha.alp2/RadSlaCha.lambda_b * RadSlaCha.T + 2*3.14*s)*exp(-4*3.14*s/RadSlaCha.T*RadSlaCha.S_2)/s for s in 1:10) "correction factor for the floor heating according to Multizone Building modeling with Type56 and TRNBuild (see documentation). 
    If tabs is used, corr=0";
  final parameter Modelica.SIunits.Length pipeDiaInt = RadSlaCha.d_a - 2*RadSlaCha.s_r
    "Pipe internal diameter";
  parameter Boolean from_dp = false
    "= true, use m_flow = f(dp) else dp = f(m_flow)"
    annotation (Evaluate=true, Dialog(tab="Advanced"));

  //Reynold number Re = ( (m_flow / rho / A) * D * rho )  / mu / numParCir.
  Modelica.SIunits.ReynoldsNumber rey=
    m_flow/nParCir/A_pipe*pipeDiaInt/mu_default "Reynolds number";

  // specific mass flow rates
  Real m_flowSp(unit="kg/(m2.s)")=port_a.m_flow/A_floor
    "mass flow rate per unit floor area";
  Real m_flowSpLimit = IDEAS.Utilities.Math.Functions.smoothMax(m_flow_nominal/A_floor/100,m_flowSp, m_flow_nominal/A_floor/100);

  Modelica.SIunits.ThermalInsulance R_t
    "Total equivalent specific resistivity as defined by Koschenz in eqn 4-59";

  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b heatPortEmb
    "Port to the core of a floor heating/concrete activation"
    annotation (Placement(transformation(extent={{-10,90},{10,110}}),
        iconTransformation(extent={{-10,90},{10,110}})));

protected
  final parameter Modelica.SIunits.Length L_r=A_floor/RadSlaCha.T/nParCir
    "Length of one of the parallel circuits";
  final parameter Modelica.SIunits.Area A_pipe=
    Modelica.Constants.pi/4*pipeDiaInt^2
    "Pipe internal cross section surface area";
  final parameter Medium.ThermodynamicState sta_default=
     Medium.setState_pTX(T=Medium.T_default, p=Medium.p_default, X=Medium.X_default);
  final parameter Modelica.SIunits.Density rho_default = Medium.density(sta_default);
  final parameter Modelica.SIunits.DynamicViscosity mu_default = Medium.dynamicViscosity(sta_default)
    "Dynamic viscosity at nominal condition";
  final parameter Modelica.SIunits.SpecificHeatCapacity cp_default = Medium.specificHeatCapacityCp(sta_default)
    "Heat capacity at nominal condition";
  final parameter Modelica.SIunits.MassFlowRate m_flow_nominal_pos = abs(m_flow_nominal)
    "Absolute value of nominal flow rate";
  final parameter Modelica.SIunits.MassFlowRate m_flow_turbulent =  mu_default*pipeDiaInt/4*Modelica.Constants.pi*reyHi
    "Turbulent flow if |m_flow| >= m_flow_turbulent";
  final parameter Modelica.SIunits.Pressure dp_nominal_pos = abs(dp_nominal)
    "Absolute value of nominal pressure";
  final parameter Modelica.SIunits.ReynoldsNumber reyLo=2700
    "Reynolds number where transition to turbulence starts";
  final parameter Modelica.SIunits.ReynoldsNumber reyHi=4000
    "Reynolds number where transition to turbulence ends";

public
  FixedResistances.ParallelFixedResistanceDpM res(
    redeclare package Medium = Medium,
    m_flow_nominal=m_flow_nominal,
    dp_nominal=dp_nominal,
    nParCir=nParCir,
    final use_dh=true,
    dh=pipeDiaInt,
    ReC=reyHi,
    allowFlowReversal=allowFlowReversal,
    from_dp=from_dp,
    homotopyInitialization=homotopyInitialization,
    linearized=linearized,
    dp(nominal=100))
               annotation (Placement(transformation(extent={{20,-10},{40,10}})));
  parameter Boolean homotopyInitialization = true "= true, use homotopy method"
    annotation(Evaluate=true, Dialog(tab="Advanced"));
  parameter Boolean linearized = false
    "= true, use linear relation between m_flow and dp for any flow rate"
    annotation(Evaluate=true, Dialog(tab="Advanced"));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heatFlowSolid
    annotation (Placement(transformation(extent={{-40,70},{-20,90}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow heatFlowWater
    annotation (Placement(transformation(extent={{-40,30},{-20,50}})));
  Sensors.Temperature senTemIn(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-110,20},{-90,40}})));
  Modelica.Blocks.Sources.RealExpression Q(y=
        IDEAS.Utilities.Math.Functions.spliceFunction(
        x=m_flow - m_flow_nominal/100,
        pos=(senTemIn.T - heatPortEmb.T)/R_t*A_floor,
        neg=0,
        deltax=m_flow_nominal/100))
    annotation (Placement(transformation(extent={{-100,50},{-72,70}})));
  Modelica.Blocks.Math.Gain negate(k=-1)
    annotation (Placement(transformation(extent={{-56,36},{-48,44}})));
initial equation
   assert(m_flowMin/A_floor*Medium.specificHeatCapacityCp(sta_default)*(R_w_val_min + R_r_val + R_x_val) >= 0.5,
     "Model is not valid for the set nominal and minimal mass flow rate, discretisation in multiple parts is required");
  if RadSlaCha.tabs then
    assert(RadSlaCha.S_1 > 0.3*RadSlaCha.T, "Thickness of the concrete or screed layer above the tubes is smaller than 0.3 * the tube interdistance. 
    The model is not valid for this case");
    assert(RadSlaCha.S_2 > 0.3*RadSlaCha.T, "Thickness of the concrete or screed layer under the tubes is smaller than 0.3 * the tube interdistance. 
      The model is not valid for this case");
  else
    assert(RadSlaCha.alp2 < 1.212, "In order to use the floor heating model, RadSlaCha.alp2 need to be < 1.212");
    assert(RadSlaCha.d_a/2 < RadSlaCha.S_2, "In order to use the floor heating model, RadSlaCha.alp2RadSlaCha.d_a/2 < RadSlaCha.S_2 needs to be true");
    assert(RadSlaCha.S_1/RadSlaCha.T <0.3, "In order to use the floor heating model, RadSlaCha.S_1/RadSlaCha.T <0.3 needs to be true");
  end if;
equation

  if useSimplifiedRt then
    R_t = IDEAS.Utilities.Math.Functions.inverseXRegularized(2*m_flowSpLimit*cp_default, 1e-8) + R_w_val + R_r_val + R_x_val;
  else
    R_t = IDEAS.Utilities.Math.Functions.inverseXRegularized(m_flowSpLimit*cp_default*(1-exp(-1/((R_w_val+R_r_val+R_x_val+R_c)*m_flowSpLimit*cp_default))), 1e-8)-R_c;
  end if;

  connect(res.port_a, vol.ports[2]) annotation (Line(
      points={{20,0},{-54,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(res.port_b, port_b) annotation (Line(
      points={{40,0},{100,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(heatFlowSolid.port, heatPortEmb) annotation (Line(
      points={{-20,80},{0,80},{0,100}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(heatFlowWater.port, vol.heatPort) annotation (Line(
      points={{-20,40},{-20,10},{-44,10}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(senTemIn.port, port_a) annotation (Line(
      points={{-100,20},{-100,0}},
      color={0,127,255},
      smooth=Smooth.None));
  connect(Q.y, heatFlowSolid.Q_flow) annotation (Line(
      points={{-70.6,60},{-60,60},{-60,80},{-40,80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(negate.y, heatFlowWater.Q_flow) annotation (Line(
      points={{-47.6,40},{-40,40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(negate.u, heatFlowSolid.Q_flow) annotation (Line(
      points={{-56.8,40},{-60,40},{-60,80},{-40,80}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}),
            graphics),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}),
         graphics={Line(
          points={{-90,0},{-80,0},{-80,-60},{-60,-60},{-60,80},{-40,80},{-40,
              -60},{-20,-60},{-20,80},{0,80},{0,-60},{20,-60},{20,80},{40,80},{
              40,-60},{60,-60},{60,80},{80,80},{80,0},{100,0}},
          color={0,0,255},
          smooth=Smooth.None)}),
    Documentation(info="<html>
<p><b>Description</b> </p>
<p>Dynamic model of an embedded pipe for a concrete core activation or a floor heating element. This&nbsp;model&nbsp;is&nbsp;based&nbsp;on&nbsp;the&nbsp;norm&nbsp;prEN&nbsp;15377&nbsp;for&nbsp;the&nbsp;nomenclature&nbsp;but&nbsp;relies&nbsp;more&nbsp;on&nbsp;the&nbsp;background&nbsp;as&nbsp;developed&nbsp;in&nbsp;(Koschenz,&nbsp;2000).&nbsp;The R_x_val for the floor heating is calculated according to the TRNSYS guide lines (TRNSYS, 2007)  There&nbsp;is&nbsp;one&nbsp;major&nbsp;deviation:&nbsp;instead&nbsp;of&nbsp;calculating&nbsp;R_z&nbsp;(to&nbsp;get&nbsp;the&nbsp;mean&nbsp;water&nbsp;temperature&nbsp;in&nbsp;the&nbsp;tube&nbsp;from&nbsp;the&nbsp;supply&nbsp;temperature&nbsp;and&nbsp;flowrate),&nbsp;this&nbsp;mean&nbsp;water&nbsp;temperatue&nbsp;is&nbsp;modelled&nbsp;specifically,&nbsp;based&nbsp;on&nbsp;the&nbsp;mass&nbsp;of&nbsp;the&nbsp;water&nbsp;in&nbsp;the&nbsp;system.<code><font style=\"color: #006400; \">&nbsp;&nbsp;</font></code></p>
<p>The water&nbsp;mass is lumped&nbsp;to&nbsp;TOut.&nbsp;&nbsp;This&nbsp;seems&nbsp;to&nbsp;give&nbsp;the&nbsp;best&nbsp;results (see validation).&nbsp;&nbsp;When&nbsp;lumping&nbsp;to&nbsp;TMean&nbsp;there&nbsp;are&nbsp;additional&nbsp;algebraic&nbsp;constraints&nbsp;to&nbsp;be&nbsp;imposed&nbsp;on&nbsp;TOut&nbsp;which&nbsp;is&nbsp;not so elegant.&nbsp;For most simulations, TOut&nbsp;is&nbsp;more&nbsp;important&nbsp;(influences&nbsp;efficiencies&nbsp;of&nbsp;heat&nbsp;pumps,&nbsp;storage&nbsp;stratification&nbsp;etc.)</p>
<p>This model gives&nbsp;exactly&nbsp;the&nbsp;same&nbsp;results&nbsp;as&nbsp;the&nbsp;norm&nbsp;in&nbsp;both&nbsp;dynamic&nbsp;and&nbsp;static&nbsp;results,&nbsp;but&nbsp;is&nbsp;also&nbsp;able&nbsp;to&nbsp;cope&nbsp;with&nbsp;no-flow&nbsp;conditions.</p>
<p><h4>Assumptions and limitations </h4></p>
<p><ol>
<li>Dynamic model, can be used for situations with variable or no flow rates and sudden changes in supply temperauture&nbsp;</li>
<li>Nomenclature&nbsp;from&nbsp;EN&nbsp;15377.</li>
<li>Water mass is lumped to TOut</li>
</ol></p>
<p><h4>Model use</h4></p>
<p>The embeddedPipe model is to be used together with a <a href=\"modelica://IDEAS.Thermal.Components.Emission.NakedTabs\">NakedTabs</a> in a <a href=\"modelica://IDEAS.Thermal.Components.Emission.Tabs\">Tabs</a> model, or the heatPortEmb can be used to couple it to a structure similar to NakedTabs in a building element. </p>
<p>Configuration&nbsp;of&nbsp;the&nbsp;model&nbsp;can&nbsp;be&nbsp;tricky:&nbsp;the&nbsp;speed&nbsp;of&nbsp;the&nbsp;fluid&nbsp;(flowSpeed)&nbsp;is&nbsp;influencing&nbsp;the&nbsp;convective&nbsp;resistance&nbsp;(R_w_val)&nbsp;and&nbsp;therefore&nbsp;these&nbsp;2&nbsp;configurations&nbsp;are&nbsp;NOT&nbsp;the&nbsp;same:</p>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.&nbsp;10&nbsp;m2&nbsp;of&nbsp;floor&nbsp;with&nbsp;100&nbsp;kg/h&nbsp;flowrate&nbsp;(m_flowSp&nbsp;=&nbsp;10&nbsp;kg/h/m2)</p>
<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.&nbsp;1&nbsp;m2&nbsp;of&nbsp;floor&nbsp;heating&nbsp;with&nbsp;10&nbsp;kg/h&nbsp;flowrate&nbsp;(m_flowSp&nbsp;=&nbsp;10&nbsp;kg/h/m2)</p>
<p><br/>The following parameters have to be set:</p>
<p><ul>
<li>RadSlaCha is a record with all the parameters of the geometry, materials and even number of discretization layers in the nakedTabs model. Attention, this record also specifies the <u>floor surface</u>.</li>
<li>mFlow_min is used to check the validity of the operating conditions.</li>
</ul></p>
<p><br/>The validity range of the model is largely checked by assert() statements. When the mass flow rate is too low, discretization of the model is a solution to obtain models in the validity range again. </p>
<p>About&nbsp;the&nbsp;number&nbsp;of&nbsp;elements&nbsp;in&nbsp;the&nbsp;floor&nbsp;construction&nbsp;(see <a href=\"modelica://IDEAS.Thermal.Components.Emission.NakedTabs\">IDEAS.Thermal.Components.Emission.NakedTabs</a>):&nbsp;this&nbsp;seems&nbsp;to&nbsp;have&nbsp;an&nbsp;important&nbsp;impact&nbsp;on&nbsp;the&nbsp;results.&nbsp;&nbsp;1&nbsp;capacity&nbsp;above&nbsp;and&nbsp;below&nbsp;is&nbsp;clearly&nbsp;not&nbsp;enough.&nbsp;&nbsp;No&nbsp;detailed&nbsp;sensitivity&nbsp;study&nbsp;made,&nbsp;but&nbsp;it&nbsp;seems&nbsp;that&nbsp;3&nbsp;capacities&nbsp;on&nbsp;each&nbsp;side&nbsp;were&nbsp;needed&nbsp;in&nbsp;my&nbsp;tests&nbsp;to&nbsp;get&nbsp;good&nbsp;results.</p>
<p><h4>Validation </h4></p>
<p>Validation&nbsp;of&nbsp;the&nbsp;model&nbsp;is&nbsp;not&nbsp;evident&nbsp;with&nbsp;the&nbsp;data&nbsp;in&nbsp;(Koschenz,&nbsp;2000):</p>
<p><ul>
<li>4.5.1&nbsp;is&nbsp;very&nbsp;strange:&nbsp;the&nbsp;results&nbsp;seem&nbsp;to&nbsp;be&nbsp;obtained&nbsp;with&nbsp;1m2&nbsp;and&nbsp;12&nbsp;kg/h&nbsp;total&nbsp;flowrate,&nbsp;but&nbsp;this&nbsp;leads&nbsp;to&nbsp;very&nbsp;low&nbsp;flowSpeed&nbsp;value&nbsp;(although&nbsp;Reynolds&nbsp;number&nbsp;is&nbsp;still&nbsp;high)&nbsp;and&nbsp;an&nbsp;alpha&nbsp;convection&nbsp;of&nbsp;only&nbsp;144&nbsp;W/m2K&nbsp;==&GT;&nbsp;I&nbsp;exclude&nbsp;this&nbsp;case&nbsp;explicitly&nbsp;with&nbsp;an&nbsp;assert&nbsp;statement&nbsp;on&nbsp;the&nbsp;flowSpeed</li>
<li>4.6&nbsp;is&nbsp;ok&nbsp;and&nbsp;I&nbsp;get&nbsp;exactly&nbsp;the&nbsp;same&nbsp;results,&nbsp;but&nbsp;this&nbsp;leads&nbsp;to&nbsp;extremely&nbsp;low&nbsp;supply&nbsp;temperatures&nbsp;in&nbsp;order&nbsp;to&nbsp;reach&nbsp;20&nbsp;W/m2</li>
<li>4.5.2&nbsp;not&nbsp;tested</li>
</ul></p>
<p><br/>A specific report of this validation can be found in IDEAS/Specifications/Thermal/ValidationEmbeddedPipeModels_20111006.pdf</p>
<p><h4>Examples</h4></p>
<p>See combinations with NakedTabs in a Tabs model. </p>
<p><h4>References</h4></p>
<p>[Koshenz, 2000] - Koschenz, Markus, and Beat Lehmann. 2000. <i>Thermoaktive Bauteilsysteme - Tabs</i>. D&uuml;bendorf: EMPA D&uuml;bendorf. </p>
<p>[TRNSYS, 2007] - Multizone Building modeling with Type 56 and TRNBuild.</p>
</html>", revisions="<html>
<p><ul>
<li>2014 March, Filip Jorissen: IDEAS baseclasses</li>
<li>2013 May, Roel De Coninck: documentation</li>
<li>2012 April, Roel De Coninck: rebasing on common Partial_Emission</li>
<li>2011, Roel De Coninck: first version and validation</li>
</ul></p>
</html>"));
end EmbeddedPipe;
