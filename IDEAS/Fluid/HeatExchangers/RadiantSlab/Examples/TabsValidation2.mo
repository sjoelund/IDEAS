within IDEAS.Fluid.HeatExchangers.RadiantSlab.Examples;
model TabsValidation2
  extends TabsValidation(
    Tsurface(table={{0,293.1311914},{3412.889973,293.4716074},{6967.518295,
          293.8872577},{11749.64694,294.3599756},{17882.52664,294.8710168},{
          20464.24445,294.7971306},{26615.61176,294.4053601},{33011.93992,
          294.0513352},{39406.34228,293.7913532},{50838.61856,293.5151945},{
          61899.98716,293.3516947},{71976.89049,293.2629157}}),
    doubleRamp(interval=5*3600),
    Tcore(table={{0,293.1500642},{3510.335088,294.7130376},{7061.111824,
          295.3167737},{10737.83541,295.7701053},{17974.19438,296.3945757},{
          21685.58223,295.1551354},{25134.29195,294.7431442},{32270.89485,
          294.2390358},{42109.77019,293.7739569},{57475.28566,293.4246181},{
          71977.27565,293.2441071}}));
  annotation (__Dymola_Commands(file=
          "modelica://IDEAS/Resources/Scripts/Dymola/Fluid/HeatExchangers/RadiantSlab/Examples/TabsValidation2.mos"
        "Simulate and plot"), Documentation(info="<html>
<p>Validation of concrete core activation (TABS) model based on Koschenz.</p>
</html>", revisions="<html>
<ul>
<li>
January 30, 2015 by Filip Jorissen:<br/>
First implementation.
</li>
</ul>
</html>"));
end TabsValidation2;
