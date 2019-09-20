# YAGFES
<em>YAGFES</em>, o <em>Yet Another Groundwater Flow Equation Solver</em>, es un paquete desarrollado en Python 3 que asiste al usuario con la elaboración de modelos computacionales de flujo de agua en medios porosos bidimensionales para el cálculo de distribuciones de <em>carga hidráulica</em>, <em>velocidades de Darcy</em> y de <em>conos de abatimiento</em> con la presencia de pozos. <em>YAGFES</em> se desarrollo como una herramienta libre, de código abierto y gratuito. Al tratarse de un paquete desarrollado en Python 3, <em>YAGFES</em> puede ser empleado en diferentes sistemas operativos sin mayores complicaciones.

<em>YAGFES</em> está diseñado para resolver la ecuación de flujo en un medio poroso bidimensional con las siguientes consideraciones:
<ol>
<li> Se estudia el flujo en un acuífero confinado.</li>
<li> El acuífero es verticalmente homogéneo y cada sección vertical del acuífero se encuentra en equilibrio hidrostático.</li>
</ol>

Para elaborar un modelo el usuario debe proporcionar tres clases de archivos principales: archivos de configuración (<em>*.init</em>), archivos de malla (<em>*.msh</em>/<em>*.fdmsh</em>) y archivos de pozos (<em>*.well</em>). Las funciones incluidas en los módulos de elemento y diferencias finitas asistirán al usuario en la creación de archivos de lectura auxiliares (por ejemplo: los archivos <em>*.bc</em> que contienen las condiciones de frontera).

El paquete está compuesto por tres módulos:
<ol>
<li><em>fem_m.py</em>: El módulo que contiene las funciones requeridas para dar solución a la ecuación de flujo empleando el método de elemento finito.</li>
<li><em>fdm_m.py</em>: El módulo que contiene las funciones requeridas para dar solución a la ecuación de flujo empleando el método de diferencias finitas.</li>
<li><em>aux_f.py</em>: El módulo que contiene funciones adicionales para visualizar resultados y para asistir a lectura de archivos.</li>
</ol>
