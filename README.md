# YAGFES
*YAGFES*, o *Yet Another Groundwater Flow Equation Solver*, es un paquete desarrollado en Python 3 que asiste al usuario con la elaboración de modelos computacionales de flujo de agua en medios porosos bidimensionales para el cálculo de distribuciones de *carga hidráulica*, *velocidades de Darcy* y de *conos de abatimiento*. *YAGFES* se desarrollo como una herramienta libre, de código abierto y gratuito. Al tratarse de un paquete desarrollado en Python 3, *YAGFES* puede ser empleado en diferentes sistemas operativos sin mayores complicaciones.

*YAGFES* está diseñado para resolver la ecuación de flujo en un medio poroso bidimensional con las siguientes consideraciones:
1. Se estudia el flujo en un acuífero confinado.
2. El acuífero es verticalmente homogéneo y cada sección vertical del acuífero se encuentra en equilibrio hidrostático.

Para elaborar un modelo el usuario debe proporcionar tres clases de archivos principales: archivos de configuración (*.init*), archivos de malla (*.msh*/*.fdmsh*) y archivos de pozos (*.well*). Las funciones incluidas en los módulos de elemento y diferencias finitas asistirán al usuario en la creación de archivos de lectura auxiliares (por ejemplo: los archivos *.bc* que contienen las condiciones de frontera).

El paquete está compuesto por tres módulos:
1. *fem_m.py*: El módulo que contiene las funciones requeridas para dar solución a la ecuación de flujo empleando el método de elemento finito.
2. *fdm_m.py*: El módulo que contiene las funciones requeridas para dar solución a la ecuación de flujo empleando el método de diferencias finitas.
3. *aux_f.py*: El módulo que contiene funciones adicionales para visualizar resultados y para asistir a lectura de archivos.
