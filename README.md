# Simulacion-modelo-XY-Transicion-de-K-T
Código utilizado para producir los resultados del trabajo sobre la transición de Kosterlitz-Thouless en el modelo XY bidimensional

El código se compone esencialmente de 3 archivos.
- **maincorr.cpp**: Es el archivo principal que simula el sistema.
    - **Compilado**:  Para compilar usar el comando: `g++ -O3 -fopenmp -std=c++17 -o simulationcorr maincorr.cpp` 
    - **Ejecución**:  Para ejecutar usar: `./simulationcorr <L> <T_min> <T_max> <Pasos de Temperatura> <Repositorio> <Bins de termalización> <Bins de medida> <Modo>` . Para medir la correlación usar `<Modo>=1`. Para medir solo configuraciones cualquier otro número.
    - **Outputs**:    La ejecución de este código producirá una carpeta de nombre <Repositorio> en la que se guardarán los siguientes archivos:
        - La última configuración de cada temperatura en archivos separados *configatT=<temperatura>.data*,
        - La función de correlación medida en todos los pasos de medición para cada temperatura en *corratT=<temperatura>.data*,
        - La energía de la configuración medida en todos los pasos de medición para cada temperatura en *EnergyatT=<temperatura>*,
        - La lista ordenada de temperaturas en *variables.data*
- **Media correlación.ipynb**: Jupyter notebook en el que se calcula la correlación media de todos pasos de Monte Carlo y se guardan en archivos de nombre *mean_corratT=<temperatura>.data*
- **Media energías.ipynb**: Otro notebook en el que se calcula la energía media por sitio de red de todos los pasos de Monte Carlo para todas las temperaturas y se guardan en archivos de nombre *mean_energy_<directorio>.data*
- **Gráficos.ipynb**: Jupyter notebook para producir los plots mostrados en el trabajo
