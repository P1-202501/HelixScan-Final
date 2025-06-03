# 🧬Proyecto Épico Científico: Simulador de Procesos Biológicos Moleculares🧬
## 👫🏻Laura Mariela Parra y Samuel Gallego

## 🔎¡Explorando los Secretos del ADN y las Proteínas para el Entendimiento Científico!

✅¡Bienvenido/a al Proyecto Épico Científico! Este programa es una herramienta integral diseñada para simular y analizar procesos fundamentales de la biología molecular: la transcripción de ADN a ARN y la traducción de ARN a proteínas. Nuestro objetivo es crear un "universo épico" basado en la realidad científica, integrando no solo las transformaciones biológicas, sino también potentes capacidades de análisis de datos y monitoreo profesional para ofrecer una visión completa y robusta de las secuencias genéticas.

---

## 🎯Objetivo General del Proyecto

Finalizar el desarrollo del universo épico basado en un problema científico real, integrando nuevas funcionalidades programadas, herramientas de monitoreo profesional en la nube, análisis de datos estructurados, y buenas prácticas de colaboración con Git y PEP8.

---

## ✨Funcionalidades Principales

El simulador ofrece las siguientes características clave:

* **Transformación de Secuencias Biológicas:**
  
  **🧫Limpieza de ADN:** Prepara las secuencias de ADN eliminando formatos innecesarios (ej. `5'-` y `-3'`) y normalizándolas a mayúsculas.
  
  **🧪Validación de ADN:** Asegura que la secuencia de ADN solo contenga las bases nitrogenadas válidas (A, T, C, G).
  
  **👩🏻‍🔬Transcripción (ADN a ARN):** Convierte una secuencia de ADN en su correspondiente secuencia de ARN, reemplazando las Timinas (T) por Uracilos (U).
  
  **🧬Traducción (ARN a Proteína):** Traduce la secuencia de ARN en una cadena de aminoácidos, utilizando la tabla de codones estándar. La traducción se detiene al encontrar un codón "STOP".
  
  **🔬Cadena Complementaria de ADN (¡NUEVO!):** Genera la cadena de ADN complementaria (A se empareja con T, y C con G), una función esencial en la replicación y reparación del ADN.

* **Detección y Monitoreo de Anomalías Biológicas:**
  
  **🦠Codones Desconocidos:** Alerta si se encuentran secuencias de tres nucleótidos (codones) que no corresponden a ningún aminoácido conocido en la tabla genética, indicando posibles errores o mutaciones.
  
  **🤒Longitud Anormal de Proteínas:** Monitorea y advierte si las proteínas resultantes son inusualmente cortas o largas, lo cual podría implicar una traducción incompleta o errónea.
  
  **☑️Reporte a Sentry:** Todas las advertencias y errores críticos son enviados a Sentry para un monitoreo profesional en la nube.

* **Análisis de Datos con Archivos y Pandas:**
  
  **🗳️Guardado de Resultados:** Cada secuencia procesada se guarda automáticamente en un archivo `datos_adn.csv`, acumulando un registro histórico de todas las operaciones.
      
  **📊Análisis Estadístico:** Permite cargar y analizar el archivo `datos_adn.csv` para obtener estadísticas descriptivas (promedio, mínimo, máximo, desviación estándar) sobre la longitud de las proteínas, filtrar secuencias por criterios específicos (ej. proteínas largas), y contar secuencias únicas de ADN limpio.
      
  **📉Análisis de Frecuencia de Aminoácidos (¡NUEVO!):** Calcula y muestra el porcentaje de cada aminoácido presente en una secuencia de proteína dada, ofreciendo una caracterización bioquímica útil.

* **Monitoreo Profesional en la Nube (Sentry):**
  
  **🔒Registro Centralizado:** Los errores críticos, advertencias y anomalías biológicas detectadas son enviadas a Sentry, una plataforma de monitoreo en la nube, para una gestión y depuración profesional.
      
  **💻Información Detallada:** Sentry captura el contexto completo de los eventos, facilitando la identificación y resolución de problemas en el código o anomalías en los datos biológicos.
      
  **🕹️Inclusión de excepción:** El proyecto incluye al menos una excepción capturada.

* **Seguridad y Buenas Prácticas:**
  
  **🔑Variables de Entorno (`.env`):** Las claves sensibles, como el DSN de Sentry, se gestionan de forma segura utilizando un archivo `.env` y no se exponen directamente en el código ni en el repositorio público.
      
  **👾Control de Versiones (`GitFlow`):** El equipo ha trabajado siguiendo la metodología GitFlow para la gestión de ramas (`main`, `develop`, `feature`, `hotfix`) y el ciclo de desarrollo, asegurando una colaboración estructurada y un historial de commits claro.
      
  **💊Modularidad con Clases:** El código está organizado en al menos tres clases definidas por el equipo (`RegistradorSecuencias`, `GestorADN`, `AnalizadorProteinas`), cada una con atributos y métodos propios, y utilizada en el programa principal, mejorando la legibilidad, mantenibilidad y reusabilidad.

---

## 🛠️Cómo Empezar

Sigue estos pasos para poner en marcha el simulador en tu entorno local.

### 📋Requisitos Previos

Asegúrate de tener instalado:

* [Python 3.8+](https://www.python.org/downloads/)
* `pip` (gestor de paquetes de Python)

### 📦Instalación de Dependencias

   **1️⃣Clona el Repositorio:**
    ```bash
    git clone <URL_DE_TU_REPOSITORIO>
    cd <nombre-de-tu-carpeta-proyecto> # ej. cd HELISCAN-PROYECTO-FINAL-SG-LP
    ```
   
   **2️⃣Crea el entorno virtual (recomendado):**
    ```bash
    python -m venv venv
    # En Windows
    .\venv\Scripts\activate
    # En macOS/Linux
    source venv/bin/activate
    ```
   
   **3️⃣Instala las dependencias:**
    ```bash
    pip install -r requirements.txt
    ```

### 🔑Configuración de Sentry (`.env`)

Para activar el monitoreo profesional, necesitas configurar tu DSN de Sentry:

**1️⃣Crea un archivo `.env`** en la misma carpeta que `main.py`.

**2️⃣Añade tu DSN de Sentry** en este archivo, de la siguiente manera:
    ```
    SENTRY_DSN="[https://ejemplo@sentry.io/123456](https://ejemplo@sentry.io/123456)"
    ```
    *(Reemplaza el valor con tu DSN real obtenido de tu proyecto en Sentry.)*

**3️⃣Verifica tu `.gitignore`:** Es fundamental que `.env` esté incluido en tu archivo `.gitignore` para prevenir que tu clave privada sea subida al repositorio público.

---

## 🚀Ejecución del Programa

Para iniciar el simulador, simplemente ejecuta el script `main.py` desde la raíz de tu proyecto:

🎈El programa te presentará un menú interactivo:

```bash
python main.py

¡Bienvenido al Procesador de Secuencias Biológicas!

--- Menú Principal ---
1. Procesar una nueva secuencia de ADN
2. Analizar datos de secuencias (desde 'datos_adn.csv')
3. Obtener cadena de ADN complementaria
4. Analizar frecuencia de aminoácidos de una proteína
5. Salir
Seleccione una opción (1, 2, 3, 4 o 5):
```

*😀Opción 1:* Procesar una nueva secuencia de ADN: Te pedirá que ingreses una secuencia de ADN (ej. 5'-ATCGTGC-3'). El programa la limpiará, validará, transcribirá a ARN, traducirá a proteína, y reportará cualquier anomalía. Los resultados se guardarán en datos_adn.csv.

```bash
--- Resultados del Procesamiento ---
ADN ingresado:
ADN limpio: 
ARN: 
Proteína: 
Longitud de Proteína: 
------------------------------------
```

*😋Opción 2:* Analizar datos de secuencias: Cargará y analizará los datos de todas las secuencias previamente procesadas y guardadas en datos_adn.csv, mostrando estadísticas clave de longitud de proteínas y búsquedas específicas.

```bash
--- Análisis de Datos de Secuencias ---
Total de secuencias procesadas: #

Estadísticas de longitud de proteínas:
  Promedio: # aminoácidos
  Mínimo: # aminoácidos
  Máximo: # aminoácidos
  Desviación estándar: # aminoácidos

--- Búsqueda de Proteínas Largas (Longitud > 10) ---
Si/No hay proteínas con longitud mayor a 10.

Número de secuencias de ADN limpias únicas: x
--- Fin del Análisis ---
```

*🤗Opción 3:* Obtener cadena de ADN complementaria: Te pedirá una secuencia de ADN y generará su cadena complementaria.

```bash
--- Cadena Complementaria ---
ADN original: 
ADN complementario:
------------------------------
```

*☺️Opción 4:* Analizar frecuencia de aminoácidos de una proteína: Te pedirá una secuencia de proteína y mostrará la frecuencia porcentual de cada aminoácido.

```bash
Ingrese la secuencia de proteína para analizar su frecuencia (ej. FLSYWSTOP): ej: MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNELTEFAKTCVADESAENCDKSLHTLFGDKLCT
VATLRETYGEMADCCAKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQK

--- Frecuencia de Aminoácidos ---
  A: 28 (11.43%)
  L: 27 (11.02%)
  K: 26 (10.61%)
  E: 21 (8.57%)
  F: 19 (7.76%)
  D: 16 (6.53%)
  S: 13 (5.31%)
  R: 13 (5.31%)
  C: 12 (4.90%)
  T: 10 (4.08%)
  P: 9 (3.67%)
  V: 8 (3.27%)
  Y: 8 (3.27%)
  Q: 8 (3.27%)
  G: 7 (2.86%)
  H: 6 (2.45%)
  N: 5 (2.04%)
  I: 4 (1.63%)
  W: 3 (1.22%)
  M: 2 (0.82%)
```

*😉Opción 5:* Salir: Terminará el programa.

```bash
Gracias por usar el Procesador de Secuencias Biológicas. ¡Hasta pronto!
```

## ☁️Monitoreo en Sentry

Este proyecto se integra con Sentry para el monitoreo de errores. Cualquier excepción no manejada o errores específicos capturados por SequenceLogger serán reportados a tu panel de control de Sentry en tiempo real.

![Imagen de WhatsApp 2025-06-01 a las 21 20 19_b63b3bb2](https://github.com/user-attachments/assets/c17da92c-afa1-4704-88b9-0e1f744f662f)



## 📄Archivo de Registro (programa.log)
Todos los eventos internos del programa (información, advertencias, errores, etc.) se registran en el archivo programa.log ubicado en la misma carpeta raíz del proyecto. Este archivo es útil para la depuración y para revisar el flujo de ejecución del programa de forma detallada.

## 🤝Colaboración y Estilo

🐱GitFlow: El equipo ha trabajado siguiendo la metodología GitFlow para la gestión de ramas (main, develop, feature, hotfix) y el ciclo de desarrollo, asegurando una colaboración estructurada y un historial de commits claro.

🦜PEP8: El código se adhiere a las directrices de estilo PEP8 para garantizar la legibilidad y la consistencia en todo el proyecto.

## 💡Próximos Pasos / Ideas para el Futuro

👀Visualizaciones de Datos: Integrar librerías como matplotlib o seaborn para crear gráficos dinámicos (ej. histogramas de longitud de proteínas, gráficos de barras de frecuencia de aminoácidos) directamente desde la opción de análisis de datos.

🐝Simulación de Mutaciones: Añadir una funcionalidad para simular mutaciones en las secuencias y observar su impacto en las proteínas resultantes.
Interfaz Gráfica de Usuario (GUI): Desarrollar una interfaz gráfica de usuario para una experiencia aún más intuitiva y accesible para usuarios no técnicos.
