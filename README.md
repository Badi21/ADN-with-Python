# Análisis de Secuencias de ADN en Python

Este programa en Python se encarga de analizar secuencias de ADN que se encuentran en archivos con formato FASTA. Utiliza expresiones regulares para realizar distintas operaciones, como la generación de proteínas a partir de las secuencias de ADN, y ofrece una funcionalidad para buscar genes específicos utilizando expresiones regulares en los nombres de los genes.

## Funcionalidades Principales

1. **Interpretación de Archivos FASTA:**
    - El programa interpreta archivos que siguen el formato FASTA, los cuales contienen información sobre varias secuencias de ADN.
2. **Generación de Proteínas:**
    - Transforma las secuencias de ADN en proteínas, aplicando el código genético estándar.
3. **Búsqueda mediante Expresiones Regulares:**
    - Facilita la búsqueda de genes específicos a través de expresiones regulares aplicadas a los nombres de los genes.
4. **Generación de Resultados:**
    - Produce un archivo de texto que contiene las proteínas generadas a partir de las secuencias de ADN.

## Uso del Programa

1. **Ejecución del Programa:**
    - Ejecuta el programa proporcionando como entrada el nombre del archivo FASTA que deseas analizar. El programa procesará las secuencias y generará un archivo de salida con las proteínas resultantes.

    ```bash
    python nombre_del_programa.py
    ```

2. **Búsqueda mediante Expresiones Regulares:**
    - Después de la ejecución, el programa te permitirá ingresar expresiones regulares para buscar genes específicos dentro del conjunto de datos.
3. **Resultados:**
    - El programa generará un archivo de texto (`Protein.txt`) que contendrá los nombres de los genes y las secuencias de proteínas correspondientes.

## Ejemplo de Uso

1. **Ejecución del Programa:**
    - Selecciona el archivo FASTA (`nombre_del_archivo_DNA.txt`) cuando se te solicite.
2. **Búsqueda mediante Expresiones Regulares:**
    - Tras la ejecución, introduce expresiones regulares para buscar genes específicos.
3. **Resultados:**
    - Revisa el archivo `Protein.txt` para obtener las proteínas generadas y la información asociada a los genes.

## Notas Importantes

- Asegúrate de que el archivo FASTA tenga un formato correcto para garantizar una interpretación precisa.
- La búsqueda mediante expresiones regulares en los nombres de los genes distingue entre mayúsculas y minúsculas.
