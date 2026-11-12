###################################
# **Proyecto glucuronidasas**
#Este proyecto tiene como objetivo presentar datos preliminares para la solicitud de plan nacional.
###################################

**Requisitos**
### Para replicar el análisis se necesitan los siguientes recursos:

- Descargarse las secuencias de NCBI estudio (https://doi.org/10.1002/1873-3468.13527)
- Obtener el catalogo de GUS ()
- Modules: FastQC, Trimmomatic, Bowtie2


```bash
module load Miniconda
conda activate mi_proyecto

```

## Estructura del Pipline
1. Descargar archivos brutos. 
2. Limpieza de las secuencias. (01_QC_host_removal_SE.sbatch|02_ResumenQC.sbatch)
3. Preparación del catálogo.
