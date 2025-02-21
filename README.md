[![DOI](https://zenodo.org/badge/905216080.svg)](https://doi.org/10.5281/zenodo.14905733)
# RIF spatial transcriptomics with Nanostring GeoMX
![image](https://github.com/user-attachments/assets/29c381a9-8f0a-4d41-9c35-5e6cc5aa64e9)



## Running the analysis
### With singularity

1. Install [`Singularity`](https://docs.sylabs.io/guides/3.8/user-guide/)

2. Build the singularity container:
    ```console
     sudo singularity build runAnalysis.img Singularity
    ```
3. Run the analysis and render the html notebooks with a single command:
    ```console
     ./runAnalysis.img
    ```

### Alternatively using RScript

1.	Install the needed R packages
    ```console
     RScript install/install.R
    ```
2.	Run the analysis and render the html notebooks
    ```console
     RScript runAnalysis.R
    ```
